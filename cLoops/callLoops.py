#!/usr/bin/env python2.7
#--coding:utf-8 --

#python library
import os, time, sys, shutil, gzip, copy
from glob import glob
from datetime import datetime
from logging import getLogger
from fractions import Fraction

#3rd library
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#cLoops
from cLoops.settings import *
from cLoops.io import readPets, isIntraPetClass
from cLoops.cDBSCAN import cDBSCAN as DBSCAN
from cLoops.ests import estFragSize, estIntSelCutFrag
from cLoops.cPlots import plotIntSelCutFrag


def singleDBSCAN(petfile,
                 petclass,
                 eps,
                 minPts,
                 logfile_name,
                 cut=0,
                 strict_intra=True,
                 anchor_ratio=1):
                 anchor_ratio=1):
    """
    Run DBSCAN to detect interactions for one chromosome.
    mat is list, every item is [ pointId,x,y ]
    @para anchor_ratio: Loop anchors width ratio for assymmetric data
    @para anchor_ratio: Loop anchors width ratio for assymmetric data
    """
    logger = getLogger(logfile_name)
    dataI, readI, dataS, readS, dis, dss = [], [], [], [], [], []
    #put off distance filter
    key, mat, intra_flag = readPets(
        petfile, petclass, cut=0, strict_intra=strict_intra)
    #self-ligation could only happen in intra pet class
    if intra_flag and cut > 0:
        d = abs(mat[:, 2] - mat[:, 1])
        dss.extend(list(d[d < cut]))
        mat = mat[d >= cut, :]
    if len(mat) == 0:
        return key, petclass, dataI, dataS, list(dis), list(dss)
    #data for interaction records, read for readId
    report = "Clustering %s using eps as %s, minPts as %s, pre-set distance cutoff as > %s" % (
        petclass, eps, minPts, cut)
    logger.info(report)
    anchor_ratio = Fraction(anchor_ratio)
    anchor_ratio = Fraction(anchor_ratio)
    coor_factor = [1, 1]
    if anchor_ratio > 1:
    if anchor_ratio > 1:
        anchor_ratio = anchor_ratio.limit_denominator(1000)
        anchor_ratio = anchor_ratio.limit_denominator(1000)
        coor_factor = [anchor_ratio.denominator, anchor_ratio.numerator]
        coor_factor = [anchor_ratio.denominator, anchor_ratio.numerator]
        mat[:, 1:3] = mat[:, 1:3] * coor_factor
        eps = eps * max(coor_factor)
    elif anchor_ratio < 1:
    elif anchor_ratio < 1:
        anchor_ratio = (1 / anchor_ratio).limit_denominator(1000)
        anchor_ratio = (1 / anchor_ratio).limit_denominator(1000)
        coor_factor = [anchor_ratio.numerator, anchor_ratio.denominator]
        coor_factor = [anchor_ratio.numerator, anchor_ratio.denominator]
        mat[:, 1:3] = mat[:, 1:3] * coor_factor
        eps = eps * max(coor_factor)
    db = DBSCAN(mat, eps, minPts)
    labels = pd.Series(db.labels)
    #mat = np.array(mat)
    mat = pd.DataFrame(
        mat[:, 1:].astype("float"), index=mat[:, 0], columns=["X", "Y"])
    nlabels = set(labels.values)
    #collect clusters
    for label in nlabels:
        los = list(labels[labels == label].index)
        sub = mat.loc[los, :]
        #BEDPE format, +1 to escape the error that exact the same start and end
        #2017-05-18, changed to remove such interactions
        if int(np.min(sub["X"])) == int(np.max(sub["X"])) or int(
                np.min(sub["Y"])) == int(np.max(sub["Y"])):
            continue
        r = [
            key[0],
            int(np.min(sub["X"]) / coor_factor[0]),
            int(np.max(sub["X"]) / coor_factor[0]),
            key[1],
            int(np.min(sub["Y"]) / coor_factor[1]),
            int(np.max(sub["Y"]) / coor_factor[1]),
            #sub.shape[0],
            #",".join(map(str, los)),
            #los
        ]
        r += key[2:4]
        if not intra_flag or r[2] < r[4] or r[1] > r[5]:
            dataI.append(r)
            readI.extend(los)
        else:
            dataS.append(r)
            readS.extend(los)
    report = "Clustering %s finished. Estimated %s self-ligation reads and %s inter-ligation reads" % (
        petclass, len(readS), len(readI))
    logger.info(report)
    if intra_flag:
        if len(dataI) > 0:
            dis = abs(mat.loc[readI, "Y"] - mat.loc[readI, "X"])
        if len(dataS) > 0:
            dss.extend(list(abs(mat.loc[readS, "Y"] - mat.loc[readS, "X"])))
    return key, petclass, dataI, dataS, list(dis), list(dss)


def runDBSCAN(petfile,
              total_petclass,
              eps,
              minPts,
              logfile_name,
              cut=0,
              cpu=1,
              strict_intra=True,
              anchor_ratio=1):
              anchor_ratio=1):
    """
    Run DBSCAN to detect interactions for all chromosomes.
    """
    ds = Parallel(n_jobs=cpu)(delayed(
        singleDBSCAN)(petfile, petclass, eps, minPts, logfile_name, cut,
                      strict_intra, anchor_ratio) for petclass in total_petclass)
                      strict_intra, anchor_ratio) for petclass in total_petclass)
    dataI, dataS, dis, dss = {}, [], [], []
    for d in ds:
        if len(d[2]) == 0:
            continue
        dataI[d[0]] = {"petclass": d[1], "records": d[2]}
        dataS.extend(d[3])
        dis.extend(d[4])
        dss.extend(d[5])
    return dataI, dataS, dis, dss


def checkSameLoop(ra, rb):
    """
    check if two anchors are exact same
    """
    if ra[1] == rb[1] and ra[2] == rb[2] and ra[4] == rb[4] and ra[5] == rb[5]:
        return True
    else:
        return False


def combineTwice(dataI, dataI_2):
    """
    Combine multiple clustering result.
    """
    for key in dataI_2.keys():
        if key not in dataI:
            dataI[key] = {
                "petclass": dataI_2[key]["petclass"],
                "records": dataI_2[key]["records"]
            }
        else:
            for ra in dataI_2[key]["records"]:
                flag = 1
                for rb in dataI[key]["records"]:
                    if checkSameLoop(ra, rb):
                        flag = 0
                        break
                if flag:
                    dataI[key]["records"].append(ra)
    return dataI


def filterClusterByDis(data, cut, strict_intra=True):
    """
    Filter cis inter-ligation clusters by distances
    """
    for key in data:
        if not isIntraPetClass(key, strict_intra):
            continue
        nr = []
        for r in data[key]["records"]:
            d = (r[4] + r[5]) / 2 - (r[1] + r[2]) / 2
            if d >= cut:
                nr.append(r)
        data[key]["records"] = nr
    return data


def callLoops(petfile,
              fout,
              total_petclass,
              eps,
              minPts,
              cut=0,
              strict_intra=True,
              symmetrical=True,
              anchor_ratio=1,
              anchor_ratio=1,
              logfile_name='',
              plot=0,
              cpu=1,
              ds=[]):
    logger = getLogger(logfile_name)
    #3. run the DBSCAN
    if eps == 0:
        #3.0. run the DBSCAN only one time
        frags = estFragSize(ds)
        eps = [frags * 2]
    #3.1. run DBSCAN for multiple times
    dataI = {}
    dis, dss = [], []
    cuts = [
        cut,
    ]
    for ep in eps:
        #multiple minPts, for minPts like 50,40,30,20, the bigger minPts always get smaller distance cutoff
        for m in minPts:
            dataI_2, dataS_2, dis_2, dss_2 = runDBSCAN(
                petfile, total_petclass, ep, m, logfile_name, cut, cpu,
                strict_intra, anchor_ratio)
                strict_intra, anchor_ratio)
            if len(dataI_2) == 0:
                logger.info(
                    "ERROR: no inter-ligation PETs or self-ligation PETs detected for eps %s minPts %s, can't model the distance cutoff,continue anyway"
                    % (ep, m))
                continue
            cut_2, frags = estIntSelCutFrag(np.array(dis_2), np.array(dss_2))
            if plot and cut_2 != 0 and frags != 0:
                plotIntSelCutFrag(
                    dis_2,
                    dss_2,
                    cut_2,
                    frags,
                    prefix=fout + "_eps%s_minPts%s_disCutoff" % (ep, m))
            logger.info(
                "Estimated inter-ligation and self-ligation distance cutoff as %s for eps=%s, minPts=%s."
                % (cut_2, ep, m))
            #experimental
            cuts.append(cut_2)
            #cut = max(cuts)
            cut = cut_2
            dataI_2 = filterClusterByDis(dataI_2, cut, strict_intra)
            dataI = combineTwice(dataI, dataI_2)
    cuts = [c for c in cuts if c > 0]
    if len(cuts) > 0:
        cut = np.min(cuts)
    else:
        cut = 0
    #cut = np.max(cuts)
    return dataI, cut
