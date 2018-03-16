#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
pipe.py
2017-03-22: finished basically
2017-06-26: runsecond model refined
2017-06-27: modified runsecond model, as extend the second time result
2017-06-30: modified the combining method
2017-07-20: modifed some output.
2017-07-28: modified overlaps, change to the samll one.
2017-08-03: re-deisgn the underlying data structure, improve the speed as 2 fold and reduce the memory as 1/2
2017-08-07: modified pipeline, as accecpt multiple eps.
2017-08-10: to do, re-design data structure as using HDF5 as botoom structure to reduce memory usage
2017-08-18: changed combine method for loops when runing a series of eps.
2018-03-13: modified combining method for multiple eps
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#python library
import os, time, sys, shutil, gzip, copy
from glob import glob
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

#cLoops
from cLoops.settings import *
from cLoops.utils import getLogger, callSys, cFlush, mainHelp
from cLoops.io import parseRawBedpe, txt2jd, parseJd, loops2washU, loops2juice
from cLoops.cDBSCAN import cDBSCAN as DBSCAN
from cLoops.ests import estFragSize, estIntSelCutFrag
from cLoops.cPlots import plotFragSize, plotIntSelCutFrag
from cLoops.cModel import getIntSig, markIntSig, markIntSigHic

#global settings
global logger


def singleDBSCAN(f, eps, minPts, cut=0):
    """
    Run DBSCAN to detect interactions for one chromosome.
    #mat is list, every is [ pointId,x,y ]
    """
    dataI, readI, dataS, readS, dis, dss = [], [], [], [], [], []
    key, mat = parseJd(f, cut)
    if len(mat) == 0:
        return key, f, dataI, dataS, list(dis), list(dss)
    #data for interaction records, read for readId
    report = "Clustering %s and %s using eps as %s, minPts as %s,pre-set distance cutoff as > %s" % (
        key[0], key[1], eps, minPts, cut)
    logger.info(report)
    db = DBSCAN(mat, eps, minPts)
    labels = pd.Series(db.labels)
    mat = np.array(mat)
    mat = pd.DataFrame(
        mat[:, 1:].astype("float"), index=mat[:, 0], columns=["X", "Y"])
    nlabels = set(labels.values)
    #collect clusters
    for label in nlabels:
        los = list(labels[labels == label].index)
        sub = mat.loc[los, :]
        #BEDPE format,+1 to escape the error that exact the same start and end
        #2017-05-18, changed to remove such interactions
        if int(np.min(sub["X"])) == int(np.max(sub["X"])) or int(
                np.min(sub["Y"])) == int(np.max(sub["Y"])):
            continue
        r = [
            key[0],
            int(np.min(sub["X"])),
            int(np.max(sub["X"])),
            key[1],
            int(np.min(sub["Y"])),
            int(np.max(sub["Y"])),
            #sub.shape[0],
            #",".join(map(str, los)),
            #los
        ]
        if r[2] < r[4]:
            dataI.append(r)
            readI.extend(los)
        else:
            dataS.append(r)
            readS.extend(los)
    report = "Clustering %s and %s finished. Estimated %s self-ligation reads and %s inter-ligation reads" % (
        key[0], key[1], len(readS), len(readI))
    logger.info(report)
    if len(dataI) > 0:
        dis = mat.loc[readI, "Y"] - mat.loc[readI, "X"]
    if len(dataS) > 0:
        dss = mat.loc[readS, "Y"] - mat.loc[readS, "X"]
    return key, f, dataI, dataS, list(dis), list(dss)


def runDBSCAN(fs, eps, minPts, cut=0, cpu=1):
    """
    Run DBSCAN to detect interactions for all chromosomes.
    """
    ds = Parallel(n_jobs=cpu)(delayed(singleDBSCAN)(f, eps, minPts, cut)
                              for f in fs)
    dataI, dataS, dis, dss = {}, [], [], []
    for d in ds:
        if len(d[2]) == 0:
            continue
        dataI[d[0]] = {"f": d[1], "records": d[2]}
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
    return False


def combineTwice(dataI, dataI_2):
    """
    Combining multiple clustering result.
    """
    for key in dataI_2.keys():
        if key not in dataI:
            dataI[key] = {
                "f": dataI_2[key]["f"],
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


def runStat(dataI, minPts, cut, cpu, fout, hichip=0):
    """
    Calling p-values of interactions for all chromosomes.
    """
    logger.info(
        "Starting estimate significance for interactions using distance cutoff as %s"
        % cut)
    ds = Parallel(n_jobs=cpu)(delayed(getIntSig)(
        dataI[key]["f"], dataI[key]["records"], minPts, cut)
                              for key in dataI.keys())
    ds = [d for d in ds if d is not None]
    if len(ds) == 0:
        logger.error("Something wrong, no loops found, sorry, bye.")
        return 1
    ds = pd.concat(ds)
    try:
        if hichip:
            ds = markIntSigHic(ds)
        else:
            ds = markIntSig(ds)
        ds.to_csv(fout + ".loop", sep="\t", index_label="loopId")
    except:
        logger.warning(
            "Something wrong happend to significance estimation, only output called loops"
        )
        ds.to_csv(fout + "_raw.loop", sep="\t", index_label="loopId")
    return 0


def pipe(fs,
         fout,
         eps,
         minPts,
         chroms="",
         cpu=1,
         tmp=0,
         hic=0,
         washU=0,
         juice=0,
         cut=0,
         org="hg38"):
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))
    #1.pre-processing data as merge replicates, seperate into cis and trans PETs
    #cis and trans PETs files
    if os.path.isdir(fout):
        r = "working directory %s exists, return." % fout
        logger.error(r)
        return
    os.mkdir(fout)
    #2.read in PETs, collect cis PETs and convert to .jd format, if not assigned eps (eps=0),estimate eps
    cfs, ds = parseRawBedpe(fs, fout, chroms, cut, logger)
    cfs = Parallel(n_jobs=cpu)(delayed(txt2jd)(f) for f in cfs)
    #3. run the DBSCAN only one time
    if eps == 0:
        frags = estFragSize(ds)
        eps = frags * 2
        dataI, dataS, dis, dss = runDBSCAN(cfs, eps, minPts, cut, cpu)
        #4.estimate cutoff of self-ligation and inter-ligation
        cut, frags = estIntSelCutFrag(np.array(dis), np.array(dss))
        plotIntSelCutFrag(dis, dss, cut, frags, prefix=fout + "_disCutoff")
        logger.info(
            "Estimated inter-ligation and self-ligation distance cutoff as %s"
            % cut)
    else:
        #3.1. run DBSCAN for multiple times
        dataI = {}
        cuts = []
        for ep in eps:
            dataI_2, dataS_2, dis_2, dss_2 = runDBSCAN(cfs, ep, minPts, cut,
                                                       cpu)
            if len(dataI_2) == 0 or len(dataS_2) == 0:
                logger.info(
                    "ERROR: no inter-ligation PETs or self-ligation PETs detected for eps %s,can't model the distance cutoff,continue anyway"
                    % ep)
                continue
            cut_2, frags = estIntSelCutFrag(np.array(dis_2), np.array(dss_2))
            plotIntSelCutFrag(
                dis_2,
                dss_2,
                cut_2,
                frags,
                prefix=fout + "_eps%s_disCutoff" % ep)
            logger.info(
                "Estimated inter-ligation and self-ligation distance cutoff as %s for eps=%s."
                % (cut_2, ep))
            cuts.append(cut_2)
            #combine the first round and second round result
            dataI = combineTwice(dataI, dataI_2)
        #cut = min(cuts)
        cut = np.max(cuts)
        #cut = np.median(cuts)
    #5.estimate the significance
    e = runStat(dataI, minPts, cut, cpu, fout, hic)
    if e:
        shutil.rmtree(fout)
        return
    #6.convert loops and cis-PETs into washU visulization format or juicebox
    if washU:
        loops2washU(fout + ".loop", fout + "_loops_washU.txt", logger)
    if juice:
        loops2juice(fout + ".loop", fout + "_loops_juicebox.txt", logger)
    #7.remove temple files
    if tmp == False:
        shutil.rmtree(fout)


def main():
    start = datetime.now()
    global logger
    fn = os.path.join(os.getcwd(), "cLoops.log")
    logger = getLogger(fn)
    op = mainHelp()
    report = "Command line: cLoops -f {} -o {} -m {} -eps {} -minPts {} -p {} -w {} -j {} -s {} -c {} -hic {} -cut {}".format(
        op.fnIn, op.fnOut, op.mode, op.eps, op.minPts, op.cpu, op.washU,
        op.juice, op.tmp, op.chroms, op.hic, op.cut)
    logger.info(report)
    if op.mode == 0:
        if "," in str(op.eps):
            eps = map(int, op.eps.split(","))
        else:
            eps = int(op.eps)
            if eps != 0:
                eps = [eps]
        if op.minPts == 0:
            logger.error("minPts not assigned!")
            return
        else:
            minPts = op.minPts
        hic = op.hic
    if op.mode == 1:
        eps = [500, 1000, 2000]
        minPts = 5
        hic = 0
    if op.mode == 2:
        eps = [1000, 2000, 5000]
        minPts = 5
        hic = 0
    if op.mode == 3:
        eps = [1000, 2000, 4000, 6000, 8000, 10000]
        minPts = 20
        hic = 1
    report = "mode:%s\t eps:%s\t minPts:%s\t hic:%s\t" % (op.mode, eps, minPts,
                                                          hic)
    logger.info(report)
    pipe(
        op.fnIn.split(","), op.fnOut, eps, minPts, op.chroms, op.cpu, op.tmp,
        hic, op.washU, op.juice, op.cut)
    usedtime = datetime.now() - start
    r = "cLoops finished. Used CPU time: %s Bye!\n\n\n" % usedtime
    logger.info(r)


if __name__ == "__main__":
    main()