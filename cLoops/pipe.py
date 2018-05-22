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
2018-03-18: modified pre-processing bedpe files.
2018-03-27: modifed merging clustering by remove close PETs
2018-03-30: multiple minPts mode added
2018-03-31: modified self-ligation and inter-ligation distance cutoff selection.
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
from cLoops.io import parseRawBedpe, loops2washU, loops2juice
from cLoops.cModel import getIntSig, markIntSig, markIntSigHic
from cLoops.callLoops import callLoops

#global settings
global logger


def runStat(petfile, dataI, minPts, cut, strict_intra, cpu, fout,
            hichip=False):
    """
    Calling p-values of interactions for all chromosomes.
    """
    logger.info(
        "Starting estimate significance for interactions using distance cutoff as %s"
        % cut)
    ds = Parallel(n_jobs=cpu)(
        delayed(getIntSig)(petfile, dataI[key]["petclass"], dataI[key][
            "records"], minPts, cut, strict_intra) for key in dataI.keys())
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
         hic=False,
         washU=0,
         juice=0,
         cut=0,
         plot=0,
         strict_intra=True,
         anchor_ratio=1,
         symmetrical=True,
         cis_only=True,
         use_strand=False,
         logfile_name=''):
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))
    #1.pre-processing data as merge replicates, seperate into cis and trans PETs
    #cis and trans PETs files
    petfile = fout + ".pet"
    if os.path.isfile(petfile):
        r = "Temperory file %s exists, return." % petfile
        logger.error(r)
        return
    #2.read in PETs, collect cis PETs and convert to .jd format, if not assigned eps (eps=0),estimate eps
    if eps == 0:
        cfs, ds = parseRawBedpe(
            fs,
            petfile,
            chroms,
            cut,
            strict_intra,
            logger,
            symmetrical,
            cis_only,
            use_strand,
            rmRep=True,
            cal_ds=True)
    else:
        cfs, ds = parseRawBedpe(
            fs,
            petfile,
            chroms,
            cut,
            strict_intra,
            logger,
            symmetrical,
            cis_only,
            use_strand,
            rmRep=False,
            cal_ds=False)
    #3. call loops
    dataI, cut = callLoops(petfile, fout, cfs, eps, minPts, cut, strict_intra,
                           symmetrical, anchor_ratio, logfile_name, plot, cpu,
                           ds)

    #4.estimate the significance
    e = runStat(petfile, dataI, minPts, cut, strict_intra, cpu, fout, hic)
    if e:
        os.remove(fout)
        return
    #5.convert loops and cis-PETs into washU visulization format or juicebox
    if washU:
        loops2washU(fout + ".loop", fout, logger)
    if juice:
        loops2juice(fout + ".loop", fout, logger)
    #6.remove temple files
    if tmp == False:
        os.remove(petfile)


def main():
    start = datetime.now()
    global logger
    fn = os.path.join(os.getcwd(), "cLoops.log")
    logger = getLogger(fn)
    op = mainHelp()
    report = "Command line: cLoops -f {} -o {} -m {} -eps {} -minPts {} -p {} -w {} -j {} -s {} -c {} -cut {} -plot {}".format(
        op.fnIn, op.fnOut, op.mode, op.eps, op.minPts, op.cpu, op.washU,
        op.juice, op.tmp, op.chroms, op.cut, op.plot)
    if op.use_strand != '00':
        report += " -usestrand %s" % op.use_strand
    if op.hic:
        report += " -hic"
    if not op.strict_intra_flag:
        report += " -notStrictIntra"
    if not op.symmetrical_flag:
        report += " -assymmetric"
    if not op.cis_only_flag:
        report += " -includeTrans"
    logger.info(report)
    if op.use_strand == '00':
        op.use_strand = [False, False]
    elif op.use_strand == '10':
        op.use_strand = [True, False]
    elif op.use_strand == '01':
        op.use_strand = [False, True]
    elif op.use_strand == '11':
        op.use_strand = [True, True]
    if op.mode == 0:
        #parse eps
        if "," in str(op.eps):
            eps = map(int, op.eps.split(","))
            eps.sort(reverse=False)
        else:
            eps = int(op.eps)
            if eps != 0:
                eps = [eps]
        #parse minPts
        if "," in str(op.minPts):
            minPts = map(int, op.minPts.split(","))
            minPts.sort(reverse=True)
        else:
            minPts = int(op.minPts)
            if minPts != 0:
                minPts = [minPts]
            else:
                logger.error("minPts not assigned!")
                return
        hic = op.hic
    if op.mode == 1:
        eps = [500, 1000, 2000]
        minPts = [5]
        hic = False
        op.anchor_ratio = 1
    if op.mode == 2:
        eps = [1000, 2000, 5000]
        minPts = [5]
        hic = False
        op.anchor_ratio = 1
    if op.mode == 3:
        eps = [5000, 7500, 10000]
        minPts = [50, 40, 30, 20]
        hic = True
        op.anchor_ratio = 1
    if op.mode == 4:
        eps = [2500, 5000, 7500, 10000]
        minPts = [30, 20]
        hic = True
        op.anchor_ratio = 1
    if op.mode == 5:
        eps = [1500]
        minPts = [5]
        hic = False
        op.symmetrical_flag = False
        op.cis_only_flag = False
        op.use_strand = [True, False]
        op.anchor_ratio = 0.05
    report = "mode:%s\t eps:%s\t minPts:%s\t hic:%s\t" % (op.mode, eps, minPts,
                                                          hic)
    if op.mode in [1, 2, 3, 4]:
        op.symmetrical_flag = True
        op.cis_only_flag = True
        op.use_strand = [False, False]
    else:
        if op.use_strand[0] or op.use_strand[1]:
            report += " use_strand:%d%d" % (op.use_strand[0], op.use_strand[1])
        if not op.symmetrical_flag:
            report += " assymmetric"
        if not op.cis_only_flag:
            report += " includeTrans"
    logger.info(report)
    pipe(
        op.fnIn.split(","), op.fnOut, eps, minPts, op.chroms, op.cpu, op.tmp,
        hic, op.washU, op.juice, op.cut, op.plot, op.strict_intra_flag,
        op.anchor_ratio, op.symmetrical_flag, op.cis_only_flag, op.use_strand,
        fn)
    usedtime = datetime.now() - start
    r = "cLoops finished. Used CPU time: %s Bye!\n\n\n" % usedtime
    logger.info(r)


if __name__ == "__main__":
    main()
