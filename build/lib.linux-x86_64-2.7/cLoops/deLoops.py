#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
deLoops
Calling differentially enriched loops between one TF for different cells/conditions or similar situation. 
2017-06-16: basically finished, poisson test is good, pseduo add 1 for both treatment and control.
2017-06-19: small bugs for np.mean([]) will get a nan, fixed
2017-06-21: improved efficiency.
2017-08-03: modified data structure
2018-03-22: to be continued by CHEN Zhaoxiong
"""
__date__ = "2017-06-16"
__modified__ = ""
__email__ = "caoyaqiang@picb.ac.cn chenxingwei@picb.ac.cn"

#general library
import os

#3rd library
import numpy as np
import pandas as pd
from scipy.stats import poisson
from joblib import Parallel, delayed

#cLoops
from cLoops.utils import getLogger, cFlush, deloopHelp
from cLoops.io import readPet, parseIv
from cLoops.cModel import getGenomeCoverage, getCounts, getNearbyPairRegions, getPETsforRegions, getBonPvalues

#global settings
global logger


def preDs(f, h5file, chroms=[], ivac=6, ivbc=7):
    """
    Prepare input datasets, f is the .loop file and h5file is .pet file contain parsed PETs.
    """
    records = {}  #organized by chromosomes.
    if len(chroms) > 0:
        for c in chroms:
            records[c] = {"rs": {}, "f": ""}
    for i, line in enumerate(open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        #only using significant loops
        if float(line[-1]) < 1:
            continue
        iva = parseIv(line[ivac])
        ivb = parseIv(line[ivbc])
        if len(chroms) > 0 and iva[0] not in chroms and ivb[0] not in chroms:
            continue
        if iva[0] not in records:
            records[iva[0]] = {"rs": {}, "f": ""}
        records[iva[0]]["rs"][line[0]] = iva + ivb
    for chrom in records.keys():
        if len(records[chrom]["rs"]) == 0:
            del records[chrom]
            continue
        f = os.path.join(d, "%s-%s.jd" % (chrom, chrom))
        if os.path.isfile(f):
            records[chrom]["f"] = f
        else:
            logger.warning(
                "%s not found, however there are loops in that chromosome." %
                f)
            del records[chrom]
    return records


def getPermutatedBg(ivas, ivbs, model):
    rabs = []
    for na in ivas:
        try:
            nra = set(np.abs(list(getCounts(na, model))))
        except:
            continue
        nralen = float(len(nra))
        if nralen == 0:
            continue
        for nb in ivbs:
            try:
                nrb = set(np.abs(list(getCounts(nb, model))))
            except:
                continue
            if len(nrb) == 0:
                continue
            nrab = len(nra.intersection(nrb))
            #collect the value for poisson test
            rabs.append(nrab)
    #to escape warning RuntimeWarning: Mean of empty slice, this will leads to nan
    if len(rabs) == 0:
        mrabs = 0.0
    else:
        mrabs = float(np.mean(rabs))
    return mrabs


def estSigOneLoop(iva, ivb, modelt, modelc, normratio, win=5):
    """
    Estimating the significance for a loop.
    """
    rat, rbt, rabt = getPETsforRegions(iva, ivb, modelt)
    rac, rbc, rabc = getPETsforRegions(iva, ivb, modelc)
    ivas, ivbs = getNearbyPairRegions(iva, ivb, win=win)
    mrabt = getPermutatedBg(ivas, ivbs, modelt)
    mrabc = getPermutatedBg(ivas, ivbs, modelc)
    #in case for all 0 p-values, this make the p-values higer-estimated.
    lams = (np.array([mrabc, rabc]) + 1.0) * normratio
    lam = np.max(lams)
    pop = poisson.sf(rabt - 1.0, lam)
    fc = rabt / lam
    pop = max([pop, 1e-300])
    return pop, fc


def estSigTvsC(rs, modelt, Nt, modelc, Nc, pre):
    """
    Estimation of significance for treat vs control
    """
    normratio = float(Nt) / float(Nc)
    ds = {}
    i = 0
    for key, r in rs.items():
        i += 1
        if i % 100 == 0:
            report = "Estimating %s loops for %s" % (i, pre)
            cFlush(report)
        chrom = r[0]
        iva = [r[1], r[2]]
        ivb = [r[4], r[5]]
        p, fc = estSigOneLoop(iva, ivb, modelt, modelc, normratio)
        niva = "%s:%s-%s" % (chrom, iva[0], iva[1])
        nivb = "%s:%s-%s" % (chrom, ivb[0], ivb[1])
        ds[key] = {
            "iva": niva,
            "ivb": nivb,
            "poisson_p-value": p,
            "FoldEnrichment": fc
        }
    print
    if len(ds) == 0:
        return None
    ds = pd.DataFrame(ds).T
    ds["poisson_p-value_corrected"] = getBonPvalues(ds["poisson_p-value"])
    return ds


def estSigOneChr(rst, bedpet, rsc, bedpec, pre, dis=0):
    """
    Estimating the significances for the loops in one chromosome.
    """
    #all variables with suffix t is treatment, with suffix c in control
    logger.info("Building genomic coverage model for %s" % bedpet)
    modelt, Nt = getGenomeCoverage(bedpet, dis)
    logger.info("Building genomic coverage model for %s" % bedpec)
    modelc, Nc = getGenomeCoverage(bedpec, dis)
    logger.info("Starting estimation for %s of treatment." % pre)
    dst = estSigTvsC(rst, modelt, Nt, modelc, Nc, pre)
    logger.info("Estimation for %s of treatment finished." % pre)
    logger.info("Starting estimation for %s of control." % pre)
    dsc = estSigTvsC(rsc, modelc, Nc, modelt, Nt, pre)
    logger.info("Estimation for %s of control finished." % pre)
    return dst, dsc


def callDeLoops(ra, rb, prea, preb, dis=0, cpu=1):
    logger.info("Calling differentially enriched loops for %s vs %s" % (prea,
                                                                        preb))
    ds = Parallel(n_jobs=cpu)(delayed(estSigOneChr)(ra[key]["rs"], ra[key][
        "f"], rb[key]["rs"], rb[key]["f"], key, dis) for key in ra.keys())
    dsa, dsb = [], []
    for t in ds:
        if t[0] is not None:
            dsa.append(t[0])
        if t[1] is not None:
            dsb.append(t[1])
    dsa, dsb = pd.concat(dsa), pd.concat(dsb)
    dsa.to_csv(prea + ".deloop", sep="\t", index_label="loopId")
    dsb.to_csv(preb + ".deloop", sep="\t", index_label="loopId")


def main():
    global logger
    fn = os.path.join(os.getcwd(), "deLoops.log")
    logger = getLogger(fn)
    op = deloopHelp()
    if op.chroms == "":
        chroms = []
    else:
        chroms = set(op.chroms.split(","))
    cfa, cfb, dfa, dfb = op.fa, op.fb, op.da, op.db
    ra = preDs(cfa, dfa, chroms)
    rb = preDs(cfb, dfb, chroms)
    prea = os.path.split(dfa)[1]
    preb = os.path.split(dfb)[1]
    keys = set(ra.keys()).intersection(set(rb.keys()))
    for key in ra.keys():
        if key not in keys:
            del ra[key]
            logger.info("No match of %s in %s or %s" % (key, cfa, dfa))
    for key in rb.keys():
        if key not in keys:
            del rb[key]
            logger.info("No match of %s in %s or %s" % (key, cfa, dfa))
    callDeLoops(ra, rb, prea, preb, op.dis, op.cpu)


main()
