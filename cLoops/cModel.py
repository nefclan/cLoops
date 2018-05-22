#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
Stastical significance is tested for every chromosome using the local permutated background.
2018-02-01: improved data structure for genomecoverage,much faster and less memory than previouse version for significance calling,slightly changed the loops boundary.
2018-03-08: modified ChIA-PET significant loops cutoff
2018-03-16: key change, sliding step changed to the half of the mean anchor size always get more significant loops
"""
__date__ = "2017-03-15"
__modified__ = ""
__email__ = "caoyaqiang@picb.ac.cn chenxingwei@picb.ac.cn"

#general library
import gc

#3rd library
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, binom, poisson, combine_pvalues

#cLoops
from cLoops.io import readPets, parseIv, isIntraPetClass, getPetFilePara
from cLoops.utils import cFlush


def getCorLink(cs):
    """
    @param cs: [1,2,3,4], a list for the coordinates x or y
    @rtype: dic, keys is the coordinate, value is the closest next coordinate and points index in this coordinate
    """
    ts = {}
    for i, c in enumerate(cs):
        if c not in ts:
            ts[c] = []
        ts[c].append(i)
    keys = sorted(ts.keys())
    for i in xrange(len(keys) - 1):
        ni = keys[i]
        nj = keys[i + 1]
        ts[ni] = {"next": nj, "points": ts[ni]}
    ts[nj] = {"next": None, "points": ts[nj]}
    return ts


def getGenomeCoverage(petfile, petclass, cut=0, strict_intra=True):
    """
    Build the genomic model for random access. Could use a lot of memory.
    @param petfile: .pet file
    @param petclass: pet class to analysis
    @param cut: distance cutoff for self-ligation PETs.
    """
    key, mat, intra_flag = readPets(petfile, petclass, cut, strict_intra)
    j = mat.shape[0]
    if j < 2:
        return None, 0
    xs = getCorLink(mat[:, 1])
    ys = getCorLink(mat[:, 2])
    return [xs, ys], j, intra_flag


def getCounts(iv, ts):
    ps = []
    pos = None
    for i in xrange(iv[0], iv[1]):
        if i in ts:
            pos = i
            break
    while pos <= iv[1] and pos != None:
        ps.extend(ts[pos]["points"])
        pos = ts[pos]["next"]
    return set(ps)


def getPETsforRegions(iva, ivb, model, triangular_flag):
    raSource = getCounts(iva, model[0])
    rbTarget = getCounts(ivb, model[1])
    if triangular_flag:
        raTarget = getCounts(iva, model[1])
        ra = len(raSource.union(raTarget))
        rbSource = getCounts(ivb, model[0])
        rb = len(rbSource.union(rbTarget))
    else:
        ra = len(raSource)
        rb = len(rbTarget)
    rab = len(raSource.intersection(rbTarget))
    return ra, rb, rab


def getNearbyPairRegions(iva, ivb, win=5):
    """
    @param iva: [start,end] 
    Get the nearby regions for interacting two locus, win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    """
    ivas, ivbs = [], []
    ca = sum(iva) / 2
    cb = sum(ivb) / 2
    sa = (iva[1] - iva[0]) / 2
    sb = (ivb[1] - ivb[0]) / 2
    step = (sa + sb) / 2
    for i in xrange(0 - win, win + 1):
        if i == 0:
            continue
        niva = [iva[0], iva[1]]
        niva[0] = max([0, ca + i * step - sa])
        niva[1] = max([0, ca + i * step + sa])
        nivb = [ivb[0], ivb[1]]
        nivb[0] = max([0, cb + i * step - sb])
        nivb[1] = max([0, cb + i * step + sb])
        ivas.append(niva)
        ivbs.append(nivb)
    return ivas, ivbs


def getMultiplePsFdr(iva, ivb, model, N, win=5, triangular_flag=True):
    """
    for the interval a and b, searching its nearby windows to estimate FDR and p-values.
    return ra, rb, rab, es, fdr, hyp, pop, nbp
    """
    ra, rb, rab = getPETsforRegions(iva, ivb, model, triangular_flag)
    hyp = max([1e-300, hypergeom.sf(rab - 1.0, N, ra, rb)])
    ivas, ivbs = getNearbyPairRegions(iva, ivb, win=win)
    #nras is a list for storing points ids for permutated regions
    nras, nrbs = [], []
    for na in ivas:
        nraSource = getCounts(na, model[0])
        if triangular_flag:
            nraTarget = getCounts(na, model[1])
            nra = nraSource.union(nraTarget)
        else:
            nra = nraSource
        nras.append(nra)
    for nb in ivbs:
        nrbTarget = getCounts(nb, model[1])
        if triangular_flag:
            nrbSource = getCounts(nb, model[0])
            nrb = nrbSource.union(nrbTarget)
        else:
            nrb = nrbTarget
        nrbs.append(nrb)
    #caculating the permutated background
    rabs, nbps = [], []
    for nra in nras:
        nralen = float(len(nra))
        for nrb in nrbs:
            nrblen = len(nrb)
            nrab = float(len(nra.intersection(nrb)))
            if nrab > 0:
                #collect the value for poisson test
                rabs.append(nrab)
                #collect the possibility for following binomial test
                den = nrab / (nralen * nrblen)
                nbps.append(den)
            else:
                nbps.append(0.0)
                rabs.append(0.0)
    if len(rabs) == 0:
        return ra, rb, rab, np.inf, 0.0, hyp, 0.0, 1e-300, 1e-300,
    rabs = np.array(rabs)
    #local fdr
    fdr = len(rabs[rabs > rab]) / float(len(rabs))
    mrabs = float(np.mean(rabs))
    #enrichment score
    if mrabs > 0:
        es = rab / np.mean(rabs[rabs > 0])
    else:
        es = np.inf
    #simple possion test
    lam = mrabs
    pop = max([1e-300, poisson.sf(rab - 1.0, lam)])
    #simple binomial test
    bp = np.mean(nbps) * ra * rb / N
    nbp = max([1e-300, binom.sf(rab - 1.0, N - rab, bp)])
    return ra, rb, rab, es, fdr, hyp, pop, nbp



def getBonPvalues(ps):
    """
    Return the Bonferroni corrected p-values.
    """
    ps = np.array(ps)
    ps = ps * len(ps)
    ps[ps > 1.0] = 1.0
    return ps


def checkOneEndOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if ya > xb or xa > yb:
        return False
    return True


def checkOverlap(ivai, ivbi, strandi, ivaj, ivbj, strandj):
    """
    check the overlap of two anchors,ra=[chr_left,left_start,left_end,chr_right,right_start,right_end, strand]
    """
    if strandi != strandj:
        return False
    if checkOneEndOverlap(ivai[1], ivai[2], ivaj[1], ivaj[2]) and checkOneEndOverlap(
            ivbi[1], ivbi[2], ivbj[1], ivbj[2]):
        return True
    return False


def removeDup(ds, bpcut=1e-5):
    """
    Remove overlapped called loops, keep the more significant one for multiple eps result. 
    @param:ds, from getIntSig
    @param:bpcut, bionomial p-value cutoff
    """
    uniqueds = {}
    rekeys = set()
    loopscore = {}
    retable = {}
    for key in ds.keys():
        #if ds[key]["binomial_p-value"] > bpcut:
        #    continue
        #loopscore[key] = ds[key]["binomial_p-value"]
        loopscore[key] = float(ds[key]["rab"]) / ds[key]["ra"] / ds[key]["rb"]
    #significant loop keys
    keys = sorted(loopscore, key=loopscore.get, reverse=True)
    for i in xrange(len(keys)):
        keyi = keys[i]
        if keyi in rekeys:
            continue
        ivai = parseIv(ds[keyi]["iva"])
        ivbi = parseIv(ds[keyi]["ivb"])
        if "strand" in ds[keyi]:
            strandi = ds[keyi]["strand"]
        else:
            strandi = ''
        reds = retable.pop(keyi, [])
        for j in xrange(i + 1, len(keys)):
            keyj = keys[j]
            if keyj in rekeys:
                continue
            ivaj = parseIv(ds[keyj]["iva"])
            ivbj = parseIv(ds[keyj]["ivb"])
            if "strand" in ds[keyj]:
                strandj = ds[keyj]["strand"]
            else:
                strandj = ''
            if checkOverlap(ivai, ivbi, strandi,
                ivaj, ivbj, strandj):
                reds.append(keyj)
        if len(reds) == 0:
            uniqueds[keyi] = ds[keyi]
        elif ds[keyi]["binomial_p-value"] <= bpcut:
            uniqueds[keyi] = ds[keyi]
            for keyj in reds:
                rekeys.add(keyj)
                retable.pop(keyj, [])
        else:
            for keyj in reds:
                if keyj in retable:
                    retable[keyj].append(keyi)
                else:
                    retable[keyj] = [keyi]
    return uniqueds

def removeDup_bak(ds, bpcut=1e-5):
    """
    Remove overlapped called loops, keep the more significant one for multiple eps result. 
    @param:ds, from getIntSig
    @param:bpcut, bionomial p-value cutoff
    """
    uniqueds = {}
    reds = {}
    rekeys = set()
    keys = ds.keys()
    for i in xrange(len(keys)):
        keyi = keys[i]
        if keyi in rekeys:
            continue
        ivai = parseIv(ds[keyi]["iva"])
        ivbi = parseIv(ds[keyi]["ivb"])
        if "strand" in ds[keyi]:
            strandi = ds[keyi]["strand"]
        else:
            strandi = ''
        #1 means unique loops
        flag = 1
        #collect overlapped loops
        for j in xrange(i + 1, len(keys)):
            keyj = keys[j]
            if keyj in rekeys:
                continue
            ivaj = parseIv(ds[keyj]["iva"])
            ivbj = parseIv(ds[keyj]["ivb"])
            if "strand" in ds[keyj]:
                strandj = ds[keyj]["strand"]
            else:
                strandj = ''
            flagj = checkOverlap(ivai, ivbi, strandi,
                ivaj, ivbj, strandj)
            #there is overlapped loops,collect them
            if flagj:
                if keyi not in reds:
                    reds[keyi] = [keyi]
                    rekeys.add(keyi)
                reds[keyi].append(keyj)
                rekeys.add(keyj)
                flag = 0
        #collect unique loops
        if flag:
            uniqueds[keyi] = ds[keyi]
    #for overlapped loops, choose the more significant ones
    for key in reds.keys():
        ts = {}
        for t in reds[key]:
            if ds[t]["binomial_p-value"] > bpcut:
                continue
            #ts[t] = ds[t]["binomial_p-value"]
            #first select the significant loops, then select the loops with  smaller anchors and higher density
            ts[t] = float(ds[t]["rab"]) / ds[t]["ra"] / ds[t]["rb"]
            """
            Used for debugging
                iva = parseIv(ds[t]["iva"])
                ivb = parseIv(ds[t]["ivb"])
                a = iva[2]-iva[1]+ivb[2]-ivb[1]
                b = float(ds[t]["rab"])/ds[t]["ra"]/ds[t]["rb"]
                c = float(ds[t]["rab"])/a
                print t
                print a,b,c,ds[t]["rab"],ds[t]["ra"],ds[t]["rb"],ds[t]["ES"],ds[t]["binomial_p-value"],ds[t]["poisson_p-value"]
            print 
            """
        if len(ts) == 0:
            continue
        ts = pd.Series(ts)
        ts.sort_values(inplace=True, ascending=False)
        uniqueds[ts.index[0]] = ds[ts.index[0]]
    return uniqueds



def getIntSig(petfile, petclass, records, minPts, discut, strict_intra=True):
    """
    @param:discut, distance cutoff determined for self-ligation pets.
    """
    print "Starting estimate significance for %s candidate interactions in %s" % (len(records), petclass)
    model, N, intra_flag = getGenomeCoverage(petfile, petclass, discut, strict_intra)
    print "Genomic coverage model built from %s" % petclass
    if N == 0:
        print "No PETs parsed as requiring distance cutoff >%s from %s" % (
            discut, petclass)
        return None
    ds = {}
    i = 0
    petclass_key = petclass.split("_")
    para = getPetFilePara(petfile)
    #determine whether contact matrix is a triangular matrix
    triangular_flag = para['symmetrical'] and isIntraPetClass(petclass_key, True)
    for r in records:
        chroms = [r[0], r[3]]
        key = "%s_%s" % (petclass, i)
        #key = "%s_%s" % (petclass, i)
        iva = [max(0, r[1]), r[2]]
        ivb = [max(0, r[4]), r[5]]
        #filter loops
        distance = 0
        if intra_flag:
            distance = abs(sum(ivb) / 2.0 - sum(iva) / 2.0)
            if distance < discut:
                continue
        ra, rb, rab = getPETsforRegions(iva, ivb, model, triangular_flag)
        if rab < max(minPts):
            continue
        i += 1
        if i % 100 == 0:
            cFlush("%s interaction p-values estimated for %s" % (i, petclass))
        ra, rb, rab, es, fdr, hyp, pop, nbp = getMultiplePsFdr(
            iva, ivb, model, N, 5, triangular_flag)
        #this part should be furthur modified, as for most ideable data, there are no noise, so the es should be inf, however, not possible
        ds[key] = {
            "distance": distance,
            "ra": ra,
            "rb": rb,
            "rab": rab,
            "ES": es,
            "FDR": fdr,
            "hypergeometric_p-value": hyp,
            "poisson_p-value": pop,
            "binomial_p-value": nbp,
            "iva": "%s:%s-%s" % (chroms[0], iva[0], iva[1]),
            "ivb": "%s:%s-%s" % (chroms[1], ivb[0], ivb[1]),
        }
        strandinfo = '_'.join(petclass_key[2:4])
        if strandinfo != '':
            ds[key]["strand"] = strandinfo

    #memory usage
    del model
    gc.collect()
    print
    if len(ds.keys()) == 0:
        return None
    ds = removeDup(ds)
    if len(ds.keys()) == 0:
        return None
    ds = pd.DataFrame(ds).T
    ds["poisson_p-value_corrected"] = getBonPvalues(ds["poisson_p-value"])
    ds["binomial_p-value_corrected"] = getBonPvalues(ds["binomial_p-value"])
    ds["hypergeometric_p-value_corrected"] = getBonPvalues(
        ds["hypergeometric_p-value"])
    return ds


def markIntSig(ds,
               escut=2.0,
               fdrcut=1e-2,
               bpcut=1e-5,
               ppcut=1e-5,
               hypcut=1e-10):
    """
    gpcut is general p-value cutoff for binomial test, poisson test and hypergeometric test.
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller hypergeometric result
    c = ds.loc[b.index, "hypergeometric_p-value"]
    c = c[c <= hypcut]
    #smaller poisson and binomial
    d = ds.loc[c.index, "poisson_p-value"]
    d = d[d <= ppcut]
    e = ds.loc[d.index, "binomial_p-value"]
    e = e[e <= bpcut]
    rs = e.index
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[rs] = 1.0
    ds["significant"] = ns
    return ds


def markIntSigHic(ds, escut=2.0, fdrcut=0.01, bpcut=1e-5, ppcut=1e-5):
    """
    For HiChIP/HiC data, hypergeometric test is not working, poisson and binomial works well. For mouse data, pcut=1e-3 maybe better
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller poisson and binomial result
    c = ds.loc[b.index, "poisson_p-value"]
    c = c[c <= ppcut]
    d = ds.loc[b.index, "binomial_p-value"]
    d = d[d <= bpcut]
    e = c.index.intersection(d.index)
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[e] = 1.0
    ds["significant"] = ns
    return ds



def markStripSig(ds, escut=2.0, fdrcut=0.1, ppcut=1e-5, es_cut=0.2):
    """
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller poisson
    c = ds.loc[b.index, "poisson_p-value"]
    c = c[c <= ppcut]
    #higher enrichment of es_a, es_b
    d = ds.loc[c.index, "ES_ra"]
    d = d[d >= es_cut]
    e = ds.loc[c.index, "ES_rb"]
    e = e[e >= es_cut]
    f = d.index.union(e.index)
    #higher PETs
    #rs = (ds["ES"] >= escut) & (ds["FDR"] <= fdrcut) & (ds["poisson_p-value"] <= ppcut) & (
    #    (ds["ES_ra"] >= es_cut) | (ds["ES_rb"] >= es_cut))
    rs = f
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[rs] = 1.0
    ds["significant"] = ns
    return ds

def getNearbyPairRegionsForStrips(iva, ivb, win=5):
    """
    @param iva: [start,end] 
    """
    lena = iva[1] - iva[0]
    lenb = ivb[1] - ivb[0]
    ivas, ivbs = [], []
    ca = sum(iva) / 2
    cb = sum(ivb) / 2
    sa = (iva[1] - iva[0]) / 2
    sb = (ivb[1] - ivb[0]) / 2
    if lena > lenb:
        step = sb
        for i in xrange(0 - win, win + 1):
            if i == 0:
                continue
            nivb = [ivb[0], ivb[1]]
            nivb[0] = max([0, cb + i * step - sb])
            nivb[1] = max([0, cb + i * step + sb])
            ivas.append(iva)
            ivbs.append(nivb)
        return ivas, ivbs
    if lena < lenb:
        step = sa
        for i in xrange(0 - win, win + 1):
            if i == 0:
                continue
            niva = [iva[0], iva[1]]
            niva[0] = max([0, ca + i * step - sa])
            niva[1] = max([0, ca + i * step + sa])
            ivas.append(niva)
            ivbs.append(ivb)
        return ivas, ivbs
        
def getStripPsFdr(iva, ivb, model, N, win=5, triangular_flag=True):
    """
    for the interval a and b, searching its nearby windows to estimate FDR and p-values.  
    return ra, rb, rab, es, es_ra,es_rb, fdr, pop, nbp
    """
    ra, rb, rab = getPETsforRegions(iva, ivb, model, triangular_flag)
    ivas, ivbs = getNearbyPairRegionsForStrips(iva, ivb, win=win)
    #nras is a list for storing points ids for permutated regions
    nras, nrbs = [], []
    for na in ivas:
        nraSource = getCounts(na, model[0])
        if triangular_flag:
            nraTarget = getCounts(na, model[1])
            nra = nraSource.union(nraTarget)
        else:
            nra = nraSource
        nras.append(nra)
    for nb in ivbs:
        nrbTarget = getCounts(nb, model[1])
        if triangular_flag:
            nrbSource = getCounts(nb, model[0])
            nrb = nrbSource.union(nrbTarget)
        else:
            nrb = nrbTarget
        nrbs.append(nrb)
    #caculating the permutated background
    rabs, nbps = [], []
    for nra in nras:
        nralen = float(len(nra))
        for nrb in nrbs:
            nrblen = len(nrb)
            nrab = float(len(nra.intersection(nrb)))
            if nrab > 0:
                #collect the value for poisson test
                rabs.append(nrab)
                #collect the possibility for following binomial test
                den = nrab / (nralen * nrblen)
                nbps.append(den)
            else:
                nbps.append(0.0)
                rabs.append(0.0)
    if len(rabs) == 0:
        return ra, rb, rab, np.inf, rab / float(ra), rab / float(
            rb), 0.0, 0.0, 1e-300, 1e-300,
    rabs = np.array(rabs)
    #local fdr
    fdr = len(rabs[rabs > rab]) / float(len(rabs))
    mrabs = float(np.mean(rabs))
    #enrichment score
    if mrabs > 0:
        es = rab / np.mean(rabs[rabs > 0])
    else:
        es = np.inf
    #simple possion test
    lam = mrabs
    pop = max([1e-300, poisson.sf(rab - 1.0, lam)])
    #simple binomial test
    bp = np.mean(nbps) * ra * rb / N
    nbp = max([1e-300, binom.sf(rab - 1.0, N - rab, bp)])
    return ra, rb, rab, es, rab / float(ra), rab / float(rb), fdr, pop, nbp

def getStripSig(petfile, petclass, records):
    """
    """
    print "Starting estimate significance for %s candidate interactions in %s" % (len(records), petclass)
    model, N, intra_flag = getGenomeCoverage(petfile, petclass)
    print "Genomic coverage model built from %s" % petclass
    if N == 0:
        print "No PETs parsed as requiring distance cutoff >%s from %s" % (
            discut, petclass)
        return None
    ds = {}
    i = 0
    petclass_key = petclass.split("_")
    para = getPetFilePara(petfile)
    #determine whether contact matrix is a triangular matrix
    triangular_flag = para['symmetrical'] and isIntraPetClass(petclass_key, True)
    for r in records:
        chroms = [r[0], r[3]]
        key = "%s_%s" % (petclass, i)
        iva = [max(0, r[1]), r[2]]
        ivb = [max(0, r[4]), r[5]]
        ra, rb, rab = getPETsforRegions(iva, ivb, model, triangular_flag)
        if i % 100 == 0:
            cFlush("%s interaction p-values estimated for %s" % (i, petclass))
        ra, rb, rab, es, es_ra, es_rb, fdr, pop, nbp = getStripPsFdr(
            iva, ivb, model, N, 5, triangular_flag)
        #this part should be furthur modified, as for most ideable data, there are no noise, so the es should be inf, however, not possible
        ds[key] = {
            "ra": ra,
            "rb": rb,
            "rab": rab,
            "ES": es,
            "ES_ra": es_ra,
            "ES_rb": es_rb,
            "FDR": fdr,
            "poisson_p-value": pop,
            "binomial_p-value": nbp,
            "iva": "%s:%s-%s" % (chroms[0], iva[0], iva[1]),
            "ivb": "%s:%s-%s" % (chroms[1], ivb[0], ivb[1])
        }
        i += 1
    print
    #memory usage
    del model
    gc.collect()
    if len(ds.keys()) == 0:
        return None
    ds = pd.DataFrame(ds).T
    return ds
