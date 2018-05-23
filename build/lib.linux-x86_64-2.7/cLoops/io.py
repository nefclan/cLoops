#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
cLoops.io.py 
2018-05-22: re-design the data structure, all based on HDF5 now. 
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, random, gzip, glob

#3rd
import numpy as np
import h5py

#cLoops
from cLoops.utils import callSys
from cLoops.utils import cFlush


class PET(object):
    #cA is the center of left read
    __slots__ = [
        "chromA", "chromB", "startA", "startB", "endA", "endB", "strandA",
        "strandB", "cA", "cB", "distance"
    ]

    def __init__(self, d, symmetrical=True, use_strand=(False, False)):
        """
        d is line = line.split( "\n" )[ 0 ].split( "\t" ) from BEDPE file 
        """
        self.chromA = d[0]
        self.startA = int(d[1])
        self.endA = int(d[2])
        self.strandA = d[8]
        self.chromB = d[3]
        self.startB = int(d[4])
        self.endB = int(d[5])
        self.strandB = d[9]
        self.cA = (self.startA + self.endA) / 2
        self.cB = (self.startB + self.endB) / 2
        if self.chromA == self.chromB:
            #keep distance information
            self.distance = abs(self.cA - self.cB)
            #In symmetrical data, use_strand_num should be 0 or 2
            use_strand_num = sum(use_strand)
            if symmetrical:
                if use_strand_num == 0:
                    #make sure left end is alwasy smaller than right end
                    if self.startA + self.endA > self.startB + self.endB:
                        self.exchange()
                elif use_strand_num == 2:
                    #inter-strand interaction only have one pet class
                    #pet class: chrA_chrA_+_- or chrA_chrA_0_1
                    #'-' > '+' but 0 < 1
                    if (self.strandA > self.strandB):
                        self.exchange()
                    #intra-strand interaction
                    #left end is alwasy small than right end
                    elif (self.strandA == self.strandB and
                          self.startA + self.endA > self.startB + self.endB):
                        self.exchange()
        else:
            self.distance = None
            #smaller chrom number is always in left end if symmetrical
            if symmetrical and not chrOrdered(self.chromA, self.chromB):
                self.exchange()

    def exchange(self):
        self.startA, self.startB = self.startB, self.startA
        self.endA, self.endB = self.endB, self.endA
        self.strandA, self.strandB = self.strandB, self.strandA
        self.cA, self.cB = self.cB, self.cA


def chrOrdered(chrA, chrB):
    #determine which chromosome should be listed ahead
    baseNameA = chrA.lstrip('chr')
    baseNameB = chrB.lstrip('chr')
    flagA = baseNameA.isdigit()
    flagB = baseNameB.isdigit()
    if flagB:
        if not flagA or int(baseNameA) > int(baseNameB):
            return False
        else:
            return True
    elif not flagA:
        if baseA > baseB:
            return False
        else:
            return True
    else:
        return True


def isIntraPetClass(petclass_key, strict=True):
    #IntraPet: two ends of pets locate in same chromosome or
    #same strand if param strict set True
    if petclass_key[0] == petclass_key[1]:
        if not strict:
            return True
        n_strand = len(petclass_key) - 2
        if n_strand == 0:
            return True
        if n_strand == 2 and petclass_key[2] == petclass_key[3]:
            return True
    return False


def parseRawBedpe(fs,
                  petfile,
                  cs,
                  cut,
                  strict_intra,
                  logger,
                  symmetrical=True,
                  cis_only=True,
                  use_strand=False,
                  rmRep=False,
                  cal_ds=False):
    """
    Get the cis-PETs, organized by chromosomes. Input could be mixed PETs in bedpe.gz. Also change read id to numbers.
    @param fs: bedpe files of replicates, could be .bedpe or .bedpe.gz
    @param fout: intermediate HDF5 file name
    @param cs: chroms that wanted, list like ["chr1","chr2"]
    @param symmetrical: whether two ends of reads exchangeable
    @param cis_only: only use reads with two ends locate in same chromosome
    @param use_strand: Bool or list define whether strand infomation should
            be taken into account. If symmetrical is True, values in list
            must be same.
    @param rmRep: wehther remove redundant reads
    """
    if isinstance(use_strand, int):
        use_strand = (use_strand, use_strand)
    use_strand_num = sum(use_strand)
    if symmetrical and sum(use_strand) == 1:
        logger.error(
            "Error: symmetrical data could not only use one strand information"
        )
        return [], []
    #writing chunksize
    chunksize = 100000
    #pet data in classes
    pet_collection = {}
    #cis files
    cfs = []
    #distance between PETs mapped to different strands
    ds = []
    #create HDF5 file
    h5fh = h5py.File(petfile, 'w')
    h5pet_group = h5fh.create_group('pets')
    parameter_group = h5fh.create_group('parameters')
    parameter_group.create_dataset('symmetrical', data=symmetrical)
    parameter_group.create_dataset('cis_only', data=cis_only)
    parameter_group.create_dataset('use_strand', data=use_strand)
    i, j, = 0, 0
    for f in fs:
        r = "Parsing PETs from %s, requiring initial distance cutoff > %s" % (
            f, cut)
        logger.info(r)
        if f.endswith(".gz"):
            of = gzip.open(f, "rb")
        else:
            of = open(f)
        for line in of:
            i += 1
            if i % 100000 == 0:
                cFlush("%sk PETs processed from %s" % (i, f))
            line = line.rstrip("\n").split("\t")
            if "*" in line and "-1" in line:
                continue
            if len(line) < 6:
                continue
            try:
                pet = PET(line, symmetrical, use_strand)
            except:
                continue
            if cis_only:
                #cis reads
                if pet.chromA != pet.chromB:
                    continue
                #get PETs in required chroms
                if len(cs) > 0 and (not (pet.chromA in cs)):
                    continue
            else:
                #get PETs in required chroms
                if len(cs) > 0 and (not (pet.chromA in cs)
                                    and not (pet.chromB in cs)):
                    continue
            petclass_key = [pet.chromA, pet.chromB]
            if use_strand[0]:
                petclass_key.append(pet.strandA)
            if use_strand[1]:
                petclass_key.append(pet.strandB)
            #filtering too close cis PETs
            if cut > 0 and isIntraPetClass(
                    petclass_key, strict_intra) and pet.distance < cut:
                continue
            petclass = '_'.join(petclass_key)
            if petclass not in pet_collection:
                #c:counts, p:pets, r:uniq pet position set
                pet_collection[petclass] = {'c': 1, 'p': [], 'r': set()}
            #remove redundant reads in this petclass
            if rmRep and (pet.cA, pet.cB) in pet_collection[petclass]["r"]:
                continue
            pet_collection[petclass]['p'].append(
                [pet_collection[petclass]['c'], pet.cA, pet.cB])
            pet_collection[petclass]['c'] += 1
            if rmRep:
                pet_collection[petclass]['r'].add((pet.cA, pet.cB))
            j += 1
            if len(pet_collection[petclass]['p']) == chunksize:
                chunk = np.array(pet_collection[petclass]['p'])
                if petclass not in cfs:
                    cfs.append(petclass)
                    h5pet_group.create_dataset(
                        petclass,
                        shape=chunk.shape,
                        maxshape=(None, chunk.shape[1]),
                        chunks=chunk.shape,
                        dtype=chunk.dtype)
                    h5pet_group[petclass][:] = chunk
                else:
                    h5pet_group[petclass].resize(
                        pet_collection[petclass]['c'], axis=0)
                    h5pet_group[petclass][
                        pet_collection[petclass]['c'] - chunksize:] = chunk
                pet_collection[petclass]['p'] = []
            if cal_ds and pet.chromA == pet.chromB and pet.strandA != pet.strandB:
                ds.append(pet.distance)
        of.close()
    for petclass in pet_collection:
        if len(pet_collection[petclass]['p']) > 0:
            chunk = np.array(pet_collection[petclass]['p'])
            if petclass not in cfs:
                cfs.append(petclass)
                h5pet_group.create_dataset(petclass, data=chunk)
            else:
                h5pet_group[petclass].resize(
                    pet_collection[petclass]['c'], axis=0)
                h5pet_group[petclass][
                    pet_collection[petclass]['c'] - len(chunk):] = chunk
    h5fh.close()
    del pet_collection
    if cis_only:
        r = "Totaly %s PETs from %s, in which %s cis PETs" % (i, ",".join(fs),
                                                              j)
    else:
        r = "Totaly %s PETs from %s, in which %s useful PETs" % (i,
                                                                 ",".join(fs),
                                                                 j)
    logger.info(r)
    return cfs, ds


def getPetFilePara(fh):
    close_flag = False
    if not isinstance(fh, h5py._hl.files.File):
        fh = h5py.File(fh, 'r')
        close_flag = True
    p_group = fh['parameters']
    parameters = {}
    for key in p_group.keys():
        parameters[key] = p_group[key].value
    if close_flag:
        fh.close()
    return parameters


def getPetFileStrandSymbol(fh):
    close_flag = False
    if not isinstance(fh, h5py._hl.files.File):
        fh = h5py.File(fh, 'r')
        close_flag = True
    p_group = fh['pets']
    strand_symbol = set()
    for petclass in p_group.keys():
        strand_keys = petclass.split('_')[2:]
        for strand in strand_keys:
            strand_symbol.add(strand.encode('ascii', 'ignore'))
    if close_flag:
        fh.close()
    return strand_symbol


def readPets(petfile, petclass, cut=0, strict_intra=True):
    """
    read data from .pet file 
    """
    key = tuple(petclass.split('_'))
    intra_flag = isIntraPetClass(key, strict_intra)
    h5fh = h5py.File(petfile, 'r')
    if petclass not in h5fh['pets']:
        return key, np.ndarray(shape=(0, 3)), intra_flag
    mat = h5fh['pets'][petclass][:]
    if intra_flag and cut > 0:
        mat = mat[abs(mat[:, 1] - mat[:, 2]) >= cut, :]
    h5fh.close()
    return key, mat, intra_flag


def parseIv(iv):
    iv = [
        iv.split(":")[0],
        int(iv.split(":")[1].split("-")[0]),
        int(iv.split(":")[1].split("-")[1])
    ]
    return iv


def loops2washU(fin, fout, logger, significant=1):
    """
    Convert interaction level loop file to washU long range interactions. 
    Track format according to http://wiki.wubrowse.org/Long-range
    @param fin: interactions in loop file
    @param fout: washU  long-range interaction text file prefix
    @param significant: if set 1, only convert significant loops.
    """
    logger.info("Converting %s to washU long range interaction track." % fin)
    ss = {'0': 'minus', '1': 'plus', '-': 'minus', '+': 'plus'}
    filehandles = {}
    for i, line in enumerate(open(fin)):
        if i == 0:
            continue
        line = line.rstrip("\n").split("\t")
        #only using significant results
        if significant and float(line[-1]) < 1:
            continue

        key = line[0].split('_')
        del key[-1]
        strand_key = map(lambda x: ss[x], key[2:])
        if tuple(strand_key) not in filehandles:
            suffix_str = '_'.join([''] + strand_key)
            filehandles[tuple(strand_key)] = open(
                fout + suffix_str + '_loops_washU.txt', 'w')
        f = filehandles[tuple(strand_key)]

        #iva,ivb,ES
        nline = [line[6], line[7], line[1]]
        if nline[2] == 'inf':
            #washU wont accept inf
            nline[2] = '99999'
        iva = parseIv(nline[0])
        ivb = parseIv(nline[1])
        if iva[0] == ivb[0]:
            if iva[1] + iva[2] > ivb[1] + ivb[2]:
                nline[2] = '-' + nline[2]
        elif not chrOrdered(iva[0], ivb[0]):
            #washU won't show trans-chromosome loops, just in case
            nline[2] = '-' + nline[2]
        f.write("\t".join(map(str, nline)) + "\n")
    for f in filehandles.values():
        f.close()
    logger.info(
        "Converting %s to washU long range interaction track finished." % fin)


def loops2juice(fin, fout, logger, significant=1):
    """
    Convert interaction level loop file to Juicebox 2D annotation features. 
    The txt file format according to https://github.com/theaidenlab/juicebox/wiki/Loading-Annotations-(Annotations-menu)
    @param fin: interactions in loop file
    @param fout: juicebox long-range interaction text file prefix
    @param significant: if set 1, only convert significant loops.
    all p-values are -log10(p) transformed to escape all shown as 0 in juicebox.
    """
    logger.info("Converting %s to Juicebox 2D annotation feature." % fin)
    headerline = [
        "chromosome1", "x1", "x2", "chromosome2", "y1", "y2", "color",
        "observed", "loopId", "FDR", "EnrichmentScore", "distance",
        "-log10(binomial_p-value)", "-log10(poisson_p-value)",
        "-log10(hypergeometric_p-value)"
    ]
    ss = {'0': 'minus', '1': 'plus', '-': 'minus', '+': 'plus'}
    filehandles = {}
    for i, line in enumerate(open(fin)):
        if i == 0:
            continue
        line = line.rstrip("\n").split("\t")
        #only using significant results
        if significant and float(line[-1]) < 1:
            continue

        key = line[0].split('_')
        del key[-1]
        strand_key = map(lambda x: ss[x], key[2:])
        if tuple(strand_key) not in filehandles:
            suffix_str = '_'.join([''] + strand_key)
            filehandles[tuple(strand_key)] = open(
                fout + suffix_str + '_loops_juicebox.txt', 'w')
            f = filehandles[tuple(strand_key)]
            f.write("\t".join(headerline) + "\n")
        else:
            f = filehandles[tuple(strand_key)]

        iva = parseIv(line[6])
        ivb = parseIv(line[7])
        color = '"0,255,255"'
        #For assymetrical data loops with different direction are marked by color
        if iva[0] == ivb[0]:
            if iva[1] + iva[2] > ivb[1] + ivb[2]:
                color = '"0,0,255"'
        elif not chrOrdered(iva[0], ivb[0]):
            color = '"0,0,255"'
        try:
            nline = [
                iva[0], iva[1], iva[2], ivb[0], ivb[1], ivb[2], color,
                line[10], line[0], line[2], line[1], line[4],
                -np.log10(float(line[3])), -np.log10(float(line[8])),
                -np.log10(float(line[5]))
            ]
        except:
            continue
        f.write("\t".join(map(str, nline)) + "\n")
    for f in filehandles.values():
        f.close()
    logger.info(
        "Converting %s to Juicebox 2D annotation feature finished." % fin)


def stripes2juice(fin, fout, logger, significant=1):
    """
    Convert interaction level stripe file to Juicebox 2D annotation features. 
    The txt file format according to https://github.com/theaidenlab/juicebox/wiki/Loading-Annotations-(Annotations-menu)
    @param fin: interactions in loop file
    @param fout: juicebox long-range interaction text file prefix
    @param significant: if set 1, only convert significant loops.
    all p-values are -log10(p) transformed to escape all shown as 0 in juicebox.
    """
    logger.info("Converting %s to Juicebox 2D annotation feature." % fin)
    headerline = [
        "chromosome1",
        "x1",
        "x2",
        "chromosome2",
        "y1",
        "y2",
        "color",
        "observed",
        "stripeId",
        "FDR",
        "EnrichmentScore",
        "EnrichmentScore_X",
        "EnrichmentScore_Y",
        "-log10(binomal_p-value)",
        "-log10(poisson_p-value)",
    ]
    ss = {'0': 'minus', '1': 'plus', '-': 'minus', '+': 'plus'}
    filehandles = {}
    for i, line in enumerate(open(fin)):
        if i == 0:
            continue
        line = line.rstrip("\n").split("\t")
        if significant > 0 and float(line[12]) < 1:
            continue

        key = line[0].split('_')
        del key[-1]
        strand_key = map(lambda x: ss[x], key[2:])
        if tuple(strand_key) not in filehandles:
            suffix_str = '_'.join([''] + strand_key)
            filehandles[tuple(strand_key)] = open(
                fout + suffix_str + '_juicebox.txt', 'w')
            f = filehandles[tuple(strand_key)]
            f.write("\t".join(headerline) + "\n")
        else:
            f = filehandles[tuple(strand_key)]

        iva = parseIv(line[6])
        ivb = parseIv(line[7])
        nline = [
            iva[0],
            iva[1],
            iva[2],
            ivb[0],
            ivb[1],
            ivb[2],
            '"0,255,255"',
            line[10],
            line[0],
            line[4],
            line[1],
            line[2],
            line[3],
            -np.log10(float(line[5])),
            -np.log10(float(line[8])),
        ]
        f.write("\t".join(map(str, nline)) + "\n")
    for f in filehandles.values():
        f.close()
    logger.info(
        "Converting %s to Juicebox 2D annotation feature finished." % fin)


def pet2washU(fin, fout, cut, ext, strict_intra=True):
    """
    Convert PETs to washU long range interactions. 
    Track format according to http://wiki.wubrowse.org/Long-range
    @param fin: .pet file
    @param fout: prefix of output files
    @param cut: threshold to filter out pets
    @param ext: pet anchor region extension
    """
    print "Converting %s to washU track." % fin
    tmp = str(random.random())
    while glob.glob(tmp + '*'):
        tmp = str(random.random())
    #read pets
    h5fh = h5py.File(fin, 'r')
    h5pet_group = h5fh['pets']
    filehandles = {}
    suffix = []
    ss = {'0': 'minus', '1': 'plus', '-': 'minus', '+': 'plus'}
    for petclass in h5pet_group.keys():
        print "converting %s" % petclass
        mat = h5pet_group[petclass][:]
        class_str = petclass.encode('ascii', 'ignore')
        #key: chrA_chrB[_-][_+]
        key = class_str.split('_')
        #cutoff filtering
        if isIntraPetClass(key, strict_intra):
            mat = mat[abs(mat[:, 1] - mat[:, 2]) >= cut, :]
        strand_key = map(lambda x: ss[x], key[2:])
        if tuple(strand_key) not in filehandles:
            suffix_str = '_'.join([''] + strand_key)
            suffix.append(suffix_str)
            filehandles[tuple(strand_key)] = open(tmp + suffix_str, 'w')
        f = filehandles[tuple(strand_key)]
        for t in mat:
            a = (key[0], max([0, t[1] - ext]), t[1] + 1 + ext)
            b = (key[1], max([0, t[2] - ext]), t[2] + 1 + ext)
            score = 1
            #whether direction be presented in score
            #set color in washU browser to show direction
            #positive score direction: from left to right
            #negative score direction: from right to left
            if key[0] == key[1]:
                if t[1] > t[2]:
                    score = -1
            elif not chrOrdered(key[0], key[1]):
                #washU won't show trans-chromosome reads, just in case
                score = -1
            linea = [
                a[0], a[1], a[2],
                "%s:%s-%s,%s" % (b[0], b[1], b[2], score), t[0] * 2, "."
            ]
            lineb = [
                b[0], b[1], b[2],
                "%s:%s-%s,%s" % (a[0], a[1], a[2], score), t[0] * 2 + 1, "."
            ]
            f.write("\t".join(map(str, linea)) + "\n")
            f.write("\t".join(map(str, lineb)) + "\n")
    h5fh.close()
    for f in filehandles.values():
        f.close()
    commands = []
    for suffix_str in suffix:
        c1 = "bedtools sort -i %s > %s" % (tmp + suffix_str, fout + suffix_str)
        c2 = "rm %s" % tmp + suffix_str
        c3 = "bgzip %s" % fout + suffix_str
        c4 = "tabix -p bed %s.gz" % fout + suffix_str
        commands.extend([c1, c2, c3, c4])
    callSys(commands)
    print "Converting %s to washU random accessed track finished." % fout
    return 1


def pet2hic(fin, fout, cut, org, strict_intra=True):
    """
    Convert reads level bedpe to HIC.
    Track format according to https://github.com/theaidenlab/juicer/wiki/Pre#file-format
    @param fin: .pets file
    @param fout: prefix of output files
    """
    print "Converting %s to .hic file which could be loaded in juicebox" % fin
    tmp = str(random.random())
    while glob.glob(tmp + '*'):
        tmp = str(random.random())
    #read pets
    h5fh = h5py.File(fin, 'r')
    para = getPetFilePara(h5fh)
    h5pet_group = h5fh['pets']
    filehandles = {}
    suffix = []
    ss = {'0': 'minus', '1': 'plus', '-': 'minus', '+': 'plus'}
    for petclass in h5pet_group.keys():
        print "converting %s" % petclass
        mat = h5pet_group[petclass][:]
        class_str = petclass.encode('ascii', 'ignore')
        key = class_str.split('_')
        #cutoff filtering
        if isIntraPetClass(key, strict_intra):
            mat = mat[abs(mat[:, 1] - mat[:, 2]) >= cut, :]
        strand_key = map(lambda x: ss[x], key[2:])
        if para['symmetrical']:
            if tuple(strand_key) not in filehandles:
                suffix_str = '_'.join([''] + strand_key)
                suffix.append(suffix_str)
                filehandles[tuple(strand_key)] = open(tmp + suffix_str, 'w')
            f = filehandles[tuple(strand_key)]
            for t in mat:
                line = [0, key[0], t[1], 0, 1, key[1], t[2], 1]
                f.write("\t".join(map(str, line)) + "\n")
        else:
            if tuple(strand_key + ['forward']) not in filehandles:
                for direction in ['forward', 'backward']:
                    suffix_str = '_'.join([''] + strand_key + [direction])
                    suffix.append(suffix_str)
                    filehandles[tuple(strand_key + [direction])] = open(
                        tmp + suffix_str, 'w')
            f_forward = filehandles[tuple(strand_key + ['forward'])]
            f_backward = filehandles[tuple(strand_key + ['backward'])]
            if key[0] == key[1]:
                #intra-chromosome interaction
                for t in mat:
                    line = [0, key[0], t[1], 0, 1, key[1], t[2], 1]
                    if t[1] <= t[2]:
                        f_forward.write("\t".join(map(str, line)) + "\n")
                    else:
                        f_backward.write("\t".join(map(str, line)) + "\n")
            else:
                #inter-chromosome interaction
                if chrOrdered(key[0], key[1]):
                    f = f_forward
                else:
                    f = f_backward
                for t in mat:
                    line = [0, key[0], t[1], 0, 1, key[1], t[2], 1]
                    f.write("\t".join(map(str, line)) + "\n")
    h5fh.close()
    for f in filehandles.values():
        f.close()
    commands = []
    for suffix_str in suffix:
        if not para['cis_only']:
            c1 = "juicer_tools pre {fin} {fout} {org}".format(
                fin=tmp + suffix_str, fout=fout + suffix_str + '.hic', org=org)
        else:
            c1 = "juicer_tools pre -d {fin} {fout} {org}".format(
                fin=tmp + suffix_str, fout=fout + suffix_str + '.hic', org=org)
        c2 = "rm %s" % tmp + suffix_str
        #for debug only
        #c2 = ''
        commands.extend([c1, c2])
    callSys(commands)
    print "Converting %s to juicer's hic file finished." % fout
