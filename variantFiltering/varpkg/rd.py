#! /usr/bin/env python

from __future__ import print_function
import numpy

def altRefReadDepthRatio(sv1RC, sv2RC, geno):
    rcSamples = set(sv1RC.keys()).intersection(set(sv2RC.keys()))
    hetRC = list()
    refRC = list()
    for sp in rcSamples:
        if (sv1RC[sp] > 0) and (sv2RC[sp] > 0):
            if sum(geno[sp]) == 0:
                refRC.append(float(sv1RC[sp])/float(sv2RC[sp]))
            elif sum(geno[sp]) == 1:
                hetRC.append(float(sv1RC[sp])/float(sv2RC[sp]))
    if (len(hetRC)) and (len(refRC)):
        return numpy.median(numpy.array(hetRC))/numpy.median(numpy.array(refRC))
    else:
        return None


def rdAltRefRatio(((start1, end1), (start2, end2)), (sv1hap, sv2hap), (sv1RC, sv2RC)):
    sv1samples = set(sv1hap.keys())
    sv2samples = set(sv2hap.keys())
    altRC = list()
    refRC = list()
    for sp in set(sv1RC.keys()).intersection(set(sv2RC.keys())):
        if (sv1RC[sp] > 0) and (sv2RC[sp] > 0):
            if sp in sv1samples.intersection(sv2samples):
                if (end1-start1) > (end2-start2):
                    altRC.append(float(sv1RC[sp])/float(sv2RC[sp]))
                else:
                    altRC.append(float(sv2RC[sp])/float(sv1RC[sp]))
            elif sp not in sv1samples.union(sv2samples):
                if (end1-start1) > (end2-start2):
                    refRC.append(float(sv1RC[sp])/float(sv2RC[sp]))
                else:
                    refRC.append(float(sv2RC[sp])/float(sv1RC[sp]))
    if (len(altRC)) and (len(refRC)):
        return numpy.median(numpy.array(altRC))/numpy.median(numpy.array(refRC))
    else:
        return None



def validRdRatio(ro, readDepthRatio, checkReadDepth):
    if not checkReadDepth:
        return (True, None)
    if readDepthRatio is not None:
        if ((ro < 0.9) and (readDepthRatio <= 0.925)) or ((ro < 0.95) and (readDepthRatio >= 1.025)) or ((ro >= 0.9) and (readDepthRatio < 1.025)) or (ro > 0.95):
            if (ro < 0.9) and (readDepthRatio <= 0.925):
                return(True, "DEL")
            elif (ro < 0.95) and (readDepthRatio >= 1.025):
                return(True, "DUP")
            else:
                return(True, "INV")
        else:
            return (False, None)
    else:
        return (False, None)
