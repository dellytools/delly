#! /usr/bin/env python

from __future__ import print_function

# Returns sample carrier concordance
def carrierConcordance(sv1hap, sv2hap):
    sv1samples = set(sv1hap.keys())
    sv2samples = set(sv2hap.keys())
    intersectSamples = sv1samples.intersection(sv2samples)
    denominator = len(sv1samples.union(sv2samples))
    return float(len(intersectSamples))/float(denominator)



if __name__ == "__main__":
    import sys
    print("Enter carrier samples for SV1: s1, s2, s3, s4, ...")
    svA = dict.fromkeys([x.strip() for x in sys.stdin.readline().split(',')])
    print("Enter carrier samples for SV2: s1, s2, s3, s4, ...")
    svB = dict.fromkeys([x.strip() for x in sys.stdin.readline().split(',')])
    print(carrierConcordance(svA, svB))
