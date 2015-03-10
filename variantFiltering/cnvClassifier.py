#! /usr/bin/env python

from __future__ import print_function
from varpkg.overlap import overlapValid
import vcf
import argparse
import numpy
import banyan
import networkx
import collections
import re

def carrierConcordance(sv1hap, sv2hap):
    sv1samples = set(sv1hap.keys())
    sv2samples = set(sv2hap.keys())
    intersectSamples = sv1samples.intersection(sv2samples)
    denominator = len(sv1samples.union(sv2samples))
    return float(len(intersectSamples))/float(denominator)

def rdAltRefRatio(sv1RC, sv2RC, geno):
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


# Parse command line
parser = argparse.ArgumentParser(description='Deletion/Duplication filter.')
parser.add_argument('-v', '--vcf', metavar='cnv.vcf', required=True, dest='cnvVCF', help='deletion/duplication vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-r', '--readDepth', dest='readDepth', action='store_true', help='Filter sites according to read-depth')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Parse control deletions/duplications
sv = dict()
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ('CONTROL' not in record.INFO.keys()) or (record.INFO['CONTROL'] == 0):
            continue
        if record.INFO['CONTROL'] not in sv.keys():
            sv[record.INFO['CONTROL']] = collections.defaultdict(int)
        for call in record.samples:
            if call.called:
                sv[record.INFO['CONTROL']][call.sample] += call['RC']

# Parse deletions/duplications
cnvRegion = dict()
G = networkx.Graph()
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ('CONTROL' in record.INFO.keys()) and (record.INFO['CONTROL'] != 0):
            continue
        if (not args.siteFilter) or (len(record.FILTER) == 0):
            hap = dict()
            rc = collections.defaultdict(int)
            peCount = 0
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    hap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        peCount += call['DV']
            if len(hap):
                svStart = record.POS
                svEnd = record.INFO['END']
                svControlID = re.sub(r"^[A-Z0]*","", record.ID)
                rdRatio = rdAltRefRatio(rc, sv[int(svControlID)], hap)
                #print(record.CHROM, svStart, svEnd, record.ID, rdRatio, sep="\t")
                if rdRatio is not None:
                    if ((record.INFO['SVTYPE'] == "DEL") and (rdRatio < 0.8)) or ((record.INFO['SVTYPE'] == "DUP") and (rdRatio > 1.15)):
                        # Valid Call
                        if not cnvRegion.has_key(record.CHROM):
                            cnvRegion[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                        G.add_node(record.ID)
                        G.node[record.ID]['Score'] = peCount
                        for cnvIStart, cnvIEnd in cnvRegion[record.CHROM].overlap((svStart, svEnd)):
                            otherID = cnvRegion[record.CHROM][(cnvIStart, cnvIEnd)]
                            if overlapValid((svStart, svEnd), (cnvIStart, cnvIEnd), 0.1, 10000):
                                G.add_edge(record.ID, otherID)
                        cnvRegion[record.CHROM][(svStart, svEnd)] = record.ID

# Pick best deletion/duplication for all overlapping calls
selectedSVs = set()
for H in networkx.connected_component_subgraphs(G):
    bestScore = -1.0
    for n, d in H.nodes_iter(data=True):
        if d['Score'] > bestScore:
            bestScore = d['Score']
            bestSV = n
    if bestSV is not None:
        selectedSVs.add(bestSV)

# Extract selected calls
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if record.ID in selectedSVs:
            vcf_writer.write_record(record)
