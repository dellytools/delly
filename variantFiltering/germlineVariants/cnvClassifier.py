#! /usr/bin/env python

from __future__ import print_function
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from varpkg.overlap import overlapValid
from varpkg.concordance import carrierConcordance
from varpkg.rd import altRefReadDepthRatio
import vcf
import argparse
import numpy
import banyan
import networkx
import collections
import re

# Parse command line
parser = argparse.ArgumentParser(description='Deletion/Duplication filter.')
parser.add_argument('-v', '--vcf', metavar='cnv.vcf', required=True, dest='cnvVCF', help='deletion/duplication vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Parse deletions/duplications
cnvRegion = dict()
G = networkx.Graph()
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    for record in vcf_reader:
        if (not args.siteFilter) or (len(record.FILTER) == 0):
            hetRC = list()
            refRC = list()
            peCount = 0
            for call in record.samples:
                if call.called:
                    if (call['RC'] > 0) and (call['RCL'] + call['RCR'] > 0):
                        hap = [int(gVal) for gVal in call['GT'].split('/')]
                        if sum(hap) == 0:
                            refRC.append(float(call['RC'])/float(call['RCL'] + call['RCR']))
                        elif sum(hap) == 1:
                            hetRC.append(float(call['RC'])/float(call['RCL'] + call['RCR']))
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        peCount += call['DV']
            if (len(hetRC)) and (len(refRC)):
                svStart = record.POS
                svEnd = record.INFO['END']
                rdRatio = numpy.median(numpy.array(hetRC))/numpy.median(numpy.array(refRC))
                #print(record.CHROM, svStart, svEnd, record.ID, rdRatio, sep="\t")
                if ((record.INFO['SVTYPE'] == "DEL") and (rdRatio < 0.8)) or ((record.INFO['SVTYPE'] == "DUP") and (rdRatio >= 1.3) and (rdRatio <= 1.75)):
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
