#! /usr/bin/env python

from __future__ import print_function
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from varpkg.overlap import overlapValid
import vcf
import argparse
import numpy
import banyan
import networkx

# Parse command line
parser = argparse.ArgumentParser(description='Deletion/Duplication filter.')
parser.add_argument('-v', '--vcf', metavar='cnv.vcf', required=True, dest='cnvVCF', help='deletion/duplication vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-q', '--quality', type=float, default=0.6, metavar='0.6', required=False, dest='quality', help='required quality [0,1] (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Parse command-line
quality = args.quality

# Parse deletions/duplications
cnvRegion = dict()
G = networkx.Graph()
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    for record in vcf_reader:
        if (not args.siteFilter) or (len(record.FILTER) == 0):
            precise = False
            if 'PRECISE' in record.INFO.keys():
                precise = record.INFO['PRECISE']
            rc = list()
            hetRC = list()
            refRC = list()
            gqAlt = list()
            gqRef = list()
            ratioRef = [0]
            ratioAlt = []
            support = 0
            refpass = False
            altpass = False
            for call in record.samples:
                if call.called:
                    hap = [int(gVal) for gVal in call['GT'].split('/')]
                    if sum(hap) == 0:
                        if call['FT'] == "PASS":
                            refpass = True
                        gqRef.append(call['GQ'])
                        if precise:
                            ratioRef.append(float(call['RV'])/float(call['RR'] + call['RV']))
                        else:
                            ratioRef.append(float(call['DV'])/float(call['DR'] + call['DV']))
                        if call['RC'] > 0:
                            rc.append(call['RC'])
                            if call['RCL'] + call['RCR'] > 0:
                                refRC.append(float(call['RC'])/float(call['RCL'] + call['RCR']))
                    else:
                        if call['FT'] == "PASS":
                            altpass = True
                        gqAlt.append(call['GQ'])
                        if precise:
                            ratioAlt.append(float(call['RV'])/float(call['RR'] + call['RV']))
                            support += call['RV']
                        else:
                            ratioAlt.append(float(call['DV'])/float(call['DR'] + call['DV']))
                            support += call['DV']
                        if (sum(hap) == 1) and (call['RCL'] + call['RCR'] > 0):
                            hetRC.append(float(call['RC'])/float(call['RCL'] + call['RCR']))
            callRate = float(len(gqRef) + len(gqAlt)) / float(len(record.samples))
            if (len(gqAlt)) and (len(gqRef)) and (refpass) and (altpass) and (callRate >= 0.75):
                svStart = record.POS
                svEnd = record.INFO['END']
                altgq = numpy.median(numpy.array(gqAlt))
                refgq = numpy.median(numpy.array(gqRef))
                altratio = numpy.median(numpy.array(ratioAlt))
                # Normalize values
                if altgq > 200:
                    altgq = float(1)
                else:
                    altgq = float(altgq) / float(200)
                if refgq > 200:
                    refgq = float(1)
                else:
                    refgq = float(refgq) / float(200)
                validRdRatio = False
                qIndex = (altratio + altgq + refgq) / 3
                if (record.INFO['SVTYPE'] == "DEL") or (record.INFO['SVTYPE'] == "DUP"):
                    if (record.INFO['END'] - record.POS < 500) and (precise):
                        if (record.INFO['SRQ']>0.9):
                            # No read-depth check
                            validRdRatio = True
                    elif (len(hetRC)) and (len(refRC)):
                        rcref = numpy.median(numpy.array(rc))
                        if rcref > 500:
                            rcref = float(1)
                        else:
                            rcref = float(rcref) / float(500)
                        qIndex = (rcref + altratio + altgq + refgq) / 4
                        rdRatio = numpy.median(numpy.array(hetRC))/numpy.median(numpy.array(refRC))
                        if ((record.INFO['SVTYPE'] == "DEL") and (rdRatio < 0.8)) or ((record.INFO['SVTYPE'] == "DUP") and (rdRatio >= 1.3) and (rdRatio <= 1.75)):
                            validRdRatio = True
                else:
                    validRdRatio = True

                # Check quality
                #print(record.CHROM, svStart, svEnd, record.ID, qIndex, numpy.percentile(ratioRef, 99), altgq, refgq, altratio, sep="\t")
                if (validRdRatio) and (qIndex > quality) and (numpy.percentile(ratioRef, 99) == 0):
                    if not cnvRegion.has_key(record.CHROM):
                        cnvRegion[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    G.add_node(record.ID)
                    G.node[record.ID]['Score'] = support
                    for cnvIStart, cnvIEnd in cnvRegion[record.CHROM].overlap((svStart, svEnd)):
                        otherID = cnvRegion[record.CHROM][(cnvIStart, cnvIEnd)]
                        if (record.INFO['SVTYPE'] == "INS") or (overlapValid((svStart, svEnd), (cnvIStart, cnvIEnd), 0.1, 10000)):
                            G.add_edge(record.ID, otherID)
                    cnvRegion[record.CHROM][(svStart - 15, svEnd + 15)] = record.ID  # padding for PRECISE insertion

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
