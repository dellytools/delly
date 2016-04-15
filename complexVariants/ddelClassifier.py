#! /usr/bin/env python

from __future__ import print_function
import sys, os
from varpkg.overlap import overlapValid, overlapMetrics
from varpkg.concordance import carrierConcordance
from varpkg.rd import altRefReadDepthRatio
import vcf
import argparse
import collections
import numpy
import banyan
import networkx
import operator
import re

# Parse command line
parser = argparse.ArgumentParser(description='Double-deletion classification.')
parser.add_argument('-v', '--vcf', metavar='inv.vcf', required=True, dest='invVCF', help='inversion/deletion vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-c', '--complexSV', metavar='complexSV.vcf', required=False, dest='complexVCF', help='complex SV vcf file (optional)')
parser.add_argument('-s', '--spacerLength', metavar='10000', required=False, dest='spacerLength', help='max. spacer length (optional)')
parser.add_argument('-d', '--deletionLength', metavar='100', required=False, dest='deletionLength', help='min. deletion/spacer length (optional)')
parser.add_argument('-m', '--carrierConcordance', metavar='0.5', required=False, dest='minCarrierConcordance', help='min. carrier concordance (optional)')
parser.add_argument('-p', '--reciprocalOverlap', metavar='0.5', required=False, dest='maxReciprocalOverlap', help='max. reciprocal overlap (optional)')
parser.add_argument('-r', '--readDepth', dest='readDepth', action='store_true', help='Filter sites according to read-depth')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Parse command-line
spacerLength = 10000
if args.spacerLength:
    spacerLength = int(args.spacerLength)
deletionLength = 100
if args.deletionLength:
    deletionLength = int(args.deletionLength)
minCarrierConcordance = 0.5
if args.minCarrierConcordance:
    minCarrierConcordance = float(args.minCarrierConcordance)
maxReciprocalOverlap = 0.5
if args.maxReciprocalOverlap:
    maxReciprocalOverlap = float(args.maxReciprocalOverlap)
maxInvSize = 50000

# Parse 3to3 inversions
sv = dict()
svDups = collections.defaultdict(list)
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ('CONTROL' in record.INFO.keys()) and (record.INFO['CONTROL'] != 0):
            continue
        if ((record.INFO['CT'] == '3to3') or (record.INFO['CT'] == '3to5')) and ((not args.siteFilter) or (len(record.FILTER) == 0)) and ((record.INFO['END']-record.POS+1) <= maxInvSize):
            hap = dict()
            nonRefHap = dict()
            rc = collections.defaultdict(int)
            peCount = 0
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    hap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                        peCount += call['DV']
            if len(nonRefHap):
                if not sv.has_key(record.CHROM):
                    sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                if (record.POS, record.INFO['END']) not in sv[record.CHROM]:
                    sv[record.CHROM][(record.POS, record.INFO['END'])] = {'id': record.ID, 'pe': peCount, 'nonRefHap': nonRefHap, 'rc': rc, 'hap': hap}
                else:
                    svDups[(record.CHROM, record.POS, record.INFO['END'])].append({'id': record.ID, 'pe': peCount, 'nonRefHap': nonRefHap, 'rc': rc, 'hap': hap})

# Parse valid control regions
control = dict()
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ('CONTROL' not in record.INFO.keys()) or (record.INFO['CONTROL'] == 0):
            continue
        svSize = record.INFO['END'] - record.POS
        if (record.CHROM in sv.keys()) and (len(sv[record.CHROM].overlap((record.POS + int(0.25 * svSize), record.INFO['END'] - int(0.25 * svSize))))):
            continue
        if record.INFO['CONTROL'] not in control.keys():
            control[record.INFO['CONTROL']] = collections.defaultdict(int)
        for call in record.samples:
            if call.called:
                control[record.INFO['CONTROL']][call.sample] += call['RC']

# Parse 5to5 inversions
invRegion = dict()
G = networkx.Graph()
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ('CONTROL' in record.INFO.keys()) and (record.INFO['CONTROL'] != 0):
            continue
        if ((record.INFO['CT'] == '5to5') or (record.INFO['CT'] == '3to5')) and ((not args.siteFilter) or (len(record.FILTER) == 0)) and (sv.has_key(record.CHROM)) and ((record.INFO['END']-record.POS+1) <= maxInvSize):
            nonRefHap = dict()
            peCount = 0
            rc = dict()
            hap = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    hap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                        peCount += call['DV']
            if len(nonRefHap):
                # Collect overlapping calls
                s1 = record.POS
                e1 = record.INFO['END']
                invInfo = {'id': "NA", 'start': 0, 'end': 0, 'score': -1}
                padding = 0
                if record.INFO['SVTYPE'] == "DEL":
                    padding = spacerLength
                for s2, e2 in sv[record.CHROM].overlap((max(s1 - padding, 0), e1 + padding)):
                     for inv3to3 in svDups[(record.CHROM, s2, e2)] + [sv[record.CHROM][(s2, e2)]]:
                        if inv3to3['id'] == record.ID:
                            continue
                        spacer = min(abs(e1-s2), abs(e2-s1))
                        (recO, nestedO, recUnion, bpOffset, oLen) = overlapMetrics((s1, e1), (s2, e2))
                        cc = carrierConcordance(nonRefHap, inv3to3['nonRefHap'])
                        delLength = min(e1-s1, e2-s2)
                        if record.INFO['SVTYPE'] == "INV":
                            delLength = min(max(s1,s2)-min(s1,s2), max(e1,e2)-min(e1,e2))
                        if (spacer <= spacerLength) and (cc >= minCarrierConcordance) and (min(spacer, delLength) >= deletionLength) and (recO <= maxReciprocalOverlap):
                            svControlID1 = re.sub(r"^[A-Z0]*","", record.ID)
                            rdRatio1 = 1.0
                            if int(svControlID1) in control.keys():
                                rdRatio1 = altRefReadDepthRatio(rc, control[int(svControlID1)], hap)
                            svControlID2 = re.sub(r"^[A-Z0]*","", inv3to3['id'])
                            rdRatio2 = 1.0
                            if int(svControlID2) in control.keys():
                                rdRatio2 = altRefReadDepthRatio(inv3to3['rc'], control[int(svControlID2)], inv3to3['hap'])
                            print(record.CHROM, record.POS, record.INFO['END'], record.ID, record.CHROM, s2, e2, inv3to3['id'], spacer, delLength, cc, rdRatio1, rdRatio2)
                            if (not args.readDepth) or ((rdRatio1<0.8) and (rdRatio2<0.8)):
                                score = float(min(peCount, inv3to3['pe'])) * float(cc)
                                if score > invInfo['score']:
                                    invInfo = {'id': inv3to3['id'], 'start': min(s1, s2), 'end': max(e1, e2), 'score': score}
                if invInfo['score'] >= 0:
                    if not invRegion.has_key(record.CHROM):
                        invRegion[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    G.add_node((record.ID, invInfo['id']))
                    G.node[(record.ID, invInfo['id'])]['Score'] = invInfo['score']
                    for invIStart, invIEnd in invRegion[record.CHROM].overlap((invInfo['start'], invInfo['end'])):
                        (id1, id2) = invRegion[record.CHROM][(invIStart, invIEnd)]
                        if overlapValid((invInfo['start'], invInfo['end']), (invIStart, invIEnd), 0.1, 10000):
                            G.add_edge((record.ID, invInfo['id']), (id1, id2))
                    invRegion[record.CHROM][(invInfo['start'], invInfo['end'])] = (record.ID, invInfo['id'])

# Pick best pair of inversions out of all overlapping calls
idPairs = dict()
for H in networkx.connected_component_subgraphs(G):
    bestScore = -1.0
    for n, d in H.nodes_iter(data=True):
        if d['Score'] > bestScore:
            bestScore = d['Score']
            bestSVs = n
    idPairs[bestSVs] = bestScore

# Extract selected calls
selectedSVs = dict()
svm = "NA"
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if record.ID in set(svID for pair in idPairs.keys() for svID in pair):
            ci = max(abs(record.INFO['CIEND'][0]), abs(record.INFO['CIEND'][1]), abs(record.INFO['CIPOS'][0]), abs(record.INFO['CIPOS'][1]))
            svm = record.INFO['SVMETHOD']
            gtl = list()
            for call in record.samples:
                gtl.append(call['GT'])
            selectedSVs[record.ID] = {'chr': record.CHROM, 'start': record.POS, 'end': record.INFO['END'], 'pe': record.INFO['PE'], 'mapq': record.INFO['MAPQ'], 'ci': ci, 'gt': gtl}
            vcf_writer.write_record(record)

# Create complex SV VCF file
if args.complexVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    vcf_reader.infos['SV2INFO'] = vcf.parser._Info('SV2INFO', 4, 'String', 'Second SV information of the form chr, start, end, type')
    vcf_reader.infos['MC'] = vcf.parser._Info('MC', '.', 'String', 'Merged calls.')
    vcf_writer = vcf.Writer(open(args.complexVCF, 'w'), vcf_reader, lineterminator='\n')
    vcf_writer.close()
    with open(args.complexVCF, 'a') as f:
        idCount = 1
        for (id1, id2) in idPairs.keys():
            ci = max(selectedSVs[id1]['ci'], selectedSVs[id2]['ci'])
            ((s1, e1), (s2, e2)) = ((selectedSVs[id1]['start'], selectedSVs[id1]['end']), (selectedSVs[id2]['start'], selectedSVs[id2]['end']))

            #print(selectedSVs[id1]['chr'], min(s1,s2), max(e1,e2), "DDel" + str(idCount), sep='\t')
            svType = "DDEL"
            updSVType="DEL"
            if id1[0:3] == "INV":
                sv2info = {'start': max(s1,s2), 'end': min(e1,e2), 'type': "INV"}
            else:
                sv2info = {'start': min(e1,e2), 'end': max(s1,s2), 'type': "INS"}
                #sv2info = None

            # output VCF
            info = "IMPRECISE;"
            info += "SVTYPE=" + updSVType + ";"
            info += "PE=" + str(selectedSVs[id1]['pe'] + selectedSVs[id2]['pe']) + ";"
            info += "MAPQ=" + str((selectedSVs[id1]['mapq'] + selectedSVs[id2]['mapq'])/2) + ";"
            info += "CIEND=" + str(-ci) + "," + str(ci) + ";CIPOS=" + str(-ci) + "," + str(ci) + ";SVMETHOD=" + svm + ";"
            info += "MC=" + id1 + "," + id2 + ";"
            if sv2info is not None:
                info += "SV2INFO=" + selectedSVs[id2]['chr'] + "," + str(sv2info['start']) + "," + str(sv2info['end']) + "," + sv2info['type'] + ";"

            # Del1
            svPos = min(s1,s2)
            if id1[0:3] == "INV":
                svEnd = max(s1,s2)
            else:
                svEnd = min(e1,e2)
            info1 = info + "END=" + str(svEnd)
            print(selectedSVs[id1]['chr'], svPos, updSVType + str(idCount).zfill(6), "N", "<" + svType + ">", ".", "PASS", info1, "GT", sep="\t", file=f, end="")
            for (gt1, gt2) in zip(selectedSVs[id1]['gt'], selectedSVs[id2]['gt']):
                if (gt1 != None) and (gt1 == gt2):
                    print("\t" + gt1, end="", file=f)
                else:
                    print("\t./.", end="", file=f)
            print("", file=f)
            idCount += 1

            # Del2
            svEnd = max(e1,e2)
            if id1[0:3] == "INV":
                svPos = min(e1,e2)
            else:
                svPos = max(s1,s2)
            info2 = info + "END=" + str(svEnd)
            print(selectedSVs[id1]['chr'], svPos, updSVType + str(idCount).zfill(6), "N", "<" + svType + ">", ".", "PASS", info2, "GT", sep="\t", file=f, end="")
            for (gt1, gt2) in zip(selectedSVs[id1]['gt'], selectedSVs[id2]['gt']):
                if (gt1 != None) and (gt1 == gt2):
                    print("\t" + gt1, end="", file=f)
                else:
                    print("\t./.", end="", file=f)
            print("", file=f)
            idCount += 1

