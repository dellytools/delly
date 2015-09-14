#! /usr/bin/env python

from __future__ import print_function
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from varpkg.overlap import overlapValid, overlapMetrics
from varpkg.concordance import carrierConcordance
from varpkg.rd import rdAltRefRatio, validRdRatio
import vcf
import argparse
import collections
import banyan
import networkx
import operator


# Parse command line
parser = argparse.ArgumentParser(description='Inversion & proximal duplication classification.')
parser.add_argument('-v', '--vcf', metavar='inv.vcf', required=True, dest='invVCF', help='inversion or merged deletion/duplication vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-c', '--complexSV', metavar='complexSV.vcf', required=False, dest='complexVCF', help='complex SV vcf file (optional)')
parser.add_argument('-i', '--insOffset', metavar='170', required=False, dest='maxInsertionOffset', help='max. insertion offset (optional)')
parser.add_argument('-d', '--dupLength', metavar='150', required=False, dest='minDuplicationLength', help='min. duplication length (optional)')
parser.add_argument('-m', '--carrierConcordance', metavar='0.5', required=False, dest='minCarrierConcordance', help='min. carrier concordance (optional)')
parser.add_argument('-n', '--nestedOverlap', metavar='0.75', required=False, dest='minNestedOverlap', help='min. required nested overlap (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Parse command-line
maxInsertionOffset = 170
if args.maxInsertionOffset:
    maxInsertionOffset = int(args.maxInsertionOffset)
minDuplicationLength = 150
if args.minDuplicationLength:
    minDuplicationLength = int(args.minDuplicationLength)
minNestedOverlap = 0.75
if args.minNestedOverlap:
    minNestedOverlap = float(args.minNestedOverlap)
minCarrierConcordance = 0.5
if args.minCarrierConcordance:
    minCarrierConcordance = float(args.minCarrierConcordance)
maxSvSize = 50000

# Parse 3to3 inversions
sv = dict()
svDups = collections.defaultdict(list)
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ((record.INFO['CT'] == '3to3') or (record.INFO['SVTYPE'] == 'DUP')) and ((not args.siteFilter) or (len(record.FILTER) == 0)) and ((record.INFO['END']-record.POS+1) <= maxSvSize):
            precise = False
            if 'PRECISE' in record.INFO.keys():
                precise = record.INFO['PRECISE']
            nonRefHap = dict()
            support = 0
            rc = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    if call.gt_type != 0:
                        if precise:
                            if call['RV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                                support += call['RV']
                        else:
                            if call['DV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                                support += call['DV']
            if len(nonRefHap):
                if not sv.has_key(record.CHROM):
                    sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                if (record.POS, record.INFO['END']) not in sv[record.CHROM]:
                    sv[record.CHROM][(record.POS, record.INFO['END'])] = {'id': record.ID, 'sup': support, 'hap': nonRefHap, 'rc': rc}
                else:
                    svDups[(record.CHROM, record.POS, record.INFO['END'])].append({'id': record.ID, 'sup': support, 'hap': nonRefHap, 'rc': rc})


# Parse 5to5 inversions
invRegion = dict()
G = networkx.Graph()
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    for record in vcf_reader:
        if ((record.INFO['CT'] == '5to5') or (record.INFO['SVTYPE'] == 'DEL')) and ((not args.siteFilter) or (len(record.FILTER) == 0)) and (sv.has_key(record.CHROM)) and ((record.INFO['END']-record.POS+1) <= maxSvSize):
            precise = False
            if 'PRECISE' in record.INFO.keys():
                precise = record.INFO['PRECISE']
            nonRefHap = dict()
            support = 0
            rc = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    if call.gt_type != 0:
                        if precise:
                            if call['RV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                                support += call['RV']
                        else:
                            if call['DV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                                support += call['DV']
            if len(nonRefHap):
                # Collect overlapping calls
                s1 = record.POS
                e1 = record.INFO['END']
                invInfo = {'id': "NA", 'start': 0, 'end': 0, 'score': -1}
                for s2, e2 in sv[record.CHROM].overlap((s1, e1)):
                    for inv3to3 in svDups[(record.CHROM, s2, e2)] + [sv[record.CHROM][(s2, e2)]]:
                        (recO, nestedO, recUnion, bpOffset, oLen) = overlapMetrics((s1, e1), (s2, e2))
                        minBpOffset = min(abs(s2-s1), abs(e2-e1))
                        maxBpOffset = max(abs(s2-s1), abs(e2-e1))
                        cc = carrierConcordance(nonRefHap, inv3to3['hap'])
                        if (nestedO >= minNestedOverlap) and (minBpOffset < maxInsertionOffset) and (cc >= minCarrierConcordance):
                            rdRatio = rdAltRefRatio(((s1, e1), (s2, e2)), (nonRefHap, inv3to3['hap']), (rc, inv3to3['rc']))
                            valid, updSVType = validRdRatio(recO/nestedO, rdRatio)
                            if valid:
                                if (record.INFO['SVTYPE'] != 'INV') and (updSVType != "DUP"):
                                    continue
                                if (updSVType != 'DUP') or (maxBpOffset > minDuplicationLength):
                                    score = float(min(support, inv3to3['sup'])) * float(cc)
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
svtinv = "NA"
if args.invVCF:
    vcf_reader = vcf.Reader(open(args.invVCF), 'r', compressed=True) if args.invVCF.endswith('.gz') else vcf.Reader(open(args.invVCF), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if record.ID in set(reduce(operator.add, map(list, idPairs.keys()))):
            ci = max(abs(record.INFO['CIEND'][0]), abs(record.INFO['CIEND'][1]), abs(record.INFO['CIPOS'][0]), abs(record.INFO['CIPOS'][1]))
            svm = record.INFO['SVMETHOD']
            svtinv = record.INFO['SVTYPE']
            gtl = list()
            rc = dict()
            nonRefHap = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    if call.gt_type != 0:
                        if precise:
                            if call['RV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                        else:
                            if call['DV'] >= 2:
                                nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                gtl.append(call['GT'])
            selectedSVs[record.ID] = {'chr': record.CHROM, 'start': record.POS, 'end': record.INFO['END'], 'pe': record.INFO['PE'], 'mapq': record.INFO['MAPQ'], 'ci': ci, 'gt': gtl, 'rc': rc, 'nonRefHap': nonRefHap}
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
            (recO, nestedO, recUnion, bpOffset, oLen) = overlapMetrics((s1, e1), (s2, e2))
            rdRatio = rdAltRefRatio(((s1, e1), (s2, e2)), (selectedSVs[id1]['nonRefHap'], selectedSVs[id2]['nonRefHap']), (selectedSVs[id1]['rc'], selectedSVs[id2]['rc']))
            updSVType = validRdRatio(recO/nestedO, rdRatio)[1]

            # Separate the 2 SVs
            if updSVType == "DUP":
                if abs(s1-s2) > abs(e1-e2):
                    (dupStart, dupEnd, insStart, insEnd) = (min(s1, s2), max(s1, s2), min(e1, e2), max(e1, e2))
                else:
                    (insStart, insEnd, dupStart, dupEnd) = (min(s1, s2), max(s1, s2), min(e1, e2), max(e1, e2))
                svPos = dupStart
                svEnd = dupEnd
                svType = "PDUP"
                if svtinv == "INV":
                    svType = "INVDUP"
                sv2info = {'start': insStart, 'end': insEnd, 'type': "INS"}
            elif updSVType == "DEL":
                if abs(s1-s2) > abs(e1-e2):
                    (delStart, delEnd, invStart, invEnd) = (min(s1, s2), max(s1, s2), max(s1, s2), max(e1, e2))
                else:
                    (delStart, delEnd, invStart, invEnd) = (min(e1, e2), max(e1, e2), min(s1, s2), min(e1, e2))
                svPos = delStart
                svEnd = delEnd
                svType = "INVDEL"
                sv2info = {'start': invStart, 'end': invEnd, 'type': "INV"}
            elif updSVType == "INV":
                svPos = min(s1, s2)
                svEnd = max(e1, e2)
                svType = "INV"
                sv2info = None

            # output VCF
            info = "IMPRECISE;"
            info += "SVTYPE=" + updSVType + ";"
            info += "PE=" + str(selectedSVs[id1]['pe'] + selectedSVs[id2]['pe']) + ";"
            info += "MAPQ=" + str((selectedSVs[id1]['mapq'] + selectedSVs[id2]['mapq'])/2) + ";"
            info += "CIEND=" + str(-ci) + "," + str(ci) + ";CIPOS=" + str(-ci) + "," + str(ci) + ";SVMETHOD=" + svm + ";"
            info += "END=" + str(svEnd) + ";"
            info += "MC=" + id1 + "," + id2
            if sv2info is not None:
                info += ";SV2INFO=" + selectedSVs[id2]['chr'] + "," + str(sv2info['start']) + "," + str(sv2info['end']) + "," + sv2info['type']
            print(selectedSVs[id1]['chr'], svPos, updSVType + str(idCount).zfill(6), "N", "<" + svType + ">", ".", "PASS", info, "GT", sep="\t", file=f, end="")
            for (gt1, gt2) in zip(selectedSVs[id1]['gt'], selectedSVs[id2]['gt']):
                if (gt1 != None) and (gt1 == gt2):
                    print("\t" + gt1, end="", file=f)
                else:
                    print("\t./.", end="", file=f)
            print("", file=f)
            idCount += 1
