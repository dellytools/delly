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
parser = argparse.ArgumentParser(description='Proximal duplication classification.')
parser.add_argument('-v', '--delVCF', metavar='del.vcf', required=True, dest='delVCF', help='deletion vcf file (required)')
parser.add_argument('-w', '--dupVCF', metavar='dup.vcf', required=True, dest='dupVCF', help='duplication vcf file (required)')
parser.add_argument('-o', '--outDEL', metavar='out.del.vcf', required=True, dest='outDEL', help='output del vcf file (required)')
parser.add_argument('-p', '--outDUP', metavar='out.dup.vcf', required=True, dest='outDUP', help='output dup vcf file (required)')
parser.add_argument('-c', '--complexSV', metavar='complexSV.vcf', required=False, dest='complexVCF', help='complex SV vcf file (optional)')
parser.add_argument('-i', '--insOffset', metavar='170', required=False, dest='maxInsertionOffset', help='max. insertion offset (optional)')
parser.add_argument('-d', '--dupLength', metavar='150', required=False, dest='minDuplicationLength', help='min. duplication length (optional)')
parser.add_argument('-m', '--carrierConcordance', metavar='0.5', required=False, dest='minCarrierConcordance', help='min. carrier concordance (optional)')
parser.add_argument('-n', '--nestedOverlap', metavar='0.75', required=False, dest='minNestedOverlap', help='min. required nested overlap (optional)')
parser.add_argument('-r', '--readDepth', dest='readDepth', action='store_true', help='Filter sites according to read-depth')
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


# Get sample intersection
sampleSet = set()
if (args.delVCF) and (args.dupVCF):
    vcf_reader1 = vcf.Reader(open(args.delVCF), 'r', compressed=True) if args.delVCF.endswith('.gz') else vcf.Reader(open(args.delVCF), 'r', compressed=False)
    vcf_reader2 = vcf.Reader(open(args.dupVCF), 'r', compressed=True) if args.dupVCF.endswith('.gz') else vcf.Reader(open(args.dupVCF), 'r', compressed=False)
    sampleSet = set(vcf_reader1.samples).intersection(vcf_reader2.samples)

# Parse duplications
sv = dict()
svDups = collections.defaultdict(list)
if args.dupVCF:
    vcf_reader = vcf.Reader(open(args.dupVCF), 'r', compressed=True) if args.dupVCF.endswith('.gz') else vcf.Reader(open(args.dupVCF), 'r', compressed=False)
    for record in vcf_reader:
        if (record.INFO['SVTYPE'] == 'DUP') and ((not args.siteFilter) or (len(record.FILTER) == 0)) and ((record.INFO['END']-record.POS+1) <= maxSvSize):
            nonRefHap = dict()
            peCount = 0
            rc = dict()
            for call in record.samples:
                if (call.called) and (call.sample in sampleSet):
                    rc[call.sample] = call['RC']
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                        peCount += call['DV']
            if len(nonRefHap):
                if not sv.has_key(record.CHROM):
                    sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                if (record.POS, record.INFO['END']) not in sv[record.CHROM]:
                    sv[record.CHROM][(record.POS, record.INFO['END'])] = {'id': record.ID, 'pe': peCount, 'hap': nonRefHap, 'rc': rc}
                else:
                    svDups[(record.CHROM, record.POS, record.INFO['END'])].append({'id': record.ID, 'pe': peCount, 'hap': nonRefHap, 'rc': rc})

# Parse deletions VCF
dupRegion = dict()
G = networkx.Graph()
if args.delVCF:
    vcf_reader = vcf.Reader(open(args.delVCF), 'r', compressed=True) if args.delVCF.endswith('.gz') else vcf.Reader(open(args.delVCF), 'r', compressed=False)
    for record in vcf_reader:
        if (record.INFO['SVTYPE'] == 'DEL') and ((not args.siteFilter) or (len(record.FILTER) == 0)) and (sv.has_key(record.CHROM)) and ((record.INFO['END']-record.POS+1) <= maxSvSize):
            nonRefHap = dict()
            peCount = 0
            rc = dict()
            for call in record.samples:
                if (call.called) and (call.sample in sampleSet):
                    rc[call.sample] = call['RC']
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                        peCount += call['DV']
            if len(nonRefHap):
                # Collect overlapping calls
                s1 = record.POS
                e1 = record.INFO['END']
                dupInfo = {'id': "NA", 'start': 0, 'end': 0, 'score': -1}
                for s2, e2 in sv[record.CHROM].overlap((s1, e1)):
                    for dup5to3 in svDups[(record.CHROM, s2, e2)] + [sv[record.CHROM][(s2, e2)]]:
                        (recO, nestedO, recUnion, bpOffset, oLen) = overlapMetrics((s1, e1), (s2, e2))
                        minBpOffset = min(abs(s2-s1), abs(e2-e1))
                        maxBpOffset = max(abs(s2-s1), abs(e2-e1))
                        cc = carrierConcordance(nonRefHap, dup5to3['hap'])
                        if (nestedO >= minNestedOverlap) and (minBpOffset < maxInsertionOffset) and (maxBpOffset > minDuplicationLength) and (cc >= minCarrierConcordance):
                            rdRatio = rdAltRefRatio(((s1, e1), (s2, e2)), (nonRefHap, dup5to3['hap']), (rc, dup5to3['rc']))
                            if validRdRatio(recO/nestedO, rdRatio, args.readDepth)[0]:
                                score = float(min(peCount, dup5to3['pe'])) * float(cc)
                                if score > dupInfo['score']:
                                    dupInfo = {'id': dup5to3['id'], 'start': min(s1, s2), 'end': max(e1, e2), 'score': score}
                if dupInfo['score'] >= 0:
                    if not dupRegion.has_key(record.CHROM):
                        dupRegion[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    G.add_node((record.ID, dupInfo['id']))
                    G.node[(record.ID, dupInfo['id'])]['Score'] = dupInfo['score']
                    for dupIStart, dupIEnd in dupRegion[record.CHROM].overlap((dupInfo['start'], dupInfo['end'])):
                        (id1, id2) = dupRegion[record.CHROM][(dupIStart, dupIEnd)]
                        if overlapValid((dupInfo['start'], dupInfo['end']), (dupIStart, dupIEnd), 0.1, 10000):
                            G.add_edge((record.ID, dupInfo['id']), (id1, id2))
                    dupRegion[record.CHROM][(dupInfo['start'], dupInfo['end'])] = (record.ID, dupInfo['id'])

# Pick best pair
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
if args.delVCF:
    vcf_reader = vcf.Reader(open(args.delVCF), 'r', compressed=True) if args.delVCF.endswith('.gz') else vcf.Reader(open(args.delVCF), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outDEL, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if record.ID in set(reduce(operator.add, map(list, idPairs.keys()))):
            ci = max(abs(record.INFO['CIEND'][0]), abs(record.INFO['CIEND'][1]), abs(record.INFO['CIPOS'][0]), abs(record.INFO['CIPOS'][1]))
            svm = record.INFO['SVMETHOD']
            gtl = list()
            rc = dict()
            nonRefHap = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                gtl.append(call['GT'])
            selectedSVs[record.ID] = {'chr': record.CHROM, 'start': record.POS, 'end': record.INFO['END'], 'pe': record.INFO['PE'], 'mapq': record.INFO['MAPQ'], 'ci': ci, 'gt': gtl, 'rc': rc, 'nonRefHap': nonRefHap}
            vcf_writer.write_record(record)
if args.dupVCF:
    vcf_reader = vcf.Reader(open(args.dupVCF), 'r', compressed=True) if args.dupVCF.endswith('.gz') else vcf.Reader(open(args.dupVCF), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outDUP, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if record.ID in set(reduce(operator.add, map(list, idPairs.keys()))):
            ci = max(abs(record.INFO['CIEND'][0]), abs(record.INFO['CIEND'][1]), abs(record.INFO['CIPOS'][0]), abs(record.INFO['CIPOS'][1]))
            gtl = list()
            rc = dict()
            nonRefHap = dict()
            for call in record.samples:
                if call.called:
                    rc[call.sample] = call['RC']
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        nonRefHap[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                gtl.append(call['GT'])
            selectedSVs[record.ID] = {'chr': record.CHROM, 'start': record.POS, 'end': record.INFO['END'], 'pe': record.INFO['PE'], 'mapq': record.INFO['MAPQ'], 'ci': ci, 'gt': gtl, 'rc': rc, 'nonRefHap': nonRefHap}
            vcf_writer.write_record(record)

# Create complex SV VCF file
if args.complexVCF:
    vcf_reader = vcf.Reader(open(args.delVCF), 'r', compressed=True) if args.delVCF.endswith('.gz') else vcf.Reader(open(args.delVCF), 'r', compressed=False)
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
            updSVType = validRdRatio(recO/nestedO, rdRatio, args.readDepth)[1]
            if updSVType != "DUP":
                continue
            if abs(s1-s2) > abs(e1-e2):
                (dupStart, dupEnd, insStart, insEnd) = (min(s1, s2), max(s1, s2), min(e1, e2), max(e1, e2))
            else:
                (insStart, insEnd, dupStart, dupEnd) = (min(s1, s2), max(s1, s2), min(e1, e2), max(e1, e2))
            svType = "PDUP"
            sv2info = {'start': insStart, 'end': insEnd, 'type': "INS"}
            print(selectedSVs[id1]['chr'], dupStart, dupEnd, updSVType, nestedO, recO, recO/nestedO, rdRatio, sep="\t")

            # output VCF
            info = "IMPRECISE;"
            info += "SVTYPE=" + updSVType + ";"
            info += "PE=" + str(selectedSVs[id1]['pe'] + selectedSVs[id2]['pe']) + ";"
            info += "MAPQ=" + str((selectedSVs[id1]['mapq'] + selectedSVs[id2]['mapq'])/2) + ";"
            info += "CIEND=" + str(-ci) + "," + str(ci) + ";CIPOS=" + str(-ci) + "," + str(ci) + ";SVMETHOD=" + svm + ";"
            info += "END=" + str(dupEnd) + ";SVLEN=" + str(dupEnd-dupStart+1) + ";"
            info += "MC=" + id1 + "," + id2
            if sv2info is not None:
                info += ";SV2INFO=" + selectedSVs[id2]['chr'] + "," + str(sv2info['start']) + "," + str(sv2info['end']) + "," + sv2info['type']
            print(selectedSVs[id1]['chr'], dupStart, updSVType + str(idCount).zfill(6), "N", "<" + svType + ">", ".", "PASS", info, "GT", sep="\t", file=f, end="")
            for (gt1, gt2) in zip(selectedSVs[id1]['gt'], selectedSVs[id2]['gt']):
                if (gt1 != None) and (gt1 == gt2):
                    print("\t" + gt1, end="", file=f)
                else:
                    print("\t./.", end="", file=f)
            print("", file=f)
            idCount += 1
