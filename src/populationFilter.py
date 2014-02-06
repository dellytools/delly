#! /usr/bin/env python

from __future__ import print_function
import vcf
import argparse
import numpy
import banyan


#Functions
def overlapValid(s1, e1, s2, e2, reciprocalOverlap=0.5, maxOffset=500):
    if (e1 < s2) or (s1 > e2):
        return False
    overlapLen = float(min(e1, e2) - max(s1, s2))
    # Check reciprocal overlap
    if (overlapLen < 1) or (float(e1-s1)/overlapLen < reciprocalOverlap) or (float(e2-s2)/overlapLen < reciprocalOverlap):
        return False
    # Check offset
    if (abs(s2-s1) > maxOffset) or (abs(e2-e1) > maxOffset):
        return False
    return True


# Parse command line
parser = argparse.ArgumentParser(description='Filter for reliable SV sites.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-g', '--gq', metavar='30', required=False, dest='gqCut', help='min. GQ (optional)')
parser.add_argument('-m', '--minsize', metavar='500', required=False, dest='minSize', help='min. size (optional)')
parser.add_argument('-n', '--maxsize', metavar='5000000', required=False, dest='maxSize', help='max. size (optional)')
parser.add_argument('-s', '--sample', metavar='NA12878', required=False, dest='sampleID', help='required carrier sample (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Command-line args
sampleID = ""
if args.sampleID:
    sampleID = args.sampleID
gqCut = 30
if args.gqCut:
    gqCut = int(args.gqCut)
minSize = 500
if args.minSize:
    minSize = int(args.minSize)
maxSize = 5000000
if args.maxSize:
    maxSize = int(args.maxSize)

# Collect high-quality SVs
sv = dict()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    for record in vcf_reader:
        if (not args.siteFilter) or (len(record.FILTER) == 0):
            gqRef = []
            gqAlt = []
            ratioRef = []
            ratioAlt = []
            carrierSample = False
            if sampleID == "":
                carrierSample = True
            for call in record.samples:
                if (call.called) and (call['FT'] == "PASS"):
                    if call.gt_type == 0:
                        if call['DV'] == 0:
                            gqRef.append(call['GQ'])
                        else:
                            ratioRef.append(float(call['DV'])/float(call['DR'] + call['DV']))
                    if call.gt_type != 0:
                        if (not carrierSample) and (call.sample == sampleID):
                            carrierSample = True
                        gqAlt.append(call['GQ'])
                        ratioAlt.append(float(call['DV'])/float(call['DR'] + call['DV']))
            if carrierSample:
                if (len(gqRef)) and (len(gqAlt)) and (record.INFO['SVLEN'] >= minSize) and (record.INFO['SVLEN'] <= maxSize) and (numpy.median(gqRef) >= gqCut) and (numpy.median(gqAlt) >= gqCut):
                    if (len(ratioRef) == 0) and (numpy.median(ratioAlt) >= 0.3):
                        if not sv.has_key(record.CHROM):
                            sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                        sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'])

# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if (record.CHROM not in sv.keys()) or ((record.POS, record.INFO['END']) not in sv[record.CHROM].keys()):
            continue
        overlapList = sv[record.CHROM].overlap((record.POS, record.INFO['END']))
        svID, score = sv[record.CHROM][(record.POS, record.INFO['END'])]
        foundBetterHit = False
        for cStart, cEnd in overlapList:
            cSvID, cScore = sv[record.CHROM][(cStart, cEnd)]
            if svID == cSvID:
                continue
            if not overlapValid(record.POS, record.INFO['END'], cStart, cEnd, 0.1, 10000000):
                continue
            if (cScore > score) or ((cScore == score) and (cSvID < svID)):
                foundBetterHit = True
                break
        if not foundBetterHit:
            vcf_writer.write_record(record)
