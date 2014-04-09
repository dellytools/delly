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
parser.add_argument('-g', '--gqAlt', metavar='15', required=False, dest='gqAlt', help='min. GQ for carriers (optional)')
parser.add_argument('-q', '--gqRef', metavar='15', required=False, dest='gqRef', help='min. GQ for non-carriers (optional)')
parser.add_argument('-m', '--minsize', metavar='500', required=False, dest='minSize', help='min. size (optional)')
parser.add_argument('-n', '--maxsize', metavar='5000000', required=False, dest='maxSize', help='max. size (optional)')
parser.add_argument('-a', '--altaf', metavar='0.4', required=False, dest='altAF', help='min. alt. AF (optional)')
parser.add_argument('-r', '--ratioGeno', metavar='0.4', required=False, dest='ratioGeno', help='min. fraction of genotyped samples (optional)')
parser.add_argument('-s', '--sample', metavar='NA12878', required=False, dest='sampleID', help='required carrier sample (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
parser.add_argument('-p', '--paired', dest='pairedFilter', action='store_true', help='Require 3to3 and 5to5 inversion support')
args = parser.parse_args()

# Command-line args
sampleID = ""
if args.sampleID:
    sampleID = args.sampleID
gqAltCut = 15
if args.gqAlt:
    gqAltCut = int(args.gqAlt)
gqRefCut = 15
if args.gqRef:
    gqRefCut = int(args.gqRef)
minSize = 500
if args.minSize:
    minSize = int(args.minSize)
maxSize = 5000000
if args.maxSize:
    maxSize = int(args.maxSize)
altAF = 0.4
if args.altAF:
    altAF = float(args.altAF)
ratioGeno = 0.4
if args.ratioGeno:
    ratioGeno = float(args.ratioGeno)

# Collect high-quality SVs
sv = dict()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    for record in vcf_reader:
        if (record.INFO['SVLEN'] >= minSize) and (record.INFO['SVLEN'] <= maxSize) and ((not args.siteFilter) or (len(record.FILTER) == 0)):
            gqRef = []
            gqAlt = []
            ratioRef = [0]
            ratioAlt = []
            carrierSample = False
            if sampleID == "":
                carrierSample = True
            for call in record.samples:
                if (call.called):
                    if call.gt_type == 0:
                        ratioRef.append(float(call['DV'])/float(call['DR'] + call['DV']))
                        if (call['DV'] == 0):
                            gqRef.append(call['GQ'])
                    if call.gt_type != 0:
                        ratioAlt.append(float(call['DV'])/float(call['DR'] + call['DV']))
                        if (not carrierSample) and (call.sample == sampleID):
                            carrierSample = True
                        if (call['DV'] >= 2):
                            gqAlt.append(call['GQ'])
            genotypeRatio = float(len(gqAlt)+len(gqRef)) /  float(len(record.samples))
            if (carrierSample) and (genotypeRatio>ratioGeno):
                if (len(gqRef)) and (len(gqAlt)) and (numpy.median(gqRef) >= gqRefCut) and (numpy.median(gqAlt) >= gqAltCut):
                    if (numpy.percentile(ratioRef, 99) == 0) and (numpy.median(ratioAlt) >= altAF):
                        print(record.CHROM, record.POS, record.INFO['END'], record.INFO['CT'])
                        if not sv.has_key(record.CHROM):
                            sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                        sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'], record.INFO['CT'])


# Kick-out the unpaired SVs
if (args.pairedFilter):
    filteredSVs = dict()
    for chrName in sv.keys():
        for start, end in sv[chrName].keys():
            overlapList = sv[chrName].overlap((start, end))
            svID, score, ct = sv[chrName][(start, end)]
            for cStart, cEnd in overlapList:
                cSvID, cScore, cCt = sv[chrName][(cStart, cEnd)]
                if (svID!=cSvID) and (ct!=cCt) and (overlapValid(start, end, cStart, cEnd, 0.9, 100)):
                    if not filteredSVs.has_key(chrName):
                        filteredSVs[chrName] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    filteredSVs[chrName][(start, end)] = sv[chrName][(start, end)]
                    break
    sv=filteredSVs

# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if (record.CHROM not in sv.keys()) or ((record.POS, record.INFO['END']) not in sv[record.CHROM].keys()):
            continue
        overlapList = sv[record.CHROM].overlap((record.POS, record.INFO['END']))
        svID, score, ct = sv[record.CHROM][(record.POS, record.INFO['END'])]
        foundBetterHit = False
        for cStart, cEnd in overlapList:
            cSvID, cScore, cCt = sv[record.CHROM][(cStart, cEnd)]
            if svID == cSvID:
                continue
            # There should be at least 10% overlap otherwise ignore it
            if not overlapValid(record.POS, record.INFO['END'], cStart, cEnd, 0.1, 10000000):
                continue
            if (cScore > score) or ((cScore == score) and (cSvID < svID)):
                foundBetterHit = True
                break
        if not foundBetterHit:
            vcf_writer.write_record(record)
