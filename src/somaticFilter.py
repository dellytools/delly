#! /usr/bin/env python

from __future__ import print_function
import argparse
import vcf
import numpy
import re
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
parser = argparse.ArgumentParser(description='Filter for somatic SVs.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svType', help='SV type [DEL, DUP, INV] (required)')
parser.add_argument('-a', '--altaf', metavar='0.25', required=False, dest='altAF', help='min. alt. AF (optional)')
parser.add_argument('-m', '--minsize', metavar='500', required=False, dest='minSize', help='min. size (optional)')
parser.add_argument('-n', '--maxsize', metavar='500000000', required=False, dest='maxSize', help='max. size (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Command-line args
minSize = 500
if args.minSize:
    minSize = int(args.minSize)
maxSize = 500000000
if args.maxSize:
    maxSize = int(args.maxSize)
altAF = 0.25
if args.altAF:
    altAF = float(args.altAF)

# Collect high-quality SVs
sv = dict()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    for record in vcf_reader:
        rcRef = []
        rcAlt = []
        for call in record.samples:
            if (re.search(r"[Nn]ormal", call.sample) != None) and (call.called) and (call.gt_type == 0) and (call['DV'] == 0):
                rcRef.append(call['RC'])
            if (re.search(r"[Tt]umor", call.sample) != None) and (call.called) and (call.gt_type != 0) and (call['DV'] >= 2) and (float(call['DV'])/float(call['DV']+call['DR'])>=altAF):
                rcAlt.append(call['RC'])
        if (len(rcRef) > 0) and (len(rcAlt) > 0) and (record.INFO['SVLEN'] >= minSize) and (record.INFO['SVLEN'] <= maxSize):
            rdRatio = numpy.median(rcAlt)/numpy.median(rcRef)
            if (args.svType == 'INV') or (record.INFO['SVLEN'] <= 10000) or ((args.svType == 'DEL') and (rdRatio <= 0.75)) or ((args.svType == 'DUP') and (rdRatio >= 1.25)):
                if (not args.siteFilter) or (len(record.FILTER) == 0):
                    if not sv.has_key(record.CHROM):
                        sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'])

# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    vcf_reader.infos['SOMATIC'] = vcf.parser._Info('SOMATIC', 0, 'Flag', 'Somatic structural variant.')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if (record.CHROM not in sv.keys()) or ((record.POS, record.INFO['END']) not in sv[record.CHROM].keys()):
            continue
        overlapList = sv[record.CHROM].overlap((record.POS, record.INFO['END']))
        svID, score = sv[record.CHROM][(record.POS, record.INFO['END'])]
        foundBetterHit = False
        for cStart, cEnd in overlapList:
            if not overlapValid(record.POS, record.INFO['END'], cStart, cEnd):
                continue
            cSvID, cScore = sv[record.CHROM][(cStart, cEnd)]
            if svID == cSvID:
                continue
            if (cScore > score) or ((cScore == score) and (cSvID < svID)):
                foundBetterHit = True
                break
        if not foundBetterHit:
            record.INFO['SOMATIC'] = True
            vcf_writer.write_record(record)
