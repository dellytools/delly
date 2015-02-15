#! /usr/bin/env python

from __future__ import print_function
import vcf
import argparse
import numpy
import banyan
import collections


#Functions
def overlapValid((s1, e1), (s2, e2), reciprocalOverlap=0.5, maxOffset=500):
    if (e1 < s2) or (s1 > e2) or ((e1-s1) <= 0) or ((e2-s2) <= 0):
        return False
    overlapLen = float(min(e1, e2) - max(s1, s2))
    # Check reciprocal overlap
    if (overlapLen < 1) or (overlapLen/float(e1-s1) < reciprocalOverlap) or (overlapLen/float(e2-s2) < reciprocalOverlap):
        return False
    # Check offset
    if max(abs(s2-s1), abs(e2-e1)) > maxOffset:
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
parser.add_argument('-k', '--keepOverlap', dest='keepOverlap', action='store_true', help='keep overlapping PE calls')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Command-line args
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
svDups = collections.defaultdict(list)
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    for record in vcf_reader:
        if (record.INFO['SVLEN'] >= minSize) and (record.INFO['SVLEN'] <= maxSize) and ((not args.siteFilter) or (len(record.FILTER) == 0)):
            precise = False
            if 'PRECISE' in record.INFO.keys():
                precise = record.INFO['PRECISE']
            gqRef = []
            gqAlt = []
            ratioRef = [0]
            ratioAlt = []
            for call in record.samples:
                if call.called:
                    if call.gt_type == 0:
                        if precise:
                            ratioRef.append(float(call['RV'])/float(call['RR'] + call['RV']))
                        else:
                            ratioRef.append(float(call['DV'])/float(call['DR'] + call['DV']))
                        if ((precise) and (call['RV'] == 0)) or ((not precise) and (call['DV'] == 0)):
                            gqRef.append(call['GQ'])
                    if (call.gt_type != 0) and (call['FT']=="PASS"):
                        if precise:
                            ratioAlt.append(float(call['RV'])/float(call['RR'] + call['RV']))
                        else:
                            ratioAlt.append(float(call['DV'])/float(call['DR'] + call['DV']))
                        if ((precise) and (call['RV'] >= 2)) or ((not precise) and (call['DV'] >= 2)):
                            gqAlt.append(call['GQ'])
            genotypeRatio = float(len(gqAlt)+len(gqRef)) /  float(len(record.samples))
            if genotypeRatio > ratioGeno:
                if (len(gqRef)) and (len(gqAlt)) and (numpy.median(gqRef) >= gqRefCut) and (numpy.median(gqAlt) >= gqAltCut):
                    if (numpy.percentile(ratioRef, 99) == 0) and (numpy.median(ratioAlt) >= altAF):
                        #print(record.INFO['END']-record.POS, len(gqRef), len(gqAlt), numpy.median(gqRef), numpy.median(gqAlt), numpy.percentile(ratioRef, 99), numpy.median(ratioAlt), genotypeRatio, sep="\t")
                        if not sv.has_key(record.CHROM):
                            sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                        if (record.POS, record.INFO['END']) not in sv[record.CHROM]:
                            sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'], record.INFO['CT'])
                        else:
                            svDups[(record.CHROM, record.POS, record.INFO['END'])].append((record.ID, record.INFO['PE'], record.INFO['CT']))

# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        if (record.CHROM not in sv.keys()) or ((record.POS, record.INFO['END']) not in sv[record.CHROM].keys()):
            continue
        overlapList = sv[record.CHROM].overlap((record.POS, record.INFO['END']))
        foundBetterHit = False
        if not args.keepOverlap:
            for cStart, cEnd in sv[record.CHROM].overlap((record.POS, record.INFO['END'])):
                if foundBetterHit:
                    break
                for cSvID, cScore, cCt in svDups[(record.CHROM, cStart, cEnd)] + [sv[record.CHROM][(cStart, cEnd)]]:
                    if record.ID != cSvID:
                        if (cScore > record.INFO['PE']) or ((cScore == record.INFO['PE']) and (cSvID < record.ID)):
                            if overlapValid((record.POS, record.INFO['END']), (cStart, cEnd), 0.1, 10000000):
                                foundBetterHit = True
                                break

        # Output VCF record
        if not foundBetterHit:
            vcf_writer.write_record(record)
