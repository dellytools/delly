#! /usr/bin/env python

from __future__ import print_function
import argparse
import vcf
import numpy
import re
import banyan
import collections
import copy

#Functions
def overlapValid(s1, e1, s2, e2, reciprocalOverlap=0.8, maxOffset=250):
    if (e1 < s2) or (s1 > e2):
        return False
    overlapLen = float(min(e1, e2) - max(s1, s2))
    lenA=float(e1-s1)
    lenB=float(e2-s2)
    if (overlapLen<0):
        sys.exit("Sites are not overlapping.")
    if (lenA<=0) or (lenB<=0):
        sys.exit("Invalid intervals.")
    # Check reciprocal overlap
    if (overlapLen/max(lenA,lenB))<reciprocalOverlap:
        return False
    # Check offset
    if max(abs(s2-s1), abs(e2-e1))>maxOffset:
        return False
    return True


# Parse command line
parser = argparse.ArgumentParser(description='Filter for somatic SVs.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svType', help='SV type [DEL, DUP, INV, INS] (required)')
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
svDups = collections.defaultdict(list)
validRecordID = set()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    for record in vcf_reader:
        if (record.INFO['SVLEN'] >= minSize) and (record.INFO['SVLEN'] <= maxSize) and ((not args.siteFilter) or (len(record.FILTER) == 0)):
            precise = False
            if ('PRECISE' in record.INFO.keys()):
                precise = record.INFO['PRECISE']
            rcRef = []
            rcAlt = []
            for call in record.samples:
                if (call.called):
                    if (re.search(r"[Nn]ormal", call.sample) != None) and (call.gt_type == 0):
                        if ((not precise) and (call['DV'] == 0)) or ((precise) and (call['RV']==0)):
                            rcRef.append(call['RC'])
                    if (re.search(r"[Tt]umor", call.sample) != None) and (call.gt_type != 0):
                        if ((not precise) and (call['DV'] >= 2) and (float(call['DV'])/float(call['DV']+call['DR'])>=altAF)) or ((precise) and (call['RV'] >= 2) and (float(call['RV'])/float(call['RR'] + call['RV'])>=altAF)):
                            rcAlt.append(call['RC'])
            if (len(rcRef) > 0) and (len(rcAlt) > 0):
                rdRatio = 1
                if numpy.median(rcRef):
                    rdRatio = numpy.median(rcAlt)/numpy.median(rcRef)
                if (args.svType == 'INV') or (args.svType == 'INS') or (record.INFO['SVLEN'] <= 10000) or ((args.svType == 'DEL') and (rdRatio <= 0.85)) or ((args.svType == 'DUP') and (rdRatio >= 1.15)):
                    validRecordID.add(record.ID)
                    if not sv.has_key(record.CHROM):
                        sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                    if (record.POS, record.INFO['END']) in sv[record.CHROM]:
                        svDups[(record.CHROM, record.POS, record.INFO['END'])].append((record.ID, record.INFO['PE']))
                    else:
                        sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'])

# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r')
    vcf_reader.infos['SOMATIC'] = vcf.parser._Info('SOMATIC', 0, 'Flag', 'Somatic structural variant.')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        # Is it a valid SV?
        if (record.ID not in validRecordID):
            continue
        # Collect overlapping and identical calls
        overlapCalls = list()
        for cSvID, cScore in svDups[(record.CHROM, record.POS, record.INFO['END'])]:
            if (record.ID != cSvID):
                overlapCalls.append((cSvID, cScore))
        for cStart, cEnd in sv[record.CHROM].overlap((record.POS, record.INFO['END'])):
            cSvID, cScore = sv[record.CHROM][(cStart, cEnd)]
            if (record.ID != cSvID) and (overlapValid(record.POS, record.INFO['END'], cStart, cEnd)):
                overlapCalls.append((cSvID, cScore))
        # Judge wether overlapping calls are better
        foundBetterHit = False
        for cSvID, cScore in overlapCalls:
            if (cScore > record.INFO['PE']) or ((cScore == record.INFO['PE']) and (cSvID < record.ID)):
                foundBetterHit = True
                break
        if not foundBetterHit:
            record.INFO['SOMATIC'] = True
            vcf_writer.write_record(record)
