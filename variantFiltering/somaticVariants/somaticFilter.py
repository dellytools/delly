#! /usr/bin/env python

from __future__ import print_function
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from varpkg.overlap import overlapValid
import argparse
import vcf
import re
import banyan
import collections
import numpy


# Parse command line
parser = argparse.ArgumentParser(description='Filter for somatic SVs.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svType', help='SV type [DEL, DUP, INV, INS, TRA] (required)')
parser.add_argument('-a', '--altaf', metavar='0.2', required=False, dest='altAF', help='min. alt. AF (optional)')
parser.add_argument('-c', '--mincov', metavar='10', required=False, dest='minCov', help='min. coverage (optional)')
parser.add_argument('-m', '--minsize', metavar='500', required=False, dest='minSize', help='min. size (optional)')
parser.add_argument('-n', '--maxsize', metavar='500000000', required=False, dest='maxSize', help='max. size (optional)')
parser.add_argument('-r', '--ratioGeno', metavar='0.75', required=False, dest='ratioGeno', help='min. fraction of genotyped samples (optional)')
parser.add_argument('-i', '--normalContamination', metavar='0.0', required=False, dest='normalCont', help='normal contamination (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
parser.add_argument('-N', '--normal', metavar='normalID', required=False, dest='nameNormal', help='normal sample name as in VCF (optional)')
parser.add_argument('-T', '--tumor', metavar='tumorID', required=False, dest='nameTumor', help='tumor sample name as in VCF (optional)')

args = parser.parse_args()

# Command-line args
minSize = 500
if args.minSize:
    minSize = int(args.minSize)
maxSize = 500000000
if args.maxSize:
    maxSize = int(args.maxSize)
altAF = 0.2
if args.altAF:
    altAF = float(args.altAF)
minCov = 10
if args.minCov:
    minCov = int(args.minCov)
ratioGeno = 0.75
if args.ratioGeno:
    ratioGeno = float(args.ratioGeno)
normalCont = 0.0
if args.normalCont:
    normalCont = float(args.normalCont)
traWindow = 2500  # 2.5kb translocation window


# Collect high-quality SVs
sv = dict()
rdRat = dict()
svDups = collections.defaultdict(list)
validRecordID = set()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    for record in vcf_reader:
        svlen = record.INFO['END'] - record.POS
        if ((args.svType == 'TRA') or ((svlen >= minSize) and (svlen <= maxSize))) and ((not args.siteFilter) or (len(record.FILTER) == 0)):
            precise = False
            if 'PRECISE' in record.INFO.keys():
                precise = record.INFO['PRECISE']
            rcRef = []
            rcAlt = []
            nCount = 0
            tCount = 0
            refpass = False
            altpass = False
            for call in record.samples:
                if call.called:
                    if (args.nameNormal and call.sample == args.nameNormal) or ((args.nameNormal is None) and ((re.search(r"[Nn]ormal", call.sample) != None) or (re.search(r"[Cc]ontrol", call.sample) != None))):
                        nCount += 1
                        if call.gt_type == 0:
                            if not precise:
                                cov = float(call['DV']+call['DR'])
                                if float(call['DV'])/cov <= normalCont:
                                    rcRef.append(call['RC'])
                                    if call['DR'] >= minCov:
                                        refpass = True
                            else:
                                cov = float(call['RV']+call['RR'])
                                if float(call['RV'])/cov <= normalCont:
                                    rcRef.append(call['RC'])
                                    if call['RR'] >= minCov:
                                        refpass = True
                    if (args.nameTumor and call.sample == args.nameTumor) or ((args.nameTumor is None) and (re.search(r"[Tt]umo[ur]", call.sample) != None)):
                        tCount += 1
                        if not precise:
                            cov = float(call['DV']+call['DR'])
                            if (call['DV'] >= 2) and (float(call['DV'])/cov >= altAF):
                                if (record.INFO['SVTYPE'] != "TRA") or (call['DV'] >= 5):
                                    rcAlt.append(call['RC'])
                                    if cov >= minCov:
                                        altpass = True
                        else:
                            cov = float(call['RV']+call['RR'])
                            if (call['RV'] >= 2) and (float(call['RV'])/cov >= altAF):
                                if (record.INFO['SVTYPE'] != "TRA") or (call['RV'] >= 5):
                                    rcAlt.append(call['RC'])
                                    if cov >= minCov:
                                        altpass = True
            genotypeRatio = float(nCount + tCount) /  float(len(record.samples))
            if (refpass) and (altpass) and (len(rcRef) == nCount) and (len(rcAlt) > 0) and (genotypeRatio >= ratioGeno):
                rdRatio = 1
                if numpy.median(rcRef):
                    rdRatio = round(numpy.median(rcAlt)/numpy.median(rcRef), 4)
                rdRat[record.ID] = rdRatio
                validRecordID.add(record.ID)
                if not sv.has_key(record.CHROM):
                    sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                if not sv.has_key(record.INFO['CHR2']):
                    sv[record.INFO['CHR2']] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
                if args.svType == 'TRA':
                    if (record.POS - traWindow, record.POS + traWindow) in sv[record.CHROM]:
                        svDups[(record.CHROM, record.POS - traWindow, record.POS + traWindow)].append((record.ID, record.INFO['PE']))
                    else:
                        sv[record.CHROM][(record.POS - traWindow, record.POS + traWindow)] = (record.ID, record.INFO['PE'])
                    if (record.INFO['END'] - traWindow, record.INFO['END'] + traWindow) in sv[record.INFO['CHR2']]:
                        svDups[(record.INFO['CHR2'], record.INFO['END'] - traWindow, record.INFO['END'] + traWindow)].append((record.ID, record.INFO['PE']))
                    else:
                        sv[record.INFO['CHR2']][(record.INFO['END'] - traWindow, record.INFO['END'] + traWindow)] = (record.ID, record.INFO['PE'])
                else:
                    if (record.POS, record.INFO['END']) in sv[record.CHROM]:
                        svDups[(record.CHROM, record.POS, record.INFO['END'])].append((record.ID, record.INFO['PE']))
                    else:
                        sv[record.CHROM][(record.POS, record.INFO['END'])] = (record.ID, record.INFO['PE'])



# Output vcf records
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    vcf_reader.infos['RDRATIO'] = vcf.parser._Info('RDRATIO', 1, 'Float', 'Read-depth ratio of tumor vs. normal.')
    vcf_reader.infos['SOMATIC'] = vcf.parser._Info('SOMATIC', 0, 'Flag', 'Somatic structural variant.')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        # Is it a valid SV?
        if record.ID not in validRecordID:
            continue
        # Judge wether overlapping calls are better
        foundBetterHit = False
        if args.svType == 'TRA':
            countBetterHits = 0
            for (chrName, start, end) in [(record.CHROM, record.POS - traWindow, record.POS + traWindow), (record.INFO['CHR2'], record.INFO['END'] - traWindow, record.INFO['END'] + traWindow)]:
                for cStart, cEnd in sv[chrName].overlap((start, end)):
                    for cSvID, cScore in svDups[(chrName, cStart, cEnd)] + [sv[chrName][(cStart, cEnd)]]:
                        if record.ID != cSvID:
                            if (cScore > record.INFO['PE']) or ((cScore == record.INFO['PE']) and (cSvID < record.ID)):
                                countBetterHits += 1
            if countBetterHits > 2:
                foundBetterHit = True
        else:
            for cStart, cEnd in sv[record.CHROM].overlap((record.POS, record.INFO['END'])):
                for cSvID, cScore in svDups[(record.CHROM, cStart, cEnd)] + [sv[record.CHROM][(cStart, cEnd)]]:
                    if (record.ID != cSvID) and (overlapValid((record.POS, record.INFO['END']), (cStart, cEnd))):
                        if (cScore > record.INFO['PE']) or ((cScore == record.INFO['PE']) and (cSvID < record.ID)):
                            foundBetterHit = True
                            break
                if foundBetterHit:
                    break
        if not foundBetterHit:
            record.INFO['RDRATIO'] = rdRat[record.ID]
            record.INFO['SOMATIC'] = True
            vcf_writer.write_record(record)
