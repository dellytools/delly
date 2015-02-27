#! /usr/bin/env python

from __future__ import print_function
from varpkg.readfq import readfq
import vcf
import argparse
import gzip
import banyan
import re

# Parse command line
parser = argparse.ArgumentParser(description='Add read-depth control region to input VCF.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-r', '--ref', metavar='ref.fa', required=True, dest='ref', help='input reference (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
args = parser.parse_args()

# Compute all stretches of Ns in the reference
nRun = dict()
refLen = dict()
f_in = gzip.open(args.ref) if args.ref.endswith('.gz') else open(args.ref)
for seqName, seqNuc, seqQuals in readfq(f_in):
    refLen[seqName] = len(seqNuc)
    print("Processing", seqName, refLen[seqName])
    if not nRun.has_key(seqName):
        nRun[seqName] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
    for m in re.compile("([Nn]+)").finditer(seqNuc):
        nRun[seqName][m.span()] = 1

# Add read-depth control region to VCF file
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    vcf_reader.infos['CONTROL'] = vcf.parser._Info('CONTROL', 0, 'Flag', 'Control variant.')
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        svSize = record.INFO['END']-record.POS
        svControlStart = max(record.POS - svSize/2, 0)
        svControlEnd = min(record.INFO['END'] + svSize/2, refLen[record.CHROM])
        if len(nRun[record.CHROM].overlap((svControlStart, svControlEnd))):
            svControlStart = max(record.POS - svSize, 0)
            svControlEnd = record.INFO['END']
            if len(nRun[record.CHROM].overlap((svControlStart, svControlEnd))):
                svControlStart = record.POS
                svControlEnd = min(record.INFO['END'] + svSize, refLen[record.CHROM])
        record.INFO['CONTROL'] = False
        vcf_writer.write_record(record)
        record.ID = record.ID.replace('0', '9', 1)
        record.POS = svControlStart
        record.INFO['END'] = svControlEnd
        record.INFO['CONTROL'] = True
        vcf_writer.write_record(record)
