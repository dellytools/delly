#! /usr/bin/env python

from __future__ import print_function
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
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
    if 'CONTROL' not in vcf_reader.infos.keys():
        vcf_reader.infos['CONTROL'] = vcf.parser._Info('CONTROL', 1, 'Integer', 'Control variant.')
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        svSize = record.INFO['END']-record.POS
        svControlID = re.sub(r"^[A-Z0]*","", record.ID)

        # Left control region
        svLeftID = record.ID.replace('0', '8', 1)
        svControlLeftStart = max(record.POS - svSize/2, 0)
        svControlLeftEnd = record.POS
        while len(nRun[record.CHROM].overlap((svControlLeftStart, svControlLeftEnd))):
            minN=refLen[record.CHROM]
            for (nStart, nEnd) in nRun[record.CHROM].overlap((svControlLeftStart, svControlLeftEnd)):
                if (nStart<minN):
                    minN=nStart
            minN=minN-1
            svControlLeftStart = max(minN - svSize/2, 0)
            svControlLeftEnd = minN

        # Right control region
        svRightID = record.ID.replace('0', '9', 1)
        svControlRightStart = record.INFO['END']
        svControlRightEnd = min(record.INFO['END'] + svSize/2, refLen[record.CHROM])
        while len(nRun[record.CHROM].overlap((svControlRightStart, svControlRightEnd))):
            maxN=0
            for (nStart, nEnd) in nRun[record.CHROM].overlap((svControlRightStart, svControlRightEnd)):
                if (nEnd>maxN):
                    maxN=nEnd
            maxN=maxN+1
            svControlRightStart = maxN
            svControlRightEnd = min(maxN + svSize/2, refLen[record.CHROM])

        # Output control regions
        if (svControlLeftStart>0) and (svControlRightEnd<refLen[record.CHROM]):
            record.INFO['CONTROL'] = 0
            vcf_writer.write_record(record)
            record.ID = svLeftID
            record.POS = svControlLeftStart
            record.INFO['END'] = svControlLeftEnd
            record.INFO['CONTROL'] = svControlID
            vcf_writer.write_record(record)

            record.ID = svRightID
            record.POS = svControlRightStart
            record.INFO['END'] = svControlRightEnd
            record.INFO['CONTROL'] = svControlID
            vcf_writer.write_record(record)

