#! /usr/bin/env python

from __future__ import print_function
from varpkg.readfq import readfq
import vcf
import argparse
import collections
import gzip

# Parse command line
parser = argparse.ArgumentParser(description='Add reference allele to VCF.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-r', '--ref', metavar='ref.fa', required=True, dest='ref', help='input reference (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
args = parser.parse_args()

# Get all the reference positions
selectedPos = collections.defaultdict(set)
nuclDict = collections.defaultdict(set)
if args.vcfFile:
    vcf_reader = gzip.open(args.vcfFile, 'rb') if args.vcfFile.endswith('.gz') else open(args.vcfFile, 'r')
    for line in vcf_reader:
        if line.startswith('#'):
            continue
        fields = line.split('\t', 2)
        selectedPos[fields[0]].add(int(fields[1]))
        nuclDict[(fields[0], int(fields[1]))] = 'N'
    vcf_reader.close()

# Store the true reference nucleotide
f_in = gzip.open(args.ref) if args.ref.endswith('.gz') else open(args.ref)
for seqName, seqNuc, seqQuals in readfq(f_in):
    for pos in selectedPos[seqName]:
        nuclDict[(seqName, pos)] = seqNuc[(pos-1):pos]

# Replace vcf reference allele
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        record.REF = nuclDict[(record.CHROM, record.POS)]
        vcf_writer.write_record(record)
