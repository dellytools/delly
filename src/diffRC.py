#! /usr/bin/env python

from __future__ import print_function
import scipy.stats
import vcf
import argparse
import csv
import sys
import collections
import numpy

# Parse command line
parser=argparse.ArgumentParser(description='Annotates sites for differential read counts between carrier and non-carrier samples.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='outFile', help='output vcf file (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svType', help='SV type [DEL, DUP, INV] (required)')
parser.add_argument('-p', '--pvalue', metavar='1.0', required=False, dest='pValueCut', help='p-value threshold for read-depth')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Set parameters
pValueCut=1.0
if (args.pValueCut):
    pValueCut=float(args.pValueCut)

# Parse vcf File
if args.vcfFile:
    vcf_reader=vcf.Reader(open(args.vcfFile), 'r')
    vcf_reader.infos['RDpval'] = vcf.parser._Info('RDpval', 1, 'Float', 'RD p-value estimated from carrier vs. non-carrier samples.')
    vcf_writer = vcf.Writer(open(args.outFile, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        rcRef=[]
        rcAlt=[]
        for call in record.samples:
            if ((call.called) and (call['FT']=="PASS") and (call.gt_type==0)):
                rcRef.append(call['RC'])
            if ((call.called) and (call.gt_type!=0)):
                rcAlt.append(call['RC'])
        wilcoxonPval=1.0
        if ((len(rcRef)>0) and (len(rcAlt)>0)):
            wilcoxonStat=scipy.stats.ranksums(rcRef, rcAlt)
            if ((args.svType=="DEL") and (wilcoxonStat[0]>0)):
                wilcoxonPval=wilcoxonStat[1]
            elif ((args.svType=="DUP") and (wilcoxonStat[0]<0)):
                wilcoxonPval=wilcoxonStat[1]
            else:
                wilcoxonPval=wilcoxonStat[1]
        record.INFO['RDpval']=float(round(wilcoxonPval,6))
        if ((record.INFO['RDpval'] <= pValueCut) and ((not args.siteFilter) or (len(record.FILTER)==0))):
            vcf_writer.write_record(record)
