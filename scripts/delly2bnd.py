#! /usr/bin/env python

from __future__ import print_function
from .readfq import readfq
import argparse
import sys
import cyvcf2
import collections

def argument_parsing():
    """
    Parse command-line options.
    """
    parser = argparse.ArgumentParser(description='Split BND calls')
    parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF/BCF file (required)')
    parser.add_argument('-r', '--ref', metavar='ref.fa', required=True, dest='ref', help='input reference (required)')
    parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='out', help='output VCF file (required)')
    return parser.parse_args()

def delly2bnd(args):
    # Fetch breakpoint positions
    bndPos = collections.defaultdict(dict)
    oldChr = None
    if not len(bndPos):
        vcf = cyvcf2.VCF(args.vcf)
        for record in vcf:
            if record.CHROM != oldChr:
                oldChr = record.CHROM
                print("Fetching BND calls...", oldChr, file=sys.stderr)

            # Ignore multi-allelics
            if len(record.ALT) > 1:
                continue

            # Only delly BND and TRA calls
            if record.INFO.get("SVTYPE") == "BND":
                ct = record.INFO.get("CT")
                chrom2 = record.INFO.get("CHR2")
                pos2 = record.INFO.get("POS2")
            elif record.INFO.get("SVTYPE") == "TRA":
                ct = record.INFO.get("CT")
                chrom2 = record.INFO.get("CHR2")
                pos2 = record.INFO.get("END")
            else:
                continue
            chrom = record.CHROM
            pos = record.POS
            bndPos[chrom][pos] = 'N'
            bndPos[chrom2][pos2] = 'N'

        # Fetch nucleotides
        if True:
            f_in = gzip.open(args.ref) if args.ref.endswith('.gz') else open(args.ref)
            for seqName, seqNuc, seqQuals in readfq(f_in):
                if seqName in bndPos.keys():
                    print("Fetching breakpoint nucleotides...", seqName, file=sys.stderr)
                    for pos in bndPos[seqName].keys():
                        bndPos[seqName][pos] = seqNuc[(pos-1):pos]

    # Parse VCF/BCF
    vcf = cyvcf2.VCF(args.vcf)
    vcf.add_info_to_header({'ID': 'MATEID', 'Description': 'ID of mate breakends', 'Type':'String', 'Number': '.'})
    w = cyvcf2.Writer(args.out, vcf)
    oldChr = None
    for record in vcf:
        if record.CHROM != oldChr:
            oldChr = record.CHROM
            print("Processing...", oldChr, file=sys.stderr)

        # Ignore multi-allelics
        if len(record.ALT) > 1:
            continue

        # Only delly BND and TRA calls
        if record.INFO.get("SVTYPE") == "BND":
            ct = record.INFO.get("CT")
            chrom2 = record.INFO.get("CHR2")
            pos2 = record.INFO.get("POS2")
        elif record.INFO.get("SVTYPE") == "TRA":
            ct = record.INFO.get("CT")
            chrom2 = record.INFO.get("CHR2")
            pos2 = record.INFO.get("END")
        else:
            continue
        chrom = record.CHROM
        pos = record.POS
        id1 = record.ID + "_2nd"
        id2 = record.ID + "_1st"
        if ct == '3to5':
            template = bndPos[chrom][pos] + "[" + chrom2 + ":" + str(pos2) + "["
            template2 = "]" + chrom + ":" + str(pos) + "]" + bndPos[chrom2][pos2]
        elif ct == '5to3':
            template = "]" + chrom2 + ":" + str(pos2) + "]" + bndPos[chrom][pos]
            template2 = bndPos[chrom2][pos2] + "[" + chrom + ":" + str(pos) + "["
        elif ct == '3to3':
            template = bndPos[chrom][pos] + "]" + chrom2 + ":" + str(pos2) + "]"
            template2 = bndPos[chrom2][pos2] + "]" + chrom + ":" + str(pos) + "]"
        else:
            template = "[" + chrom2 + ":" + str(pos2) + "[" + bndPos[chrom][pos]
            template2 = "[" + chrom + ":" + str(pos) + "[" + bndPos[chrom2][pos2]

        # 2nd breakend
        record.ID = id1
        record.INFO['MATEID'] = id2
        record.ALT = [template]
        w.write_record(record)

        # 1st breakend
        record.CHROM = chrom2
        record.set_pos(pos2 - 1)
        record.ID = id2    
        record.INFO['MATEID'] = id1
        record.REF = bndPos[chrom2][pos2]
        record.ALT = [template2]
        record.INFO['CHR2'] = chrom
        record.INFO['POS2'] = pos
        record.INFO['END'] = pos2 + 1
        if ct == '5to3':
            record.INFO['CT'] = '3to5'
        elif ct == '3to5':
            record.INFO['CT'] = '5to3'
        w.write_record(record)
    # Close file
    w.close()


def main():
    arguments = argument_parsing()
    delly2bnd(arguments)


if __name__ == '__main__':
    main()
