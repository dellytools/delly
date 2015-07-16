#!/usr/bin/env python
# joachim.weischenfeldt@embl.de

# -*- coding: utf-8 -*-

"""
Usage:
    dellyVcf2BndVcf.py -v <vcfFile> -o <outVCF> -e <emailaddr>

Options:
  -h --help     Show this screen.
  -v --vcfFile  DELLY v0.5+
  -o --output   vcf file
  -e --emailaddr   your email address
  

Convert DELLY VCF to BREAKEND (BND) VCF format. Splits each DELLY vcf pair of a call on two lines, add nucleotide at ALT col from hg19

"""

from __future__ import print_function
from docopt import docopt
import vcf
import os,sys,csv,re,copy,urllib2,gzip
from itertools import islice
from Bio import Entrez, SeqIO


arguments = docopt(__doc__)

vcfFile = arguments['<vcfFile>']
outVCF = arguments['<outVCF>']
emailaddr = arguments['<emailaddr>']

if not re.match(r'[^@]+@[^@]+\.[^@]+', emailaddr):
    print ("\nyou need to specify a valid email address ('yourname@host.com') to make requests to ncbi ftp server. Will exit now\n")
    sys.exit(1)


# Get NCBI hg19 genome ref chrom names
ftp_file = 'hg19_ftp.txt'
ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.13.assembly.txt'
req = urllib2.Request(ftp_url)
response = urllib2.urlopen(req)
ftp_file = response.read().split('\t')
chrom_dict = dict()
for e,i in enumerate(ftp_file):
    if i == "assembled-molecule":
        chrom_dict[ftp_file[e+1]] = ftp_file[e+5]


def retrieve_nt(chrom,start,end):
    '''
    retrieves the nucleotide(s), given a chromosome, start, end coordinate
    
    execute: retrieve_nt('1',1000000,1000001)
    
    '''
    chrom = chrom.replace('chr','').replace('.fa','')
    try:
        chrom_entrez = chrom_dict[str(chrom)]
    except Exception, e:
        print (chrom, 'not in chromomosome dict')
        return '.'
    Entrez.email = emailaddr
    handle = Entrez.efetch(db="nucleotide", 
                           id=chrom_entrez , 
                           rettype="fasta", 
                           strand=1, 
                           seq_start=start, 
                           seq_stop=end)
    entrez_record = SeqIO.read(handle, "fasta")
    handle.close()
    return entrez_record.seq[0]



def fix_offset(record, pos1, pos2):
    '''
    correct breakpoint offset error in breakpoint annotated DEL and DUP (DELLY V.0.0.6)
    '''
    out = pos1, pos2, ''
    if record.INFO.has_key('CONSENSUS'):
        if record.INFO['SVTYPE'] == 'DEL':
            pos2 = int(pos2) + 1
            out = pos1, pos2, abs(pos1 - pos2)
        if record.INFO['SVTYPE'] == 'DUP':
            pos1 = int(pos1) + 1
            out = pos1, pos2, abs(pos1 - pos2)
    return out


if vcfFile:
    vcf_reader = vcf.Reader(open(vcfFile), 'r', compressed=True) if vcfFile.endswith('.gz') else vcf.Reader(open(vcfFile), 'r', compressed=False)
    vcf_reader.infos['SVCLASS'] = vcf.parser._Info('SVCLASS', '1', 'String', 'Class of structural variant')
    vcf_reader.infos['STRAND'] = vcf.parser._Info('STRAND', '1', 'String', 'Breakend strand (+/-)')
    vcf_reader.infos['MATESTRAND'] = vcf.parser._Info('MATESTRAND', '1', 'String', 'Mate breakend strand (+/-)')
    vcf_reader.infos['MATEPOS'] = vcf.parser._Info('MATEPOS', '1', 'String', 'Mate position')
    vcf_reader.infos['MATECHROM'] = vcf.parser._Info('MATECHROM', '1', 'String', 'Mate chromosome')
    vcf_reader.infos['MATEID'] = vcf.parser._Info('MATEID', '1', 'String', 'Mate ID')
    vcf_writer = vcf.Writer(gzip.open(outVCF, 'wb'), vcf_reader, lineterminator='\n') if outVCF.endswith('.gz') else vcf.Writer(open(outVCF, 'w'), vcf_reader, lineterminator='\n')
    for record in vcf_reader:
        vcfInfo = record.INFO
        chrom1 = record.CHROM
        pos1 = record.POS
        svorient = record.INFO['CT']
        chrom2 = record.INFO['CHR2']
        pos2 = record.INFO['END']
        record.QUAL = record.INFO['MAPQ']
        # Fix 5' breakpoint miss-offset (DUP: POS2+1; DEL: POS1+1)
        pos1, pos2, svlen = fix_offset(record, pos1, pos2)
        # if SV length is corrected
        if svlen:
            record.INFO['SVLEN'] = svlen
        # move DEL/DUP/INV/TRA to 'SVCLASS'
        record.INFO['SVCLASS'] = record.INFO['SVTYPE']
        # annotate as breakend
        record.INFO['SVTYPE'] = 'BND'
        refFirst = 'N'
        try:
            refFirst = retrieve_nt(chrom1,pos1,pos1)
        except Exception, e:
                pass
        refMate = 'N'
        try:
            refMate = retrieve_nt(chrom2,pos2,pos2)
        except Exception, e:
                pass
        coordMate = chrom2 + ':' + str(pos2)
        coordFirst  = chrom1 + ':' + str(pos1)
        if svorient.endswith('5to5'):
            strand1 = '-'
            strand2 = '-'
            intervalMate = '['
            intervalFirst = '['
            altFirst = intervalMate + coordMate + intervalMate + refFirst
            altMate = intervalFirst + coordFirst + intervalFirst + refMate
        elif svorient.endswith('3to3'):
            strand1 = '+'
            strand2 = '+'
            intervalMate = ']' 
            intervalFirst = ']'
            altFirst = refFirst + intervalMate + coordMate + intervalMate
            altMate = refMate + intervalFirst + coordFirst + intervalFirst
        elif svorient.startswith('5to3'):
            strand1 = '-'
            strand2 = '+'
            intervalMate = ']' 
            intervalFirst = '['  
            altFirst =  intervalMate + coordMate + intervalMate + refFirst
            altMate = refMate + intervalFirst + coordFirst + intervalFirst 
        elif svorient.startswith('3to5'):
            strand1 = '+'
            strand2 = '-'
            intervalFirst = ']'
            intervalMate = '['
            altFirst =  refFirst + intervalMate + coordMate + intervalMate
            altMate =   intervalFirst + coordFirst + intervalFirst + refMate
        record.INFO.pop('END', None)
        record.INFO.pop('CHR2', None)
        recordFirst = copy.copy(record)
        recordMate = copy.copy(record)
        recordFirst.REF = refFirst
        recordFirst.ALT = [altFirst]
        recordMate.REF = refMate
        recordMate.ALT = [altMate]
        idFirst = recordFirst.ID + '_1'
        idMate = recordMate.ID + '_2'
        recordFirst.ID = idFirst
        recordMate.ID = idMate
        # deep copy INFO dict
        recordFirst.INFO = copy.copy(record.INFO)
        recordMate.INFO = copy.copy(record.INFO)
        recordFirst.INFO['MATEID'] = idMate
        recordFirst.INFO['STRAND'] = strand1 
        recordFirst.INFO['MATESTRAND'] = strand2 
        recordFirst.INFO['MATECHROM'] = chrom2
        recordFirst.INFO['MATEPOS'] = pos2
        recordMate.INFO['MATEID'] = idFirst
        recordMate.INFO['STRAND'] = strand2
        recordMate.INFO['MATESTRAND'] = strand1
        recordMate.INFO['MATECHROM'] = chrom1
        recordMate.INFO['MATEPOS'] = pos1
        svorient_pos = svorient.split('to')
        recordMate.INFO['CT'] = svorient_pos[1] + 'to' + svorient_pos[0]
        recordFirst.POS = pos1
        recordMate.POS = pos2
        recordMate.CHROM = chrom2
        vcf_writer.write_record(recordFirst)
        vcf_writer.write_record(recordMate)
    
    vcf_writer.close()

print (outVCF)

