#! /usr/bin/env python

from __future__ import print_function
import argparse
import vcf
import pysam
import os
import collections
import errno
import numpy
import random

def rev_compl(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(basecomplement[base] for base in s[::-1])


def rp(isRead1, isReverse, isMateReverse, pos, matePos):
    if isRead1:
        if not isReverse:
            if not isMateReverse:
                if pos < matePos:
                    return 0
                else:
                    return 1
            else:
                if pos < matePos:
                    return 2
                else:
                    return 3
        else:
            if not isMateReverse:
                if pos > matePos:
                    return 2
                else:
                    return 3
            else:
                if pos > matePos:
                    return 0
                else:
                    return 1
    else:
        if not isReverse:
            if not isMateReverse:
                if pos < matePos:
                    return 1
                else:
                    return 0
            else:
                if pos < matePos:
                    return 2
                else:
                    return 3
        else:
            if not isMateReverse:
                if pos > matePos:
                    return 2
                else:
                    return 3
            else:
                if pos > matePos:
                    return 1
                else:
                    return 0

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

# Parse command line
parser = argparse.ArgumentParser(description='VCF to BED converter listing all carrier samples.')
parser.add_argument('-v', '--vcf', metavar='variants.vcf', required=True, dest='vcfFile', help='input vcf file (required)')
parser.add_argument('-s', '--samples', metavar='sample.info', required=True, dest='sampleToBam', help='sample-bam list (required)')
parser.add_argument('-u', '--unmapped', metavar='unmapped.info', required=False, dest='unmappedToBam', help='unmapped sample-bam list (optional)')
parser.add_argument('-m', '--maxsize', metavar='50000', required=False, dest='maxSize', help='max. size (optional)')
parser.add_argument('-f', '--filter', dest='siteFilter', action='store_true', help='Filter sites for PASS')
args = parser.parse_args()

# Command-line args
bpPrecision = 2000
maxSize = 10000
if args.maxSize:
    maxSize = int(args.maxSize)

# Parse sample list
sampleToBam = dict()
if args.sampleToBam:
    with open(args.sampleToBam) as f:
        for line in f:
            fields = line.strip().split('\t')
            sampleToBam[fields[0]] = fields[1]

# Parse sample list
unmappedToBam = dict()
if args.unmappedToBam:
    with open(args.unmappedToBam) as f:
        for line in f:
            fields = line.strip().split('\t')
            unmappedToBam[fields[0]] = fields[1]

# Parse VCF
sv = collections.defaultdict(list)
observedSamples = set()
if args.vcfFile:
    vcf_reader = vcf.Reader(open(args.vcfFile), 'r', compressed=True) if args.vcfFile.endswith('.gz') else vcf.Reader(open(args.vcfFile), 'r', compressed=False)
    for record in vcf_reader:
        if (record.INFO['SVLEN'] <= maxSize) and ((not args.siteFilter) or (len(record.FILTER) == 0)):
            if record.INFO['SVTYPE'] == "TRA":
                continue
            carrier = set()
            dupBounds = set()
            for call in record.samples:
                if (call.sample in sampleToBam.keys()) and (call.called) and (call.gt_type != 0):
                    carrier.add(call.sample)
                    observedSamples.add(call.sample)
            if len(carrier):
                if (record.INFO['SVTYPE'] == "DUP") or (record.INFO['SVTYPE'] == "INVDUP"):
                    dupBounds.add(record.POS)
                    dupBounds.add(record.INFO['END'])
                try:
                    sv2 = record.INFO['SV2INFO']
                except KeyError:
                    sv2 = None
                if sv2 is not None:
                    svStart = min(record.POS, int(sv2[1]))
                    svEnd = max(record.INFO['END'], int(sv2[2]))
                    if (sv2[3] == "DUP") or (sv2[3] == "INVDUP"):
                        dupBounds.add(int(sv2[1]))
                        dupBounds.add(int(sv2[2]))
                else:
                    svStart = record.POS
                    svEnd = record.INFO['END']
                sv[record.CHROM].append((svStart, svEnd, record.ID, record.INFO['SVTYPE'], carrier, dupBounds))

# Parse only bams for carriers
sampleToBam = {k: sampleToBam[k] for k in set(sampleToBam.keys()).intersection(observedSamples)}

# Get clip pileups
sampleProps = dict()
if len(sv.keys()):
    for chrom in sv.keys():
        for (start, end, svID, svt, carrier, dupBounds) in sv[chrom]:
            print("Processing: " + svID + " " + chrom + ":" + str(start) + "-" + str(end))

            # Collect clip positions
            windowStart = max(0, start-bpPrecision)
            windowEnd = end+bpPrecision
            startBp = collections.Counter()
            endBp = collections.Counter()
            for sample in carrier:
                if os.path.exists(sampleToBam[sample]):
                    isize = list()
                    rpCount = collections.Counter()
                    samfile = pysam.Samfile(sampleToBam[sample])
                    for read in samfile.fetch(chrom, windowStart, windowEnd):
                        if (not read.is_unmapped) and (not read.is_qcfail) and (not read.is_duplicate) and (read.is_paired):
                            clip = read.pos
                            for (cigarType, cigarLength) in read.cigar:
                                # M: 0, I: 1, D: 2, N: 3, S: 4, H: 5, P:6
                                if cigarType in [0, 2, 3]:
                                    clip += cigarLength
                                elif cigarType in [1, 6]:
                                    pass
                                elif cigarType in [4, 5]:
                                    startBp[clip] += 1
                            if (not read.mate_is_unmapped) and (read.is_read1):
                                rpCount[rp(read.is_read1, read.is_reverse, read.mate_is_reverse, read.pos, read.mpos)] += 1
                                isize.append(abs(read.tlen))
                    sampleProps[sample] = (numpy.mean(isize)+3*numpy.std(isize), rpCount.most_common(1)[0][0])
                    samfile.close()
            clips = numpy.array(startBp.values() + endBp.values())
            clipCut = max(1, int(numpy.mean(clips)+5*numpy.std(clips)))
            clips = list()
            for sebp in startBp.most_common(10) + endBp.most_common(10):
                if sebp[1] > clipCut:
                    clips.append(sebp[0])
            if len(clips) < 2:
                continue

            # Collect reads for assembly
            fq1 = args.vcfFile.strip().split('.')[0] + "." + svID + ".1.fastq"
            fq2 = args.vcfFile.strip().split('.')[0] + "." + svID + ".2.fastq"
            silentremove(fq1)
            silentremove(fq2)
            for sample in carrier:
                if os.path.exists(sampleToBam[sample]):
                    read12 = [dict(), dict()]
                    isizeCut = sampleProps[sample][0]
                    rpDef = sampleProps[sample][1]
                    samfile = pysam.Samfile(sampleToBam[sample])
                    for read in samfile.fetch(chrom, windowStart, windowEnd):
                        if (not read.is_unmapped) and (not read.is_qcfail) and (not read.is_duplicate) and (read.is_paired):
                            cigarTypes = set([cigarType for (cigarType, cigarLength) in read.cigar])
                            if 4 in cigarTypes:
                                # Keep all clipped reads
                                read12[int(read.is_read2)][read.qname] = (read.seq, read.qual, read.is_reverse)
                            else:
                                if (not read.mate_is_unmapped) and (read.mpos > windowStart) and (read.mpos < windowEnd):
                                    normalPair = False
                                    for clipPos in clips:
                                        if (clipPos > min(read.pos, read.mpos)) and (clipPos < max(read.aend, read.mpos + read.alen)) and (abs(read.tlen) < isizeCut) and (rp(read.is_read1, read.is_reverse, read.mate_is_reverse, read.pos, read.mpos) == rpDef):
                                            normalPair = True
                                            # Normal pair that should be removed, except for DUPs
                                            for dB in dupBounds:
                                                if (dB > min(read.pos, read.mpos)) and (dB < max(read.aend, read.mpos + read.alen)) and (random.random() <= 0.25):
                                                    normalPair = False
                                            break
                                    if not normalPair:
                                        read12[int(read.is_read2)][read.qname] = (read.seq, read.qual, read.is_reverse)
                    # Fetch skipped mapped reads
                    missingReads = set(read12[0].keys()).union(set(read12[1].keys())).difference(set(read12[0].keys()).intersection(set(read12[1].keys())))
                    if len(missingReads):
                        for read in samfile.fetch(chrom, windowStart, windowEnd):
                            if read.qname in missingReads:
                                read12[int(read.is_read2)][read.qname] = (read.seq, read.qual, read.is_reverse)
                    samfile.close()
                    # Fetch skipped unmapped reads
                    missingReads = set(read12[0].keys()).union(set(read12[1].keys())).difference(set(read12[0].keys()).intersection(set(read12[1].keys())))
                    if len(missingReads):
                        if (sample in unmappedToBam.keys()) and (os.path.exists(unmappedToBam[sample])):
                            samfile = pysam.Samfile(unmappedToBam[sample])
                            for read in samfile.fetch(until_eof=True):
                                if (read.is_unmapped) and (read.qname in missingReads):
                                    read12[int(read.is_read2)][read.qname] = (read.seq, read.qual, read.is_reverse)
                    # Keep only the paired reads
                    readStore = set(read12[0].keys()).intersection(set(read12[1].keys()))
                    with open(fq1, "a") as fastq1:
                        for qname in readStore:
                            (seq, qual, rev) = read12[0][qname]
                            if rev:
                                fastq1.write("@%s\n%s\n+\n%s\n" % (qname, rev_compl(seq), qual[::-1]))
                            else:
                                fastq1.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
                    with open(fq2, "a") as fastq2:
                        for qname in readStore:
                            (seq, qual, rev) = read12[1][qname]
                            if rev:
                                fastq2.write("@%s\n%s\n+\n%s\n" % (qname, rev_compl(seq), qual[::-1]))
                            else:
                                fastq2.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
