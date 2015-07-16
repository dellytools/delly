from __future__ import print_function
import argparse, gzip, pysam, re, sys



def extract_reference_slices(faf_ref, f_fasta, f_out, f_bed=None):
    for coord in get_coords(f_fasta, f_bed):
        key="{}:{}-{}".format(*coord)
        seq=faf_ref.fetch(region=key)
        print(">{}".format(key), file=f_out)
        while seq:
            print(seq[:60], file=f_out)
            seq = seq[60:]


def get_coords(f_fasta, f_bed=None):
    coords = []
    numSequences =0
    # read BED file if available
    if f_bed:
        for line in f_bed:
            if line.strip():
                chrom, start, end = line.split()[:3]
                coords.append((chrom, int(start), int(end)))
        numSequences = len([l for l in f_fasta if l.startswith('>')])
    # if no BED available
    else: 
        for name, seq, qual in readfq(f_fasta):
            # split on potential delimiters
            numSequences += 1
            for head in re.split(r'\s+|>|[|]|@', name):
                match = re.match(r'^([a-zA-Z0-9_]+):(\d+)-(\d+)$', head)
                if match:
                    coords.append((match.groups()[0],
                                   int(match.groups()[1]),
                                   int(match.groups()[2])))
                    break
    if len(coords) != numSequences:
        raise Exception("Number of coordinates provided does not match number "
                        "of reads. Please specify coordinates in BED format "
                        "or behind fasta headers.")
    return coords


# External code
# taken form lh3
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None       # EDIT: read whole name line
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Prepare reference slices for '
                        'assemblies to be visualized in maze. If you provide '
                        'a BED file, coordinates are taken from it in the '
                        'given order. Without BED, I look for a pattern of '
                        'chr:start-end in the fasta headers, seperated by '
                        'white spaces.')
    parser.add_argument("-c", "--coords", metavar='BED', default=None,
                        type=argparse.FileType('r'), required=False,
                        help="Bed file with coordinates.")
    parser.add_argument("-f", "--fasta", metavar='FA',
                        type=argparse.FileType('rb'), required=True,
                        help="Bed file with coordinates.")
    parser.add_argument("-r", "--ref", metavar='FA', 
                        type=pysam.FastaFile, required=True,
                        help='Reference genome; index required.')
    parser.add_argument("-o", "--out", metavar='FA', default='-',  
                        type=argparse.FileType('w'), required=False,
                        help='Output file.')
    args = parser.parse_args()

    # .gz support for fasta
    if args.fasta != sys.stdin and args.fasta.name.endswith(".gz"):
        args.fasta = gzip.open(args.fasta.name, 'r')

    extract_reference_slices(args.ref, args.fasta, args.out, args.coords)
    
