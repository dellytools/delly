#! /usr/bin/env python

import pysam
import numpy as np
import h5py
import click
import re
import os.path

def gen_outfile_name(out, bam):
    if out:
        return out
    bam = os.path.basename(bam)
    if bam.endswith('.bam'):
        return re.sub(r'bam$', 'h5', bam)
    return bam + '.h5'

def gen_sample_name(smpl, bam):
    if smpl:
        return smpl
    bam = os.path.basename(bam)
    return re.sub(r'\.bam$', '', bam)

def bam_to_h5(bamfile,
              sample=None,
              chrom=None,
              compression=None,
              outfile=None):
    out_fn = gen_outfile_name(outfile, bamfile)
    sample_name = gen_sample_name(sample, bamfile)
    f_h5 = h5py.File(out_fn)
    f_h5.attrs['sample'] = sample_name

@click.command()
@click.option('-s', '--sample', help='sample name')
@click.option('-c', '--chrom', multiple=True, help='chromosome to process')
@click.option('-d', '--compression', type=click.Choice(['gzip', 'lzf']),
              help='HDF5 file compression')
@click.option('-o', '--outfile', help='output HDF5 file name')
@click.argument('bamfile')
def cli(bamfile, sample, chrom, compression, outfile):
    """Convert BAM to coverage HDF5."""
    bam_to_h5(bamfile, sample, chrom, compression, outfile)

if __name__ == '__main__':
    cli()
