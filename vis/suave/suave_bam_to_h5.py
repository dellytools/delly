#! /usr/bin/env python

# Copyright (c) 2014 Markus Hsi-Yang Fritz

from __future__ import division
import pysam
import numpy as np
import h5py
import click
import re
import os.path
import math
import sys

try:
    import mpi4py
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpi_rank = comm.rank
    mpi_size = comm.size
except ImportError:
    mpi_rank = 0
    mpi_size = 1

CHUNK_SIZE = 100


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
              compression=None,
              outfile=None,
              force=False):
    out_fn = gen_outfile_name(outfile, bamfile)
    sample_name = gen_sample_name(sample, bamfile)

    if os.path.exists(out_fn) and not force:
        click.echo('error: file "{}" exists; use --force to overwrite'\
                   .format(out_fn), err=True)
        sys.exit(1)

    if mpi_size == 1:
        f_h5 = h5py.File(out_fn, 'w')
    else:
        f_h5 = h5py.File(out_fn, 'w', driver='mpio', comm=comm)

    if mpi_rank == 0:
        f_h5.attrs['sample'] = sample_name
    f_bam = pysam.Samfile(bamfile)

    kwargs = {}
    if compression:
        kwargs['compression'] = compression

    chroms = [(sq['SN'], sq['LN']) for sq in f_bam.header['SQ']]
    chroms_chunks = np.array_split(chroms, mpi_size)
    mpi_rank_chroms_chunk = chroms_chunks[mpi_rank]

    for chrom, chrom_len in mpi_rank_chroms_chunk:
        chrom = str(chrom)

        chrom_len = int(chrom_len)
        click.echo('[t{}] processing {}'.format(mpi_rank, chrom), err=True)

        f_h5.create_group(chrom)
        f_h5[chrom].attrs['length'] = np.uint32(chrom_len)


        nbins = int(math.ceil(chrom_len/CHUNK_SIZE))
        cnts = np.zeros(nbins, dtype='uint32')

        for read in f_bam.fetch(chrom):
            if read.is_secondary or read.is_duplicate or read.is_qcfail:
                continue
            cnts[read.pos // CHUNK_SIZE] += 1

        kwargs['data'] = cnts
        f_h5.create_dataset('{}/read_counts'.format(chrom), **kwargs)
        dset = f_h5['{}/read_counts'.format(chrom)]
        dset.attrs['min'] = np.min(cnts)
        dset.attrs['max'] = np.max(cnts)
        dset.attrs['mean'] = np.mean(cnts)
        dset.attrs['median'] = np.median(cnts)
        dset.attrs['sum'] = np.sum(cnts)

    # TODO actually, don't store this as it's easy/fast to compute on the fly
    # store genome-wide read count for normalization
    read_cnt_tot = 0
    for chrom in f_h5:
        read_cnt_tot += f_h5['{}/read_counts'.format(chrom)].attrs['sum']
    f_h5.attrs['read_counts_sum'] = np.array(read_cnt_tot, dtype=np.uint64)

    f_bam.close()


@click.command()
@click.option('-s', '--sample', help='sample name')
@click.option('-c', '--compression', type=click.Choice(['gzip', 'lzf']),
              help='HDF5 file compression')
@click.option('-o', '--outfile', help='output HDF5 file name')
@click.option('-f', '--force/--no-force', default=False,
              help='overwrite existing file')
@click.argument('bamfile')
def cli(bamfile, sample, compression, outfile, force):
    """Convert BAM to coverage HDF5."""
    bam_to_h5(bamfile, sample, compression, outfile, force)

if __name__ == '__main__':
    cli()
