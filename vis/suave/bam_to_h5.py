#! /usr/bin/env python

import pysam
import numpy as np
import h5py
import click
import re
import os.path

NUM_BINS_LIST = [10000, 100000]
median_vec = np.vectorize(np.median)


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


def set_dset_attrs(dset, arr):
    dset.attrs['min'] = np.min(arr)
    dset.attrs['max'] = np.max(arr)
    dset.attrs['mean'] = np.mean(arr)
    dset.attrs['median'] = np.median(arr)


def bam_to_h5(bamfile,
              sample=None,
              chroms_keep=None,
              compression=None,
              outfile=None):
    out_fn = gen_outfile_name(outfile, bamfile)
    sample_name = gen_sample_name(sample, bamfile)

    f_h5 = h5py.File(out_fn)
    f_h5.attrs['sample'] = sample_name

    f_bam = pysam.Samfile(bamfile)

    for sq in f_bam.header['SQ']:
        chrom = sq['SN']
        if chroms_keep is not None and chrom not in chroms_keep:
            continue
        click.echo('processing ' + chrom)

        chrom_len = sq['LN']
        f_h5.create_group(chrom)
        f_h5[chrom].attrs['length'] = np.uint32(chrom_len)
        cnts = np.zeros(chrom_len, dtype='uint16')

        for pile in f_bam.pileup(chrom):
            cnts[pile.pos] = pile.n

        kwargs = {}
        if compression:
            kwargs['compression'] = compression

        for nbins in NUM_BINS_LIST:
            nbins_final = min(chrom_len, nbins)
            medians = np.array(median_vec(np.array_split(cnts, nbins_final)),
                               dtype=np.uint16)
            kwargs['data'] = medians
            f_h5.create_dataset('{}/{}'.format(chrom, nbins), **kwargs)
            set_dset_attrs(f_h5['{}/{}'.format(chrom, nbins)], medians)

    f_bam.close()


@click.command()
@click.option('-s', '--sample', help='sample name')
@click.option('-c', '--chrom', multiple=True, help='chromosome to process')
@click.option('-d', '--compression', type=click.Choice(['gzip', 'lzf']),
              help='HDF5 file compression')
@click.option('-o', '--outfile', help='output HDF5 file name')
@click.argument('bamfile')
def cli(bamfile, sample, chrom, compression, outfile):
    """Convert BAM to coverage HDF5."""
    chroms_keep = set(chrom) if chrom else None
    bam_to_h5(bamfile, sample, chroms_keep, compression, outfile)

if __name__ == '__main__':
    cli()
