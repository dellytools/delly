#! /usr/bin/env python

# Copyright (c) 2014 Markus Hsi-Yang Fritz

from __future__ import division
import click
from flask import Flask, render_template, request
import h5py
import os
import json
from collections import defaultdict
import numpy as np
import math
import re
import pysam
import csv

app = Flask(__name__)
cfg = dict()


def get_smpl_names(h5_fns, use_fns):
    if use_fns:
        return [os.path.basename(fn) for fn in h5_fns]
    ss = []
    for h5_fn in h5_fns:
        with h5py.File(h5_fn) as f:
            ss.append(f.attrs['sample'])
    return ss


def get_common_chroms(s1, s2):
    with h5py.File(cfg['smpl_to_file'][s1]) as f:
        chroms1 = set(c for c in f)
    with h5py.File(cfg['smpl_to_file'][s2]) as f:
        chroms2 = set(c for c in f)
    return sorted(list(chroms1 & chroms2))


@app.route('/chroms/<s1>/<s2>')
def chroms(s1, s2):
    return json.dumps(get_common_chroms(s1, s2))

chr2_re = re.compile(r'CHR2=([^;\s]+)')
end_re = re.compile(r'\bEND=(\d+)')
ct_re = re.compile(r'\bCT=([^;\s]+)')

@app.route('/calls/<chrom>')
def calls(chrom):
    dataset = []
    if cfg['vcf']:
        f = pysam.Tabixfile(cfg['vcf'])
        for row in csv.reader(f.fetch(str(chrom)), delimiter='\t'):
            chrom2 = chr2_re.search(row[7]).group(1)
            if chrom2 != chrom:
                continue
            start_call = int(row[1])
            end_call = int(end_re.search(row[7]).group(1))
            m = ct_re.search(row[7])
            if m:
                ct = m.group(1)
            else:
                ct = 'none'
            sv_type = row[4].lstrip('<').rstrip('>')
            dataset.append({
                'start': start_call,
                'end': end_call,
                'type': sv_type,
                'ct': ct
            })

    return json.dumps(dataset)

@app.route('/depth/<s1>/<s2>/<c>')
def depth(s1, s2, c):
    d = {
        'chrom': c,
        'chrom_len': -1,
    }
    with h5py.File(cfg['smpl_to_file'][s1], 'r') as f1, \
            h5py.File(cfg['smpl_to_file'][s2], 'r') as f2:
        assert c in f1 and c in f2
        d['chrom_len'] = int(f1[c].attrs['length'])
        # TODO: this needs some more thought...
        s = request.args.get('start', 1, type=int)
        e = request.args.get('end', d['chrom_len'], type=int)
        n_req = request.args.get('n', 10000, type=int)

        bin_start = (s-1) // 100
        bin_end = (e-1) // 100
        n_bins = bin_end - bin_start + 1

        print s, e, n_req, bin_start, bin_end

        # TODO:
        # * benchmark compressed files
        # * what about holding everything in memory?
        x = f1['/{}/read_counts'.format(c)][bin_start:bin_end+1]
        y = f2['/{}/read_counts'.format(c)][bin_start:bin_end+1]

        if n_bins <= n_req:
            x_sum = x
            y_sum = y
        else:
            chunk_size = int(round(n_bins/n_req))
            n_chunks = int(math.ceil(n_bins / chunk_size))
            pad = -n_bins % chunk_size
            print 'x_len', n_bins, len(x), 'chunk_size', chunk_size, 'n_chunks', n_chunks, 'pad', pad
            x = np.pad(x, (0, pad), mode='constant')
            x_sum = np.sum(np.split(x, n_chunks), axis=1)
            y = np.pad(y, (0, pad), mode='constant')
            y_sum = np.sum(np.split(y, n_chunks), axis=1)

        d['ratios'] = [r if np.isfinite(r) else None
                       for r in np.log2(x_sum/y_sum).tolist()]

    return json.dumps(d)


@app.route('/')
def index():
    smpl_names = get_smpl_names(cfg['h5_fns'], cfg['use_fns'])
    cfg['smpl_to_file'] = dict()
    for smpl, fn in zip(smpl_names, cfg['h5_fns']):
        cfg['smpl_to_file'][smpl] = fn
    chroms = get_common_chroms(smpl_names[0], smpl_names[0])
    return render_template('index.html',
                           samples=smpl_names,
                           chromosomes=chroms)


@click.command()
@click.option('-p', '--port', default=5000, help='port number')
@click.option('--debug/--no-debug', default=False,
              help='run server in debug mode')
@click.option('h5s', '-f', '--file', multiple=True, required=True,
              type=click.Path(exists=True), help='HDF5 coverage file')
@click.option('--usefilenames/--no-usefilenames', default=False,
              help='use file names as sample names')
@click.option('-v', '--vcf', help='VCF file')
def cli(port, debug, h5s, usefilenames, vcf):
    global cfg
    cfg['h5_fns'] = list(h5s)
    cfg['use_fns'] = usefilenames
    cfg['vcf'] = vcf
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
