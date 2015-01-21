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


@app.route('/data/<s1>/<s2>/<c>')
def data(s1, s2, c):
    d = {
        'chrom': c,
        'chrom_len': -1,
    }
    with h5py.File(cfg['smpl_to_file'][s1], 'r') as f1, \
            h5py.File(cfg['smpl_to_file'][s2], 'r') as f2:
        assert c in f1 and c in f2
        d['chrom_len'] = int(f1[c].attrs['length'])
        # TODO this all needs some more thought...
        s = request.args.get('start', 1, type=int)
        e = request.args.get('end', d['chrom_len'], type=int)
        l = e - s + 1
        n_100bp = int(math.ceil(l / 100))
        n = min(request.args.get('n', 25000, type=int), n_100bp)

        binStart = (s-1) // 100
        binEnd = (e-1) // 100

        print(s, e, n, binStart, binEnd)

        x = f1['/{}/read_counts'.format(c)][binStart:binEnd+1]
        print(x)
        pad = -len(x) % n
        x = np.pad(x, (0, pad), mode='constant')
        x_sum = np.sum(np.split(x, n), axis=1)

        y = f2['/{}/read_counts'.format(c)][binStart:binEnd+1]
        y = np.pad(y, (0, pad), mode='constant')
        y_sum = np.sum(np.split(y, n), axis=1)

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
def cli(port, debug, h5s, usefilenames):
    global cfg
    cfg['h5_fns'] = list(h5s)
    cfg['use_fns'] = usefilenames
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
