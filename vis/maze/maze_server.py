#! /usr/bin/env python

# 2014-2015 Markus Hsi-Yang Fritz

from __future__ import division
import click
from flask import Flask, render_template
import json

app = Flask(__name__)
cfg = {}

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/data')
def data():
    import time
    time.sleep(3)
    return json.dumps({'foo': 'bar'})

@click.command()
@click.option('-p', '--port', default=5000, help='port number')
@click.option('--debug/--no-debug', default=False,
              help='run server in debug mode')
@click.option('-r', '--reference', required=True,
              help='reference (single) FASTA file')
@click.option('-q', '--query', required=True,
              help='query (multi) FASTA file')
@click.option('-c', '--coords', help='reference coordinates BED file')
@click.option('-l', '--length', default=12,
              help='match length (default: 12)')
@click.option('-m', '--matches', default='mem',
              type=click.Choice(['mem', 'mum']),
              help='match type ["mem", "mum"], default=mem')
def cli(port, debug, reference, query, coords, length, matches):
    global cfg
    cfg = {
        'reference': reference,
        'query': query,
        'length': length,
        'matches': matches
    }
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
