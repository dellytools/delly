#! /usr/bin/env python

# Copyright (c) 2014 Markus Hsi-Yang Fritz

import click
from flask import Flask, render_template

app = Flask(__name__)
h5_fns = []

@app.route('/')
def index():
    return render_template('index.html')

@click.command()
@click.option('-p', '--port', default=5000, help='port number')
@click.option('-d', '--debug', default=False, help='run server in debug mode')
@click.option('h5s', '-f', '--file', multiple=True, required=True,
              help='HDF5 coverage file')
def cli(port, debug, h5s):
    global h5_fns
    h5_fns = list(h5s)
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
