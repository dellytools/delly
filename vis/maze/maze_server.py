#! /usr/bin/env python

# 2014-2015 Markus Hsi-Yang Fritz

from __future__ import division
import click
from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@click.command()
@click.option('-p', '--port', default=5000, help='port number')
@click.option('--debug/--no-debug', default=False,
              help='run server in debug mode')
def cli(port, debug):
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
