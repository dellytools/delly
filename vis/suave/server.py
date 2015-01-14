#! /usr/bin/env python

# Copyright (c) 2014 Markus Hsi-Yang Fritz

import click
from flask import Flask

app = Flask(__name__)

@click.command()
@click.option('-p', '--port', default=5000, help='port number')
@click.option('-d', '--debug', default=False, help='run server in debug mode')
def cli(port, debug):
    app.run(port=port, debug=debug)

if __name__ == '__main__':
    cli()
