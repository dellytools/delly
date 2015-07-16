# maze: match visualizer
In a similar manner to
[dot plots](http://en.wikipedia.org/wiki/Dot_plot_%28bioinformatics%29),
`maze` highlights local sequence similarity between two DNA sequences.
In particular, maximal exact substring matches are computed with
[MUMmer3](http://mummer.sourceforge.net/) and visualised.

In the context of Delly this is most useful for exploring breakpoints
or entire SVs constructed with the local assembly pipeline,
but `maze` has also been instrumental in analysing SVs from
"long read" data (PacBio, MinION).

## Getting started
For the server component, you will need Python2.7+ and the
`Flask` and `click` modules. An easy way to set this up is with a
[Virtual Environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/):

    $ cd maze
    $ virtualenv venv
    $ source venv/bin/activate
    $ pip install -r requirements/python2.7.txt

You will also need a "modern" web browser. Chrome, Firefox and Safari 
should all work. Technically, Opera and IE9+ should be fine as well 
although these haven't been tested.

`maze` uses `MUMmer3` to compute matches, thus you need the
`mummer` binary in your `$PATH`. For breakpoint analysis, you
will also need the `lastdb`, `lastal` and `last-split` binaries 
from the [LAST](http://last.cbrc.jp/) package.


## Running the app
First off, start the server:

    $ python maze_server.py

By default, port 5000 is used, but you can specify this with the
`--port` option. Now, open your web browser and head to http://localhost:5000/.

Click on the setup icon and upload the reference and query FASTA
files via drag-and-drop. While the query file can contain 
multiple sequences, the reference file should only contain one sequence
(well, it can contain more, but only the first sequence is used...).
Both, raw and gzip files are supported.

Change the match type and length if you want to and you're ready to go.
