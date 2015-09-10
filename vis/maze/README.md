# maze: match visualizer
In a similar manner to
[dot plots](http://en.wikipedia.org/wiki/Dot_plot_%28bioinformatics%29),
`maze` highlights local sequence similarity between two DNA sequences.
In particular, maximal exact substring matches are computed with
[MUMmer3](http://mummer.sourceforge.net/) and visualised.

In the context of Delly this is most useful for exploring breakpoints
or entire SVs constructed with the local assembly pipeline,
but `maze` has also been instrumental in analysing SVs from
"long read" data (PacBio, MinION) or Sanger reads.

Our breakpoint module, which is still in an experimental state, highlights
the breakpoints of your assemblies relative to another sequence (e.g. the 
reference) using local alignments computed with [last-split](http://last.cbrc.jp/).

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


For the breakpoint module to work, you additionally need the `last` and 
`last-split` binaries in a version >= 584.

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

## Coming from Delly
Delly currently ships with a beta version of an 
[assembly pipeline](https://github.com/tobiasrausch/delly/tree/master/assembly). 
This pipeline assembles reads around SV calls using `spades` and outputs
contigs for each locus into FASTA files. See the respective readme for
more information.

If you have `samtools` installed, the same pipeline will also provide a
second FASTA per locus that contains the reference allele. Simply drag 
these two files into maze and you're there.

## Inspecting several different loci
As said before, `maze` works with two FASTA files, the one containing your 
assemblies/long reads and the other one containing the piece of reference 
to compare them to. However, if you supply two FASTA files with **equally 
many entries** you will trigger a pairwise comparison between the sequences. 	

There is a helper script in the `maze` folder that can generate these
reference slices for you. All you need to provide is the FASTA file
with your assemblies/long reads and the coordinates where they stem from. 
These coordinates can for instance be written directly in the FASTA files 
behind each header separated by a space, e.g.

    >assembly1 chr2:1234567-1239999
    ACTCATGCAT...

Then run

    $ ./extract_reference_slices.py -r hg19.fa -f assemblies.fa > assemblies.reference.fa

Replace `hg19.fa` with the refernce genome you are using and `assemblies.fa`
with your reads/assemblies. 

In case the coordinates are not included in the FASTA header you can 
also input them as a table in BED format. Just make sure that the table 
has the same order as the reads/assemblies in your FASTA file. The 
command is 

    $ ./extract_reference_slices.py -r hg19.fa -f assemblies.fa -c coordinates.bed > assemblies.reference.fa 

Finally, drag the two files `assemblies.fa` and `assemblies.reference.fa`
into the setup menu of `maze`.
