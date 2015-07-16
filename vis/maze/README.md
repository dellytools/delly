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
`mummer` binary in your `$PATH`.

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

## Coming from Delly2
Delly2 currently ships with a beta version of an assembly pipeline (see
respective README).
This pipeline assembles reads around SV calls using `spades` and outputs
a list of contigs for each locus into a FASTA file. If you have `samtools`
installed, the same pipeline will also output the corresponding refernece
slices.  Simply drag these two files into maze to inspect your assemblies.

## Inspecting several different loci
`maze` works with two FASTA files, one containing your assemblies or long 
reads and the other one containing the piece of reference to compare them 
to. However, if you specify two FASTA files with equally many entries you
will trigger a pairwise comparison. Corresponding entries must be in the 
same order. 

There is a helper script in the `maze` folder that can generate these
reference slices for you. All you need to provide is the fasta file
with your assemblies/long reads and the coordinates. These coordinates can
for instance be written directly in the fasta files behind each fasta
header in the format `chrom:start-end` (separated from the actual name by
spaces). Then run

    $ ./extract_reference_slices.py -r hg19.fa -f assemblies.fa > assemblies.reference.fa

Replace `hg19.fa` with the refernce genome you are using and `assemblies.fa`
with your reads/assemblies. In case the coordinates are not included in the
fasta header you can also specify them as a table in BED format. Just make 
sure that the table has the same order as the reads/assemblies in your fasta
file. The command is 

    $ ./extract_reference_slices.py -r hg19.fa -f assemblies.fa -c coordinates.bed > assemblies.reference.fa 

Finally, drag the two files `assemblies.fa` and `assemblies.reference.fa`
into the setup menu of `maze`.
