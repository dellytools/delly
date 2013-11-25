DELLY
=====

DELLY is an integrated structural variant prediction method that can detect deletions, tandem duplications, inversions and translocations
at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately
delineate genomic rearrangements throughout the genome.


DELLY dependencies
------------------

* BamTools, (https://github.com/pezmaster31/bamtools)
* Boost C++ Libraries, (www.boost.org)
* zlib compression library, (www.zlib.net)
* kseq library to parse FASTA/FASTQ, (http://lh3lh3.users.sourceforge.net/parsefastq.shtml)

Installing DELLY
----------------

The easiest way to get DELLY is to clone the git repository

`git clone https://github.com/tobiasrausch/delly.git`

Then you need to install Boost, BamTools and KSEQ and update the paths in the Makefile.
Building DELLY just requires

`make -B src/delly`

Alternatively, statically linked binaries for Linux 64-bit are available here [delly_v0.0.11.tar.gz](http://www.embl.de/~rausch/delly_v0.0.11.tar.gz)


Running DELLY
-------------

`./src/delly -t DEL -o del.vcf -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> <sample3.sort.bam> ...`


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.

[DELLY: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)

Bioinformatics 2012 28: i333-i339.

