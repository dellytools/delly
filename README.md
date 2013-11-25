DELLY
=====

DELLY: Structural variant discovery by integrated paired-end and split-read analysis

DELLY dependencies
==================

* BamTools, https://github.com/pezmaster31/bamtools

Please follow the BamTools installation instructions.
For static binaries make sure you use "make BamTools-static".

* Boost C++ Libraries, www.boost.org

Please follow the boost installation instructions.
DELLY does require the boost iostreams library.

* zlib compression library, www.zlib.net

Please follow the zlib installation instructions: zlib.net

* kseq library to parse FASTA/FASTQ, http://lh3lh3.users.sourceforge.net/parsefastq.shtml

Installing DELLY
================

Please update the paths to Boost, BamTools and KSEQ in the Makefile.

make -B src/delly

Running DELLY
=============

