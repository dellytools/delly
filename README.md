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

DELLY just needs one bam file for every sample and the reference genome to identify split-reads. The output is in [vcf](http://vcftools.sourceforge.net/) format.
The SV type can be DEL, DUP, INV or JMP for deletions, tandem duplications, inversions and translocations, respectively.

`./src/delly -t DEL -o del.vcf -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> <sample3.sort.bam> ...`

Each bam file is assumed to be one sample. If you do have multiple bam files for a single sample please merge these bams using tools such as [Picard](http://picard.sourceforge.net/) and tag each library with a ReadGroup. To save runtime it is advisable to exclude telomere and centromere regions. For human, DELLY ships with such an exclude list.

`./src/delly -t DEL -x human.hg19.excl.tsv -o del.vcf -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> <sample3.sort.bam> ...`

If you omit the reference sequence DELLY skips the split-read analysis. The vcf output fields are explained in the vcf header and the vcf file itself follows the [vcf specification](http://vcftools.sourceforge.net/specs.html).

`grep "^#" del.vcf`


FAQ
---
* What is the smallest SV size DELLY can call?</b><br>
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, DELLY starts to call reliable SVs >=300bp.

* Can DELLY be used on a non-diploid genome?

Yes and no. The site discovery works for any ploidy. However, the genotyping follow the hom. reference, het. and hom. alternative scheme.

* Is there any visualization of the called SVs?
No, DELLY does not produce any graphical output. However, there are many viewers such as the [IGV](http://www.broadinstitute.org/igv/) that do visualize discordantly mapped paired-ends.

* How can DELLY be used to call somatic SVs?

Run DELLY on as many cancer and control genomes you have access to and then filter the tumor SVs of a single sample against all SVs called in all control genomes. In practice, this appears to be one of the best options to derive a high quality set of somatic SVs. For copy-number variable events (CNVs) such as deletions and tandem duplications a further annotation with read-depth may be helpful, although most complex rearrangements do not necessarily show such a read-depth change.

* How can DELLY be used to call germline SVs?

ToDo

* Can DELLY be used on bwa mem alignments?</b><br>
DELLY can be used with bwa mem but you have to mark shorter split-read alignments as secondary alignments using the '-M' option, e.g. bwa mem -M ref.fa file.fq.gz.

* Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?

DELLY expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them (e.g. the longest split alignment) is a primary alignment whereas all others are marked as secondary (flag 0x0100).


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.

[DELLY: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)

Bioinformatics 2012 28: i333-i339.
