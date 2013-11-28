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

Alternatively, a statically linked binary for Linux 64-bit is available here [delly_v0.1.2.tar.gz](http://www.embl.de/~rausch/delly_v0.1.2.tar.gz)


Running DELLY
-------------

DELLY just needs one bam file for every sample and the reference genome to identify split-reads. The output is in [vcf](http://vcftools.sourceforge.net/) format.
The SV type can be DEL, DUP, INV or JMP for deletions, tandem duplications, inversions and translocations, respectively.

`./src/delly -t DEL -o del.vcf -g <ref.fa> <sample1.sort.bam> ... <sampleN.sort.bam>`

Each bam file is assumed to be one sample. If you do have multiple bam files for a single sample please merge these bams using tools such as [Picard](http://picard.sourceforge.net/) and tag each library with a ReadGroup. To save runtime it is advisable to exclude telomere and centromere regions. For human, DELLY ships with such an exclude list.

`./src/delly -t DEL -x human.hg19.excl.tsv -o del.vcf -g <ref.fa> <sample1.sort.bam> ... <sampleN.sort.bam>`

If you omit the reference sequence DELLY skips the split-read analysis. The vcf file follows the [vcf specification](http://vcftools.sourceforge.net/specs.html) and all output fields are explained in the vcf header.

`grep "^#" del.vcf`

Delly ships with a small python script to annotate differential read-depth between SV carrier and non-carrier samples. This script is primarily meant as an example of how you can filter and annotate the final DELLY vcf file further. 

`python diffRC.py -t DEL -v del.vcf -o del.rd.vcf`

There are also external packages that consume VCF files with per-sample genotype likelihoods. The [arfer](https://github.com/ekg/arfer) package annotates, for instance, Hardy-Weinberg Equilibrium and the inbreeding coefficient, which is useful for selecting polymorphic sites in the genome across a population.

`cat del.rd.vcf | ./arfer/arfer > del.rd.arfer.vcf`


FAQ
---
* What about inversions and translocations?  
Not yet supported in DELLY v0.1.1. Please use the old single sample DELLY version available [here](http://www.embl.de/~rausch/delly.html).

* What is the smallest SV size DELLY can call?  
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, DELLY starts to call reliable SVs >=300bp.

* Can DELLY be used on a non-diploid genome?  
Yes and no. The SV site discovery works for any ploidy. However, the genotyping follows the classical hom. reference, het. and hom. alternative scheme.

* Is there any visualization of the called SVs?  
No, DELLY does not produce any graphical output. However, there are many viewers such as [IGV](http://www.broadinstitute.org/igv/) that do visualize discordantly mapped paired-ends.

* How can DELLY be used to call somatic SVs?  
Run DELLY jointly on the cancer data and the matched control sequencing data. Ideally, you include many control samples in a single run because assuming that any reference mapping artifact is recurrent, multiple control samples from different patients will help you to catch these reference-biases more easily. In the end, one just filters the tumor SVs against all SVs present in any of the control genomes. For copy-number variable events (CNVs), such as deletions and tandem duplications the additional normalized read-count genotype field (RC) can help to differentiate complex rearrangements that do not necessarily show a read-depth change from simple CNVs that have a supporting read-depth signal.

* Can DELLY be used on bwa mem alignments?  
DELLY can be used with bwa mem but you have to mark shorter split-read alignments as secondary alignments using the '-M' option, e.g. bwa mem -M ref.fa file.fq.gz.

* Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?  
DELLY expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them (e.g. the longest split alignment) is a primary alignment whereas all others are marked as secondary (flag 0x0100).


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.  
[DELLY: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)  
Bioinformatics 2012 28: i333-i339.
