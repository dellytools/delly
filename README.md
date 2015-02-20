Delly
=====

Delly is an integrated structural variant prediction method that can detect deletions, tandem duplications, inversions and translocations
at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately
delineate genomic rearrangements throughout the genome.


Installing Delly
----------------

The easiest way to get Delly is to download a statically linked binary from the [Delly github release page](https://github.com/tobiasrausch/delly/releases/).
Alternatively, you can build Delly from source. Delly dependencies are included as submodules so you need to do a recursive clone. 

`git clone --recursive https://github.com/tobiasrausch/delly.git`

`cd delly/`

`make all`

There is a Delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users) for questions and a few wiki pages on some additional [tools](https://github.com/tobiasrausch/delly/wiki). 


Delly multi-threading mode
--------------------------
Delly supports parallel computing using the OpenMP API (www.openmp.org).

`make PARALLEL=1 -B src/delly`

There is also a statically linked, multi-threaded binary for Linux 64-bit available under [releases](https://github.com/tobiasrausch/delly/releases/).


You can set the number of threads using the environment variable OMP_NUM_THREADS in your shell.

`export OMP_NUM_THREADS=3`

Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. 

Running Delly
-------------

Delly just needs one bam file for every sample and the reference genome to identify split-reads. The output is in [vcf](http://vcftools.sourceforge.net/) format.
The SV type can be DEL, DUP, INV or TRA for deletions, tandem duplications, inversions and translocations, respectively.

`./src/delly -t DEL -o del.vcf -g <ref.fa> <sample1.sort.bam> ... <sampleN.sort.bam>`

Each bam file is assumed to be one sample. If you do have multiple bam files for a single sample please merge these bams using tools such as [Picard](http://picard.sourceforge.net/) and tag each library with a ReadGroup. To save runtime it is advisable to exclude telomere and centromere regions. For human, Delly ships with such an exclude list.

`./src/delly -t DEL -x human.hg19.excl.tsv -o del.vcf -g <ref.fa> <sample1.bam> ... <sampleN.bam>`

If you omit the reference sequence Delly skips the split-read analysis. The vcf file follows the [vcf specification](http://vcftools.sourceforge.net/specs.html) and all output fields are explained in the vcf header.

`grep "^#" del.vcf`

Delly ships with two small python scripts, one to filter somatic variants for tumor/normal comparisons and one to filter confident SV sites in population sequencing.

`python python/somaticFilter.py -v del.vcf -o del.somatic.vcf -t DEL -f`

`python python/populationFilter.py -v del.vcf -o del.sites.vcf`

These python scripts are primarily meant as an example of how you can filter and annotate the final Delly vcf file further. They may require some fine-tuning depending on your application.


FAQ
---
* What is the smallest SV size Delly can call?  
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, Delly starts to call reliable SVs >=300bp.

* Can Delly be used on a non-diploid genome?  
Yes and no. The SV site discovery works for any ploidy. However, Delly's genotyping model assumes hom. reference, het. and hom. alternative.

* How can Delly be used to call somatic SVs?  
Run Delly jointly on the cancer data and the matched control sequencing data. Ideally, you include many control samples in a single run because assuming that any reference mapping artifact is recurrent, multiple control samples from different patients will help you to catch these reference-biases more easily. In the end, one just filters the tumor SVs against all SVs present in any of the control genomes as exemplified in the somatic filtering python script. Do not run multiple tumor genomes together since overlapping somatic SVs might have different coordinates in different tumor genomes. One exception is if these tumor samples are from the same patient (e.g., different time points or locations of a given tumor). The general setup should be tumor.bam + control.bam(s).

* Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?  
Delly expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them (e.g. the longest split alignment) is a primary alignment whereas all others are marked as secondary or supplementary (flag 0x0100 or flag 0x0800).

* Usage/discussion mailing list?         
There is a delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users).

* Docker support?            
There is a dockerized delly available [here](https://registry.hub.docker.com/u/trausch/delly/).


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.  
[Delly: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)  
Bioinformatics 2012 28: i333-i339.
