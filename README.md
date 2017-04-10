<p align="center">
  <a href="https://academic.oup.com/bioinformatics/article/28/18/i333/245403/DELLY-structural-variant-discovery-by-integrated">
    <img height="150" src="https://raw.githubusercontent.com/dellytools/assets/master/delly-logo/delly-logo-539x600.png">
  </a>
  <h1 align="center">Delly2</h1>
</p>

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/delly/README.html)
[![Build Status](https://travis-ci.org/dellytools/delly.svg?branch=master)](https://travis-ci.org/dellytools/delly)
[![Docker Automated buil](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg?style=flat-square)](https://hub.docker.com/r/dellytools/delly/)

Delly2 is an integrated structural variant prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome. Structural variants can be visualized using [Delly-maze](https://github.com/dellytools/maze) and [Delly-suave](https://github.com/dellytools/suave).


Installing Delly2
-----------------

The easiest way to get Delly2 is to download a statically linked binary from the [Delly github release page](https://github.com/dellytools/delly/releases/).
Alternatively, you can build Delly2 from source. Delly2 dependencies are included as submodules so you need to do a recursive clone. 

`git clone --recursive https://github.com/dellytools/delly.git`

`cd delly/`

`make all`

There is a Delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users) for usage and installation questions and a dockerized [delly](https://hub.docker.com/r/dellytools/delly/).



Delly2 multi-threading mode
---------------------------
Delly2 supports parallel computing using the OpenMP API (www.openmp.org).

`make PARALLEL=1 -B src/delly`

There is also a statically linked, multi-threaded binary for Linux 64-bit available under [releases](https://github.com/dellytools/delly/releases/).


You can set the number of threads using the environment variable OMP_NUM_THREADS.

`export OMP_NUM_THREADS=3`

Delly2 primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. 

Running Delly2
--------------

Delly2 needs a sorted and indexed bam file for every sample and the reference genome to identify split-reads. The output is in [bcf](http://samtools.github.io/bcftools/) format with a csi index. Prior duplicate marking is recommended. The SV type can be DEL, DUP, INV, TRA, or INS for deletions, tandem duplications, inversions, translocations and small insertions, respectively. Delly2 supports germline and somatic SV discovery, genotyping and filtering. Because of that, Delly2 has been modularized and common workflows for germline and somatic SV calling are outlined below. If you do need VCF output you need a recent version of [bcftools](http://samtools.github.io/bcftools/) (included as a submodule in Delly2) for file conversion.

`./delly/src/bcftools/bcftools view delly.bcf > delly.vcf`


Somatic SV calling
------------------

* At least one tumor sample and a matched control sample are required. All tumor/control pairs are run separately for SV discovery:

`delly call -t DEL -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam`

* Somatic pre-filtering of every tumor/control pair using a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

`delly filter -t DEL -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf`

* Re-genotype somatic sites across a larger panel of control samples to efficiently filter false postives and germline SVs. For performance reasons, this can be run in parallel for each sample (see germline SV calling) and/or directly on a combined pre-filtered somatic site list from multiple tumor/control pairs.

`delly call -t DEL -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam`

* Post-filter for somatic SVs using all control samples.

`delly filter -t DEL -f somatic -o t1.somatic.bcf -s samples.tsv geno.bcf`



Germline SV calling
-------------------

* SV calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision

`delly call -t DEL -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam`

* Merge SV sites into a unified site list 

`delly merge -t DEL -m 500 -n 1000000 -o del.bcf -b 500 -r 0.5 s1.bcf s2.bcf ... sN.bcf`

* Re-genotype merged SV site list across all samples. This can be run in parallel for each sample.

`delly call -t DEL -g hg19.fa -v del.bcf -o s1.geno.bcf -x hg19.excl s1.bam`

`delly call -t DEL -g hg19.fa -v del.bcf -o sN.geno.bcf -x hg19.excl sN.bam`

* Merge all re-genotyped samples to get a single VCF/BCF using bcftools merge

`bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf`

* Apply the germline SV filter

`delly filter -t DEL -f germline -o germline.bcf merged.bcf`

FAQ
---
* What is the smallest SV size Delly can call?  
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, Delly starts to call reliable SVs >=300bp. Delly2 also supports calling of small InDels using soft-clipped reads only. In this mode the smallest SV size called is 15bp.

* Can Delly be used on a non-diploid genome?  
Yes and no. The SV site discovery works for any ploidy. However, Delly's genotyping model assumes hom. reference, het. and hom. alternative.

* How do I run Delly if I have multiple different libraries/bam files for a single sample?    
Merge these bams using tools such as [Picard](http://broadinstitute.github.io/picard/) and tag each library with a unique ReadGroup. 

* Delly is running too slowly what can I do?    
Exclude telomere and centromere regions. Delly2 ships with such an exclude list for human and mouse samples. In addition, you can filter input reads more stringently using -q 20 and -s 15.

* Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?  
Delly expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them (e.g. the longest split alignment) is a primary alignment whereas all others are marked as secondary or supplementary (flag 0x0100 or flag 0x0800). This is the default for bwa mem.

* What pre-processing of bam files is required?    
Bam files need to be sorted and index. If multiple libraries are present for a single sample (e.g., a long-insert mate-pair library and a short-insert paired-end library) these need to be merged in a single bam file with unique ReadGroup tags. A prior marking of duplicates is recommended.

* Usage/discussion mailing list?         
There is a delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users).

* Docker support?            
There is a dockerized delly available [here](https://hub.docker.com/r/dellytools/delly/).


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.  
[Delly: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)  
Bioinformatics 2012 28: i333-i339.


License
-------
Delly is distributed under the GPLv3. Consult the accompanying [LICENSE](https://github.com/dellytools/delly/blob/master/LICENSE) file for more details.
