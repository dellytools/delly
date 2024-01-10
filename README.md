<p align="center">
  <a href="https://academic.oup.com/bioinformatics/article/28/18/i333/245403/DELLY-structural-variant-discovery-by-integrated">
    <img height="150" src="https://raw.githubusercontent.com/dellytools/assets/main/delly-logo/delly-logo-539x600.png">
  </a>
  <h1 align="center">Delly</h1>
</p>

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/delly/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/delly/badges/downloads.svg)](https://anaconda.org/bioconda/delly)
[![C/C++ CI](https://github.com/dellytools/delly/workflows/C/C++%20CI/badge.svg)](https://github.com/dellytools/delly/actions)
[![Docker CI](https://github.com/dellytools/delly/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/dellytools/delly/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/dellytools/delly/blob/main/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/dellytools/delly.svg)](https://github.com/dellytools/delly/releases)

Delly is an integrated structural variant (SV) prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.

# Installing Delly

Delly is available as a [statically linked binary](https://github.com/dellytools/delly/releases/), a [singularity container (SIF file)](https://github.com/dellytools/delly/releases/), a [docker container](https://hub.docker.com/r/dellytools/delly/) or via [Bioconda](https://anaconda.org/bioconda/delly). You can also build Delly from source using a recursive clone and make. 

`git clone --recursive https://github.com/dellytools/delly.git`

`cd delly/`

`make all`

There is a Delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users) for usage and installation questions.


# Delly multi-threading mode

Delly supports parallel computing using the OpenMP API (www.openmp.org).

`make PARALLEL=1 src/delly`

You can set the number of threads using the environment variable OMP_NUM_THREADS.

`export OMP_NUM_THREADS=2`

Delly primarily parallelizes on the sample level. Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples. 


# Running Delly

Delly needs a sorted, indexed and duplicate marked bam file for every input sample.
An indexed reference genome is required to identify split-reads.
Common workflows for germline and somatic SV calling are outlined below.

`delly call -g hg19.fa input.bam > delly.vcf`

You can also specify an output file in [BCF](http://samtools.github.io/bcftools/) format.

`delly call -o delly.bcf -g hg19.fa input.bam`

`bcftools view delly.bcf > delly.vcf`


Example
-------

A small example is included for short-read, long-read and copy-number variant calling.

`delly call -g example/ref.fa example/sr.bam > sr.vcf`

`delly lr -g example/ref.fa example/lr.bam > lr.vcf`

`delly cnv -g example/ref.fa -m example/map.fa.gz example/sr.bam > cnv.vcf`


Somatic SV calling
------------------

* At least one tumor sample and a matched control sample are required for SV discovery

`delly call -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam`

* Somatic pre-filtering requires a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

`delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf`

* Genotype pre-filtered somatic sites across a larger panel of control samples to efficiently filter false postives and germline SVs. For performance reasons, this can be run in parallel for each sample of the control panel and you may want to combine multiple pre-filtered somatic site lists from multiple tumor samples.

`delly call -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam`

* Post-filter for somatic SVs using all control samples.

`delly filter -f somatic -o t1.somatic.bcf -s samples.tsv geno.bcf`



Germline SV calling
-------------------

* SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes

`delly call -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam`

* Merge SV sites into a unified site list 

`delly merge -o sites.bcf s1.bcf s2.bcf ... sN.bcf`

* Genotype this merged SV site list across all samples. This can be run in parallel for each sample.

`delly call -g hg19.fa -v sites.bcf -o s1.geno.bcf -x hg19.excl s1.bam`

`delly call -g hg19.fa -v sites.bcf -o sN.geno.bcf -x hg19.excl sN.bam`

* Merge all genotyped samples to get a single VCF/BCF using bcftools merge

`bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf`

* Apply the germline SV filter which requires at least 20 unrelated samples

`delly filter -f germline -o germline.bcf merged.bcf`


Delly for long reads from PacBio or ONT
---------------------------------------

Delly also supports long-reads for SV discovery.

`delly lr -y ont -o delly.bcf -g hg19.fa input.bam`

`delly lr -y pb -o delly.bcf -g hg19.fa input.bam`


Read-depth profiles
-------------------

You can generate read-depth profiles with delly. This requires a mappability map which can be downloaded here:

[Mappability Maps](https://gear-genomics.embl.de/data/delly/)

The command to count reads in 10kbp mappable windows and normalize the coverage is:

`delly cnv -a -g hg19.fa -m hg19.map input.bam`

The output file `out.cov.gz` can be plotted using [R](https://www.r-project.org/) to generate normalized copy-number profiles:

`Rscript R/rd.R out.cov.gz`


Copy-number segmentation
------------------------

Read-depth profiles can also be segmented at the same time.

`delly cnv -a -u -g hg19.fa -m hg19.map input.bam`

The segmentation is in VCF format but you can extract a BED-like file using bcftools.

`bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" cnv.bcf > segmentation.bed`

Plotting:

`Rscript R/rd.R out.cov.gz segmentation.bed`

Germline CNV calling
--------------------

Delly uses GC and mappability fragment correction to call CNVs. This requires a [mappability map](https://gear-genomics.embl.de/data/delly/).

* Call CNVs for each sample and optionally refine breakpoints using delly SV calls

`delly cnv -o c1.bcf -g hg19.fa -m hg19.map -l delly.sv.bcf input.bam`

* Merge CNVs into a unified site list

`delly merge -e -p -o sites.bcf -m 1000 -n 100000 c1.bcf c2.bcf ... cN.bcf`

* Genotype CNVs for each sample

`delly cnv -u -v sites.bcf -g hg19.fa -m hg19.map -o geno1.bcf input.bam`

* Merge genotypes using [bcftools](https://github.com/samtools/bcftools)

`bcftools merge -m id -O b -o merged.bcf geno1.bcf ... genoN.bcf`

* Filter for germline CNVs

`delly classify -f germline -o filtered.bcf merged.bcf`

* Optional: Plot copy-number distribution for large number of samples (>>100)

`bcftools query -f "%ID[\t%RDCN]\n" filtered.bcf > plot.tsv`

`Rscript R/cnv.R plot.tsv`


Somatic copy-number alterations (SCNAs)
---------------------------------------

* For somatic copy-number alterations, delly first segments the tumor genome (`-u` is required). Depending on the coverage, tumor purity and heterogeneity you can adapt parameters `-z`, `-t` and `-x` which control the sensitivity of SCNA detection.

`delly cnv -u -z 10000 -o tumor.bcf -c tumor.cov.gz -g hg19.fa -m hg19.map tumor.bam`

* Then these tumor SCNAs are genotyped in the control sample (`-u` is required).

`delly cnv -u -v tumor.bcf -o control.bcf -g hg19.fa -m hg19.map control.bam`

* The VCF IDs are matched between tumor and control. Thus, you can merge both files using [bcftools](https://github.com/samtools/bcftools).

`bcftools merge -m id -O b -o tumor_control.bcf tumor.bcf control.bcf`

* Somatic filtering requires a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

`delly classify -p -f somatic -o somatic.bcf -s samples.tsv tumor_control.bcf`

* Optional: Plot the SCNAs using bcftools and R.

`bcftools query -s tumor -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" somatic.bcf > segmentation.bed`

`Rscript R/rd.R tumor.cov.gz segmentation.bed`


FAQ
---
* Visualization of SVs      
You may want to try out [wally](https://github.com/tobiasrausch/wally) to plot candidate structural variants. The paired-end coloring is explained in [wally's README](https://github.com/tobiasrausch/wally#paired-end-view) file.

* What is the smallest SV size Delly can call?  
For short-reads, this depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, Delly starts to call reliable SVs >=300bp. Delly also supports calling of small InDels using soft-clipped reads only, the smallest SV size called is 15bp. For long-reads, delly calls SVs >=30bp.

* Can Delly be used on a non-diploid genome?  
Yes and no. The SV site discovery works for any ploidy. However, Delly's genotyping model assumes diploidy (hom. reference, het. and hom. alternative). The CNV calling allows to set the baseline ploidy on the command-line.

* Delly is running too slowly what can I do?    
You should exclude telomere and centromere regions and also all unplaced contigs (`-x` command-line option). In addition, you can filter input reads more stringently using -q 20 and -s 15. Lastly, `-z` can be set to 5 for high-coverage data.

* Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?  
Delly expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them is a primary alignment whereas all others are marked as secondary or supplementary. This is the default for bwa, minimap2 and many other aligners.

* What pre-processing of bam files is required?    
Bam files need to be sorted, indexed and ideally duplicate marked.

* Usage/discussion mailing list?         
There is a delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users).

* Docker/Singularity support?            
There is a delly [docker container](https://hub.docker.com/r/dellytools/delly/) and [singularity container (*.sif file)](https://github.com/dellytools/delly/releases) available.

* How can I compute a mappability map?               
A basic mappability map can be built using [dicey](https://github.com/gear-genomics/dicey), [samtools](https://github.com/samtools/samtools) and [bwa](https://github.com/lh3/bwa) with the below commands (as an example for the sacCer3 reference):
```
dicey chop sacCer3.fa
bwa index sacCer3.fa
bwa mem sacCer3.fa read1.fq.gz read2.fq.gz | samtools sort -@ 8 -o srt.bam -
samtools index srt.bam 
dicey mappability2 srt.bam 
gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 
```

* Bioconda support?              
Delly is available via [bioconda](http://bioconda.github.io/recipes/delly/README.html).


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.      
DELLY: structural variant discovery by integrated paired-end and split-read analysis.     
Bioinformatics. 2012 Sep 15;28(18):i333-i339.       
[https://doi.org/10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)

License
-------
Delly is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/dellytools/delly/blob/main/LICENSE) file for more details.
