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

Delly is an integrated structural variant (SV) and copy-number variant (CNV) prediction method that can discover, genotype and visualize deletions, insertions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read whole-genome sequencing data. It uses paired-ends, split-reads and read-depth to discover and genotype SVs in the genome.

# Installing Delly

Delly is available as a [pre-compiled binary](https://github.com/dellytools/delly/releases/), a [singularity container (SIF file)](https://github.com/dellytools/delly/releases/), a [docker container](https://hub.docker.com/r/dellytools/delly/) or via [Bioconda](https://anaconda.org/bioconda/delly). You can also build Delly from source: 

`git clone --recursive https://github.com/dellytools/delly.git && cd delly/ && make all`

# Running Delly

Delly needs a sorted, indexed and duplicate marked BAM or CRAM file for every input sample. An indexed reference genome is required to identify split-reads and to decode CRAM files.

## Short-read SV discovery

`delly sr -g ref.fa input.bam > delly.vcf`

You can also redirect the output to a [BCF](http://samtools.github.io/bcftools/) file.

`delly sr -o delly.bcf -g ref.fa input.bam`

`bcftools view delly.bcf > delly.vcf`


## Long-read SV discovery

For Oxford Nanopore sequencing data:

`delly lr -y ont -o delly.bcf -g hg38.fa input.bam`

For PacBio sequencing data:

`delly lr -y pb -o delly.bcf -g hg38.fa input.bam`

## Assembly-based SV discovery

For a haplotype-resolved assembly:

`delly asm -y diploid -g hg38.fa s1.h1.bam s1.h2.bam s2.h1.bam s2.h2.bam ...`

For a so-called squashed assembly:

`delly asm -y squashed -g hg38.fa s1.bam s2.bam ...`

## Somatic SV calling

Somatic SV calling is available for short-reads (subcommand: sr) and long-reads (subcommand: lr).

* At least one tumor sample and a matched control sample are required for SV discovery

`delly [sr|lr] -o t1.bcf -g hg38.fa tumor1.bam control1.bam`

* Somatic filtering requires a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

`delly filter -f somatic -o t1.somatic.bcf -s samples.tsv t1.bcf`

* You can also use a larger panel of normal for somatic SV filtering

`delly [sr|lr] -o t1.bcf -g hg38.fa tumor1.bam control1.bam ... controlN.bam`


## Germline SV calling

Germline SV calling is available for short-reads (subcommand: sr) and long-reads (subcommand: lr).

* SV discovery is done by sample for high-coverage genomes

`delly [sr|lr] -g hg38.fa -o s1.bcf sample1.bam`

* Merge SV sites into an SV site list 

`delly merge -o sites.bcf s1.bcf s2.bcf ... sN.bcf`

* Genotype this merged SV site list across all samples in parallel

`delly [sr|lr] -g hg38.fa -v sites.bcf -o s1.geno.bcf s1.bam`

`delly [sr|lr] -g hg38.fa -v sites.bcf -o sN.geno.bcf sN.bam`

* Merge all genotyped samples to get a single VCF/BCF using bcftools. This can be done in chunks if necessary.

`bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf`

* Apply the germline SV filter

`delly filter -f germline -o germline.bcf merged.bcf`


## CNV calling

`delly cnv` is a whole-genome read-depth method with GC-bias and mappability correction. Whole-exome or targeted data are not supported. You can generate read-depth profiles with delly using

`delly cnv -g hg38.fa -c out.cov.gz -o out.bcf -u out.seg.bed input.bam`

The output file `out.cov.gz` and the segmentation `out.seg.bed` can be plotted using [R](https://www.r-project.org/)

`Rscript R/rd.R out.cov.gz out.seg.bed`

Instead of the segmentation, you can also visualize the CNV calls.

`bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" out.bcf > seg.bed`

`Rscript R/rd.R out.cov.gz seg.bed`


## Germline CNV calling

Germline CNV calling can be done using short- or long-reads using this workflow.

* Call CNVs for each sample

`delly cnv -o c1.bcf -g hg38.fa input.bam`

* Merge CNVs into a unified CNV site list

`delly merge -e -o sites.bcf c1.bcf c2.bcf ... cN.bcf`

* Genotype CNVs for each sample

`delly cnv -v sites.bcf -g hg38.fa -o geno1.bcf input.bam`

* Merge genotypes using [bcftools](https://github.com/samtools/bcftools)

`bcftools merge -m id -O b -o merged.bcf geno1.bcf ... genoN.bcf`

* Apply the germline CNV filter

`delly filter -f germline -o germline.bcf merged.bcf`

* Optional: Plot copy-number distributions for large number of samples (>>100)

`bcftools query -f "%ID[\t%RDCN]\n" germline.bcf > plot.tsv`

`Rscript R/cnv.R plot.tsv`

## Somatic copy-number alterations (SCNAs)

* CNVs are first called on the tumor genome. Depending on the coverage and tumor purity you can adapt the minimum CNV size (`-z`) and the sample purity (`-p`) which affect the sensitivity of SCNA detection.

`delly cnv -o tumor.bcf -c tumor.cov.gz -g hg38.fa tumor.bam`

* Then these CNVs are genotyped in the control sample.

`delly cnv -v tumor.bcf -o control.bcf -g hg38.fa control.bam`

* You can then merge both files using [bcftools](https://github.com/samtools/bcftools).

`bcftools merge -m id -O b -o tumor_control.bcf tumor.bcf control.bcf`

* Somatic filtering requires a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

`delly filter -p -f somatic -o somatic.bcf -s samples.tsv tumor_control.bcf`

* Optional: Plot the SCNAs using bcftools and R.

`bcftools query -s tumor -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" somatic.bcf > segmentation.bed`

`Rscript R/rd.R tumor.cov.gz segmentation.bed`


## Examples

Some small examples are included for short-read, long-read and copy-number variant calling.

`delly sr -g example/ref.fa -o sr.bcf example/sr.bam`

`delly lr -g example/ref.fa -o lr.bcf example/lr.bam`

`delly cnv -g example/ref.fa -c out.cov.gz -o cnv.bcf example/sr.bam`

More in-depth tutorials for SV calling are available here:

* Short-read SV calling: [https://github.com/tobiasrausch/vc](https://github.com/tobiasrausch/vc)

* Long-read SV calling: [https://github.com/tobiasrausch/sv](https://github.com/tobiasrausch/sv)


## Alternate alignments for genome graphs

Instead of providing only one input alignment, delly supports now multiple alternate alignments on different linear reference genomes using [minimap2](https://github.com/lh3/minimap2) or pan-genome graphs using [minigraph](https://github.com/lh3/minigraph).

```
minimap2 -ax map-pb -L chm13.fa sample.fq.gz
minigraph --vc -cx lr pangenome.gfa.gz sample.fq.gz
```

If the above alignment files are then stored as `sample.chm13.bam` and `sample.gaf.gz` you can use a simple tab-delimited config file for all alternate alignments with delly.

`cat align.config`

```
sample.chm13.bam   chm13.fa
sample.gaf.gz   pangenome.gfa.gz
```

`delly lr -y pb -o delly.bcf -g hg38.fa -l align.config sample.hg38.bam`

Structural variants are still reported with respect to GRCh38 coordinates but the output will only contain SVs that are not present in any of the alternate alignments. Please note that many pangenome graphs contain only SVs greater 50bp so you need to filter the above delly output to match the size range.

`bcftools view -i '(QUAL>=300) && ( ((SVTYPE=="INS") && (INFO/SVLEN>50)) || (SVTYPE=="BND") || ((INFO/END - POS)>50) )' delly.bcf`

Please note that for inter-chromosomal translocations, delly uses `INFO/CHR2` for the second chromosome. You can convert an inter-chromosomal translocation to the two-record breakend format using:

`python scripts/delly2bnd.py -v delly.bcf -r hg38.fa -o delly.bnd.bcf`


# FAQ

* Visualization of SVs      
You may want to try out [wally](https://github.com/tobiasrausch/wally) to plot candidate structural variants. The paired-end coloring is explained in [wally's README](https://github.com/tobiasrausch/wally#paired-end-view) file.

* Can Delly be used on a non-diploid genome?  
The SV site discovery works for any ploidy. However, Delly's genotyping model assumes a ploidy of 2. The CNV calling allows to set the baseline ploidy on the command-line.

* Delly is running too slowly what can I do?    
For short-reads, you should exclude telomere and centromere regions and also all unplaced contigs (`-x` command-line option). In addition, you can filter input reads more stringently using -q 20 and -s 15. Lastly, `-z` can be set to 5 for high-coverage data.

* What pre-processing of BAM/CRAM files is required?
BAM/CRAM files need to be sorted, indexed and ideally duplicate marked.

* Usage/discussion mailing list?         
There is a delly discussion group [delly-users](http://groups.google.com/d/forum/delly-users).

# Citation

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.      
DELLY: structural variant discovery by integrated paired-end and split-read analysis.     
Bioinformatics. 2012 Sep 15;28(18):i333-i339.       
[https://doi.org/10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)

# License

Delly is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/dellytools/delly/blob/main/LICENSE) file for more details.

# Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for all genomic alignment and variant processing. [Boost](https://www.boost.org/) for various data structures and algorithms. [Claude](https://claude.com/) for bug fixes, performance improvements and code suggestions. [arfer](https://github.com/ekg/arfer) for annotating SVs and [Edlib](https://github.com/Martinsos/edlib) for pairwise alignments.
