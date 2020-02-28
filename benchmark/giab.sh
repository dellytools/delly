#!/bin/bash

# Base directory
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Download GIAB data set
if [ ! -f hs37d5.fa.gz ]
then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    samtools faidx hs37d5.fa.gz
fi
if [ ! -f HG002_SVs_Tier1_v0.6.bed ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
fi
if [ ! -f HG002_SVs_Tier1_v0.6.vcf.gz ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
fi
if [ ! -f ultra-long-ont_hs37d5_phased.bam ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/final/ultra-long-ont_hs37d5_phased.bam
    samtools index ultra-long-ont_hs37d5_phased.bam
fi
if [ ! -f HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
    samtools index HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
fi
if [ ! -f HG002_PacBio_GRCh37.bam ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh37.bam
    samtools index HG002_PacBio_GRCh37.bam
fi

# Install tools
make test
if [ $? -ne 0 ]
then
    echo ""
    echo "Please install dependencies first using 'make all'"
    exit -1;
fi

# Activate environment
export PATH=${BASEDIR}/bin/bin/:${PATH}
source activate truvari

# Delly for long reads
rm -f delly.vcf*
../bin/dellyLR call -g hs37d5.fa.gz ultra-long-ont_hs37d5_phased.bam
bcftools view -i '%QUAL>400' sv.bcf | grep -v 'INV\|BND\|DUP' > delly.vcf
bgzip delly.vcf
tabix delly.vcf.gz

# truvari
rm -rf reportDelly
truvari --includebed HG002_SVs_Tier1_v0.6.bed --giabreport --passonly --no-ref a -p 0.00 -f hs37d5.fa.gz -b HG002_SVs_Tier1_v0.6.vcf.gz -c delly.vcf.gz -o reportDelly

# Done
source deactivate
