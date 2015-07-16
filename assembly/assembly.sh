#!/bin/bash

if [ $# -lt 3 ]
then
    echo "**********************************************************************"
    echo "Beta breakpoint assembly pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo "Requirement 1: Paired-end data!"
    echo "Requirement 2: SPADES assembler!"
    echo ""
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <genome.fa> <variants.vcf> <sampleToBam.list> [<unmappedToBam.list>]"
    echo ""
    exit -1
fi


# check samtools. Not cool that we need samtools but the whole assembly script is gonna be changed in the future anyways.
USE_SAMTOOLS=1
samtools 2> /dev/null
if [ $? -gt 1 ]; then echo "In order to automatically generate reference files for maze, please provide samtools faidx"; USE_SAMTOOLS=0; fi

SPADES=/g/solexa/home/build.big-al/SPAdes-3.5.0-Linux/bin/spades.py

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
GENOME=${1}
VARIANTS=${2}
BAMLIST=${3}
if [ "$#" -gt 3 ]
then
    UMAPPED="-u "${4}
else
    UMAPPED=""
fi

# Get a VCF identifier
SID=`echo ${VARIANTS} | sed 's/^.*\///' | sed 's/\..*$//'`

# Extract assembly reads
python ${BASEDIR}/extractAssemblyReads.py -v ${VARIANTS} -s ${BAMLIST} ${UMAPPED} 2> coords.table

# Run the assembly
for FQ1 in ${SID}.*.1.fastq
do
    FQ2=`echo ${FQ1} | sed 's/.1.fastq/.2.fastq/'`
    SVID=`echo ${FQ1} | sed 's/.[12].fastq$//' | sed 's/^.*\.//'`
    mkdir assembly.${SID}.${SVID}
    python ${SPADES} --pe1-1 ${FQ1} --pe1-2 ${FQ2} -o assembly.${SID}.${SVID}
    if [ -f assembly.${SID}.${SVID}/contigs.fasta ]
    then
	cat assembly.${SID}.${SVID}/contigs.fasta | gzip -c > contigs.${SID}.${SVID}.fasta.gz
    if [ $USE_SAMTOOLS -ge 1 ]; 
	then
		echo "grep $FQ1 coords.table" $(grep $FQ1 coords.table)
		grep $FQ1 coords.table | cut -f2 | xargs samtools faidx $GENOME | gzip -c > contigs.${SID}.${SVID}.reference.fasta.gz; fi
    fi
    rm -rf assembly.${SID}.${SVID} ${FQ1} ${FQ2}
done
rm coords.table
