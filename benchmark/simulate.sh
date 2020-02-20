#!/bin/bash

# Base directory
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Download GIAB ground truth
if [ ! -f chr18.fa ]
then
    if [ ! -f hs37d5.fa.gz ]
    then
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	samtools faidx hs37d5.fa.gz
    fi
    samtools faidx hs37d5.fa.gz 18 > chr18.fa
    samtools faidx chr18.fa
fi
if [ ! -f HG002_SVs_Tier1_v0.6.vcf.gz ]
then
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
fi

# Install tools
make -W .delly test
if [ $? -ne 0 ]
then
    echo ""
    echo "Please install dependencies first using 'make all'"
    exit -1;
fi

# Activate environment
export PATH=${BASEDIR}/bin/bin/:${PATH}
source activate sv

# Deletions
echo -e "nOption\tcoverage\taccuracy\trecall\tprecision\tgtconc" > summary.stats.tsv
#for COV in 10 15 20
for COV in 10
do
    #for ACC in 0.85 0.9 0.95
    for ACC in 0.95
    do
	# Output directory
	#rm -rf sim/ && mkdir -p sim

	if [ ! sim/deletions.bam ]
	then
	    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n" HG002_SVs_Tier1_v0.6.vcf.gz | grep "DEL" | grep "^18" | awk '$3-$2>=30' > sim/dels.bed
	    bedtools cluster -i sim/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | cut -f 5 | sort | uniq -u > sim/fetchDels
	    bedtools cluster -i sim/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | grep --color -w -Ff sim/fetchDels | sed 's/DEL\tCluster[0-9]*/deletion/' | awk '{print $0"\tNone\t0\t"int(rand()+1.5);}' > sim/deletions.bed
	    cat sim/deletions.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$7;}' > sim/dels.igv.bed
	    VISOR HACk -g chr18.fa -bed sim/deletions.bed -o sim/hack1
	    cat sim/deletions.bed | grep -P "\t2$" > sim/hap2.dels.bed
	    VISOR HACk -g chr18.fa -bed sim/hap2.dels.bed -o sim/hack2
	    
	    # Draw reads
	    cat chr18.fa.fai  | awk '{print $1"\t0\t"$2"\t100.0\t100.0";}' > sim/simulate.bed
	    VISOR LASeR -c ${COV} -a ${ACC} -g chr18.fa -s sim/hack1/ -bed sim/simulate.bed -o sim/laser1
	    VISOR LASeR -c ${COV} -a ${ACC} -g chr18.fa -s sim/hack2/ -bed sim/simulate.bed -o sim/laser2
	    
	    
	    # Fix me, prefix option
	    samtools view -h sim/laser1/sim.srt.bam | sed 's/^S1_/HAP1_S1/' | samtools view -bT chr18.fa > sim/laser1/sim.srt.tmp.bam
	    mv sim/laser1/sim.srt.tmp.bam sim/laser1/sim.srt.bam
	    samtools index sim/laser1/sim.srt.bam
	    samtools view -h sim/laser2/sim.srt.bam | sed 's/^S1_/HAP2_S1/' | samtools view -bT chr18.fa > sim/laser2/sim.srt.tmp.bam
	    mv sim/laser2/sim.srt.tmp.bam sim/laser2/sim.srt.bam
	    samtools index sim/laser2/sim.srt.bam
	    
	    
	    samtools merge sim/deletions.bam sim/laser1/sim.srt.bam sim/laser2/sim.srt.bam
	    samtools index sim/deletions.bam
	fi

	# Call deletions using various parameters
	for N in 50 75 100 200 500
	do
	    ../bin/dellyLR call -n ${N} -o sim/dels.bcf -g chr18.fa sim/deletions.bam
	    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n' sim/dels.bcf > sim/called.bed
	    bedtools intersect -a sim/deletions.bed -b sim/called.bed  -wao | awk '$8!="." && ((($9-$2) + ($3-$10))>-50) && ((($9-$2) + ($3-$10))<50)' > sim/true_positives.bed
	    TOTTRUTH=`cat sim/deletions.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
	    TOTCALLED=`cat sim/called.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
	    TP=`cut -f 1-4 sim/true_positives.bed | sort | uniq | wc -l | cut -f 1`
	    RECALL=`echo "${TP} / ${TOTTRUTH}" | bc -l`
	    PREC=`echo "${TP} / ${TOTCALLED}" | bc -l`
	    CONCHOM=`cut -f 7,12 sim/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t2\t1/1" | cut -f 1`
	    CONCHET=`cut -f 7,12 sim/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t1\t0/1" | cut -f 1`
	    GTCONC=`echo "(${CONCHOM} + ${CONCHET}) / ${TOTCALLED}" | bc -l`
	    echo -e "${N}\t${COV}\t${ACC}\t${RECALL}\t${PREC}\t${GTCONC}" >> summary.stats.tsv
	done

	# Clean-up
	#rm -rf sim/
    done
done

# Done
conda deactivate
