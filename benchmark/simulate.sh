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
for COV in 10 15 20
do
    for ACC in 0.85 0.9 0.95
    do
	# Simulate anew
	if [ ! -d sim_cov${COV}_acc${ACC} ]
	then
	    mkdir sim_cov${COV}_acc${ACC}

	    # Simulate haplotypes
	    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n" HG002_SVs_Tier1_v0.6.vcf.gz | grep "DEL" | grep "^18" | awk '$3-$2>=30' > sim_cov${COV}_acc${ACC}/dels.bed
	    bedtools cluster -i sim_cov${COV}_acc${ACC}/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | cut -f 5 | sort | uniq -u > sim_cov${COV}_acc${ACC}/fetchDels
	    bedtools cluster -i sim_cov${COV}_acc${ACC}/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | grep --color -w -Ff sim_cov${COV}_acc${ACC}/fetchDels | sed 's/DEL\tCluster[0-9]*/deletion/' | awk '{print $0"\tNone\t0\t"int(rand()+1.5);}' > sim_cov${COV}_acc${ACC}/deletions.bed
	    cat sim_cov${COV}_acc${ACC}/deletions.bed | grep -P "\t2$" > sim_cov${COV}_acc${ACC}/hap2.dels.bed
	    cat sim_cov${COV}_acc${ACC}/deletions.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$7;}' > sim_cov${COV}_acc${ACC}/dels.igv.bed
	    VISOR HACk -g chr18.fa -bed sim_cov${COV}_acc${ACC}/deletions.bed sim_cov${COV}_acc${ACC}/hap2.dels.bed -o sim_cov${COV}_acc${ACC}/hack

	    # Draw reads
	    cat chr18.fa.fai  | awk '{print $1"\t0\t"$2"\t100.0\t100.0";}' > sim_cov${COV}_acc${ACC}/simulate.bed
	    VISOR LASeR --addprefix -c ${COV} -a ${ACC} -g chr18.fa -s sim_cov${COV}_acc${ACC}/hack/ -bed sim_cov${COV}_acc${ACC}/simulate.bed -o sim_cov${COV}_acc${ACC}/laser
	fi

	# Call deletions using various parameters
	for N in 25 50 75 100 200 500
	do
	    ../bin/dellyLR call -n ${N} -o sim_cov${COV}_acc${ACC}/dels.bcf -g chr18.fa sim_cov${COV}_acc${ACC}/laser/sim.srt.bam
	    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n' sim_cov${COV}_acc${ACC}/dels.bcf > sim_cov${COV}_acc${ACC}/called.bed
	    bedtools intersect -a sim_cov${COV}_acc${ACC}/deletions.bed -b sim_cov${COV}_acc${ACC}/called.bed  -wao | awk '$8!="." && ((($9-$2) + ($3-$10))>-50) && ((($9-$2) + ($3-$10))<50)' > sim_cov${COV}_acc${ACC}/true_positives.bed
	    TOTTRUTH=`cat sim_cov${COV}_acc${ACC}/deletions.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
	    TOTCALLED=`cat sim_cov${COV}_acc${ACC}/called.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
	    TP=`cut -f 1-4 sim_cov${COV}_acc${ACC}/true_positives.bed | sort | uniq | wc -l | cut -f 1`
	    RECALL=`echo "${TP} / ${TOTTRUTH}" | bc -l`
	    PREC=`echo "${TP} / ${TOTCALLED}" | bc -l`
	    CONCHOM=`cut -f 7,12 sim_cov${COV}_acc${ACC}/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t2\t1/1" | cut -f 1`
	    CONCHET=`cut -f 7,12 sim_cov${COV}_acc${ACC}/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t1\t0/1" | cut -f 1`
	    GTCONC=`echo "(${CONCHOM} + ${CONCHET}) / ${TOTCALLED}" | bc -l`
	    echo -e "${N}\t${COV}\t${ACC}\t${RECALL}\t${PREC}\t${GTCONC}" >> summary.stats.tsv
	done
    done
done

# Done
conda deactivate
