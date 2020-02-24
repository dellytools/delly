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

# Benchmark SV calling
echo -e "svtype\tnOption\tmode\tcoverage\taccuracy\treadlen\trecall\tprecision\tf1\tgtconc" > summary.stats.tsv
for SVT in INS DEL
do
    for MODE in ONT PB
    do
	if [ ${MODE} == "ONT" ]
	then
	    SUBINDEL="45:25:30"
	else
	    SUBINDEL="15:50:35"
	fi
	#for COV in 10 15 20  # haplotype coverage
	for COV in 15
	do
	    #for ACC in 0.95 0.9 0.85
	    for ACC in 0.9
	    do
		#for LEN in 1000 5000 9000
		for LEN in 5000
		do
		    # Simulate anew
		    if [ ! -d sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN} ]
		    then
			mkdir sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}
		
			# Simulate haplotypes
			if [ ${SVT} == "INS" ]
			then
			    bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%ALT\n" HG002_SVs_Tier1_v0.6.vcf.gz  | grep -w "${SVT}" | grep "^18" | awk 'length($4)>50 && length($4)<10000' | sed 's/INS/insertion/' | awk '{ if ($2-OLD>10000) {print $1"\t"$2"\t"($2+1)"\t"$3"\t"$4"\t0\t"int(rand()+1.5);}; OLD=$2;}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap1.${SVT}.bed
			    cat sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap1.${SVT}.bed | awk '{print $1"\t"$2"\t"($2+length($5))"\t"$4"\tNone\t"$6"\t"$7;}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed
			else
			    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n" HG002_SVs_Tier1_v0.6.vcf.gz | grep -w "${SVT}" | grep "^18" | awk '$3-$2>=30' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/dels.bed
			    bedtools cluster -i sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | cut -f 5 | sort | uniq -u > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/fetchDels
			    bedtools cluster -i sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/dels.bed  | sed 's/DEL\t/DEL\tCluster/' | grep --color -w -Ff sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/fetchDels | sed 's/DEL\tCluster[0-9]*/deletion/' | awk '{print $0"\tNone\t0\t"int(rand()+1.5);}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed
			    cp sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap1.${SVT}.bed
			fi
			cat sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap1.${SVT}.bed | grep -P "\t2$" > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap2.${SVT}.bed
			cat sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$7;}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.igv.bed
			VISOR HACk -g chr18.fa -bed sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap1.${SVT}.bed sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hap2.${SVT}.bed -o sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hack
			
			# Draw reads
			cat chr18.fa.fai  | awk '{print $1"\t0\t"$2"\t100.0\t100.0";}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/simulate.bed
			VISOR LASeR --addprefix -c ${COV} -a ${ACC} -l ${LEN} -r ${SUBINDEL} -g chr18.fa -s sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/hack/ -bed sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/simulate.bed -o sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/laser
		    fi
		    
		    # Call SVs using various parameters
		    for N in 75 150 300
		    do
			for C in 25 50 100
			do
			    for E in 0.5 0.6 0.7
			    do
				for F in 50 100 250
				do
				    for A in 0.8 0.85 0.9
				    do
					for S in "3,-2,-3,-1" "5,-2,-5,-2" "5,-2,-3,-1"
					do
					    ../bin/dellyLR call -c ${C} -n ${N} -e ${E} -f ${F} -a ${A} -s ${S} -o sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bcf -g chr18.fa sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/laser/sim.srt.bam
					    if [ ${SVT} == "INS" ]
					    then
						bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\t%ID[\t%GT]\n' sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bcf | awk '{print $1"\t"$2"\t"($2+$3)"\t"$4"\t"$5;}' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/called.bed
					    else
						bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n' sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bcf > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/called.bed
					    fi
					    bedtools intersect -a sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed -b sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/called.bed  -wao | awk '$8!="." && ((($9-$2) + ($3-$10))>-50) && ((($9-$2) + ($3-$10))<50)' > sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/true_positives.bed
					    TOTTRUTH=`cat sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/${SVT}.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
					    TOTCALLED=`cat sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/called.bed | cut -f 1-4 | sort | uniq | wc -l | cut -f 1`
					    TP=`cut -f 1-4 sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/true_positives.bed | sort | uniq | wc -l | cut -f 1`
					    RECALL=`echo "${TP} / ${TOTTRUTH}" | bc -l`
					    PREC=`echo "${TP} / ${TOTCALLED}" | bc -l`
					    F1=`echo "2 * (${RECALL} * ${PREC}) / (${RECALL} + ${PREC})" | bc -l`
					    CONCHOM=`cut -f 7,12 sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t2\t1/1" | cut -f 1`
					    CONCHET=`cut -f 7,12 sim_svt${SVT}_${MODE}_cov${COV}_acc${ACC}_len${LEN}/true_positives.bed | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | grep -P "\t1\t0/1" | cut -f 1`
					    GTCONC=`echo "(${CONCHOM} + ${CONCHET}) / ${TOTCALLED}" | bc -l`
					    echo -e "${SVT}\t${C}\t${N}\t${E}\t${F}\t${A}\t${S}\t${MODE}\t${COV}\t${ACC}\t${LEN}\t${RECALL}\t${PREC}\t${F1}\t${GTCONC}" >> summary.stats.tsv
					done
				    done
				done
			    done
			done
		    done
		done
	    done
	done
    done
done

# Done
conda deactivate
