#!/bin/bash

DELLY=`pwd`/../src/delly
TYPE=DEL
OUT=sv.vcf

if [ -f ${DELLY} ]
then
    # delly
    ${DELLY} -t ${TYPE} -o ${OUT} -g ${TYPE}.fa ${TYPE}.bam
    DC=$?
    if [ ${DC} -ne 0 ]
    then
	echo "Test failed!"
	exit ${DC}
    fi
    
    # diff
    diff <(cut -f 1,2,4- ${OUT} | grep -v "^#" | sed 's/SVMETHOD=[^;]*;//' | sort | uniq) <(cut -f 1,2,4- ${TYPE}.vcf | grep -v "^#" | sed 's/SVMETHOD=[^;]*;//' | sort | uniq)
    DC=$?
    rm ${OUT}
    if [ ${DC} -ne 0 ]
    then
	echo "Test failed!"
	exit ${DC}
    fi

    echo "Test run was successful!"
fi
