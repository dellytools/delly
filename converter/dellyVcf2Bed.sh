#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: $0 <delly.vcf.gz>"
    exit -1
fi


vcf-sort ${1} | grep -v "^#" | grep -v "SVTYPE=TRA" | cut -f 1,2,3,8 | sed 's/\t[^\t]*;END=/\t/' | sed 's/;.*$//' | awk '{print $1"\t"$2"\t"$4"\t"$3;}'
vcf-sort ${1} | grep -v "^#" | grep "SVTYPE=TRA" | cut -f 1,2,3,8 | sed 's/\t[^\t]*CHR2=/\t/' | sed 's/;[^\t]*END=/\t/' | sed 's/;.*$//' | awk '{print $1"\t"($2-500)"\t"($2+500)"\t"$3"A\n"$4"\t"($5-500)"\t"($5+500)"\t"$3"B";}'
