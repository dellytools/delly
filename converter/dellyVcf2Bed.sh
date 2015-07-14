#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: $0 <delly.vcf.gz>"
    exit -1
fi


vcf-sort ${1} | grep -v "^#" | cut -f 1,2,3,8 | sed 's/\t[^\t]*;END=/\t/' | sed 's/;.*$//' | awk '{print $1"\t"$2"\t"$4"\t"$3;}'

