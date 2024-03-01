# Convert Delly inter-chromosomal translocations into breakend format with multiple records

`python3 delly2bnd.py -v delly.bcf -r ref.fa -o out.vcf`

`bcftools sort -O b -o bnd.bcf out.vcf`
