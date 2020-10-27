#!/bin/bash

module load bcftools
module load bedtools

out_dir="/path/to/output/files/"
out_dir_final="/path/to/final/output/files/"
gtf="[gencode_gtf].gtf"
file=${1}
sample=${2}
window=10000 # define window around gene

# Create directory on node scratch if does not already exist
mkdir -p ${out_dir}

# 1. Extract allele frequency column for variants +/- 10kb around genes in GTF and save to file
bcftools view -c1[:minor] -a -s $sample $file | bcftools view -f PASS - | \
 bedtools window -header -w $window -a stdin -b $gtf | bcftools query -f "%gnomAD_AF\t%cadd_phred\t%cadd_raw\n" > ${out_dir}/${file/.vcf.gz/_${sample}_intermediate_af.bed} 

wait

# 2. Extract gene name for variants +/- 10kb around genes in GTF and save to file
bcftools view -c1[:minor] -a -s $sample $file | bcftools view -f PASS - | \
  bedtools window -w $window -a stdin -b $gtf | awk -v OFS='\t' '{gsub(/\"|;/,"",$19)}{print $19,$1,$2}' > ${out_dir}/${file/.vcf.gz/_${sample}_intermediate_gene_name.bed} 

wait

# 3. Concatenate 1. and 2.
paste ${out_dir}/${file/.vcf.gz/_${sample}_intermediate_gene_name.bed} ${out_dir}/${file/.vcf.gz/_${sample}_intermediate_af.bed} > \
	${out_dir}/${file/.vcf.gz/_${sample}_10kb-overlap.bed}

wait

# 4. Copy 3. to final directory
mv ${out_dir}/${file/.vcf.gz/_${sample}_10kb-overlap.bed} ${out_dir_final}

wait

# 4. Clean up intermediate files
rm ${out_dir}/${file/.vcf.gz/_${sample}*.bed}
