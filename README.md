# i2QTL outlier analysis

Scripts required for i2QTL outlier calling pipeline

**Resources required (versions indicated as used in paper, if relevant):**
* gnomAD (r2.0.2)
* CADD (1.3)
* SAMtools (1.11)
* bedtools 

## Overview  

### 1. Genetic data processing
* 1.1 VCF annotation 
* 1.2 Extract variants by sample

### 2. RNA-seq data processing
* 2.1 Data correction

### 3. Outlier calling
* 3.1 Count rare variants linked to outliers/non-outliers

## 1.1 Variant annotation
* Add gnomAD allele frequency and CADD score for each variant in VCF
* Assumes VCFs are split by chromosome and variant type (i.e. SNP, indel)

```
gnomAD_dir="[location_of_gnomAD_reference_vcf]"
cadd_dir="[location_of_cadd_reference_vcf]"

for var in snp indel; do
  for chr in `seq 1 22`; do
    # Create vcfanno annotation file
    echo '[[annotation]]' > conf${chr}.toml 
    echo 'file="'${gnomAD_dir}/gnomad.genomes.r2.0.2.sites.chr${chr}.vcf.gz'"' >> conf${chr}.toml 
    echo 'fields=["AF"]' >> conf${chr}.toml 
    echo 'names=["gnomAD_AF"]' >> conf${chr}.toml 
    echo 'ops=["self"]' >> conf${chr}.toml
    echo "">> conf${chr}.toml 
    echo "">> conf${chr}.toml
    echo '[[annotation]]'>> conf${chr}.toml
    echo 'file="'$cadd_dir/CADD.v1.3.SNP.INDEL.sorted.vcf.gz'"'>> conf${chr}.toml
    echo 'fields=["phred", "raw"]'>> conf${chr}.toml
    echo 'names=["cadd_phred", "cadd_raw"]'>> conf${chr}.toml
    echo 'ops=["mean", "mean"]'>> conf${chr}.toml

    vcfanno conf${chr}.toml ${in_out_dir}/IPSCORE_HipSci_1kG.${chr}.*${var}.vcf.gz | \
        gzip -c > IPSCORE_HipSci_1kG.${chr}_${var}_gnomadAF_CADD.vcf.gz
  done
done

```

## 1.2 Extract variants by sample
* Extract variants marked as PASS and falling within 10Kb of any protein-coding or long non-coding RNA gene
* Write separate file per sample
* Concatenate output across chromosome

```
for var in snp indel; do
  for file in IPSCORE_HipSci_1kG.*_${var}_gnomadAF_CADD.vcf.gz; do 
    while read sample; do
      sbatch intersect_gtf_vcf.sh $file $sample
    done < ../GTE_WGS_combined_sample_name.txt
  done
done

parallel --verbose -j 50 --xapply -a [list_of_sample_names].txt \
  'cat IPSCORE_HipSci_1kG.*_snp_gnomadAF_CADD_{1}*.bed > ../IPSCORE_HipSci_1kG_gnomadAF_CADD_snp.{1}.bed'

parallel --verbose -j 50 --xapply -a [list_of_sample_names].txt \
  'cat IPSCORE_HipSci_1kG.*_indel_gnomadAF_CADD_{1}*.bed > ../IPSCORE_HipSci_1kG_gnomadAF_CADD_indel.{1}.bed'

```

## 2.1 Data correction
* Read uncorrected RNA-seq counts, filter genes to protein-coding and long non-coding RNA only
* Run PEER, remove global expression outlier samples

```
[see scripts/gene_exp_filter_zscore.r]
```

## 3. Outlier calling
* Intersect genetic data with gene expression outlier/non-outliers, per individual

```
exp_dat="[iPSC_corrected_zscore].txt" # Corrected RNA-seq matrix (output from 2.1)
iterator="iterator.txt" # Define different Z-score, MAF, and CADD thresholds (comma separated)

for var in snp indel; do
  while read sample; do
    sbatch outlier_rare_var.R $sample $var $exp_dat $iterator
  done < [list_of_sample_names].txt
done
```




