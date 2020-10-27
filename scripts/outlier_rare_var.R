#!/usr/bin/env Rscript

## Find expression outliers, intersect rare variants

library(data.table)

## <----------- MAIN
## Store command line argument to variable
args <- commandArgs(trailingOnly=TRUE)
sample_id <- as.character(args[1])
var_class <- as.character(args[2])
exp_dat <- as.character(args[3])
iterator <- as.character(args[4])

iterator_file <- fread(iterator, sep=",", header=T)

## Load corrected expression data z-scores
dat_zscore <- fread(exp_dat, sep="\t", header=T)
dat_zscore_df <- data.frame(dat_zscore[, -1])
rownames(dat_zscore_df) <- as.character(unlist(dat_zscore[, 1]))

## Load and process variant data for sample
ids <- fread("[list_of_sample_names].txt", header=TRUE) # load sample data
variant_file <- paste0("IPSCORE_HipSci_1kG_gnomadAF_CADD_", var_class, ".", sample_id, ".bed") # variant file prefix and suffix

## Read and process variant data for current sample
var_list <- list()
tmp <- fread(paste0("location/of/processed/variant/files/", variant_file), header=FALSE, sep="\t")
colnames(tmp) <- c("gene_id", "chr", "pos", "maf", "cadd_phred", "cadd_raw")
tmp <- tmp[-grep(",", tmp$maf), ] # remove multi-allelic calls
tmp$maf[which(tmp$maf==".")] <- -1 # convert singleton to numeric
tmp$maf <- as.numeric(tmp$maf) # convert maf column to numeric
tmp$cadd_phred[which(tmp$cadd_phred==".")] <- -1 # convert singleton to numeric
tmp$cadd_phred <- as.numeric(tmp$cadd_phred) # convert maf column to numeric

rna_id <- ids$RNA[grep(strsplit(variant_file, '[.]')[[1]][2], ids$WGS_ID)] 
var_list[[rna_id]] <- tmp

## Subset to samples with WGS
dat_filter_zscore <- dat_zscore_df[grep(names(var_list), as.character(unlist(rownames(dat_zscore_df)))), ]

message(nrow(dat_filter_zscore))

if (nrow(dat_filter_zscore)>0) {
	for (n in 1:nrow(iterator_file)) {
		z_thresh <- iterator_file$zscore[n]
		maf_min <- iterator_file$min_maf[n]
		maf_max <- iterator_file$max_maf[n]
		cadd_thresh <- iterator_file$cadd[n]

		## Find expression outliers and match with MAF to compute ratio of proportions of
		## outliers with variants and non-outliers with variants
		if (z_thresh<0) {
			outliers <- colnames(dat_filter_zscore)[which(dat_filter_zscore<z_thresh)]
        		extreme_out <- apply(dat_zscore_df, 2, function(x) names(which.min(x))) # which sample is most-extreme under-exp outlier
			extreme_out_sample <- names(extreme_out[grep(rna_id, extreme_out)])
        		outliers_final <- outliers[which(outliers %in% extreme_out_sample)]
		} else {
			outliers <- colnames(dat_filter_zscore)[which(dat_filter_zscore>z_thresh)]
                	extreme_out <- apply(dat_zscore_df, 2, function(x) names(which.max(x))) # which sample is most-extreme over-exp outlier
                	extreme_out_sample <- names(extreme_out[grep(rna_id, extreme_out)])
                	outliers_final <- outliers[which(outliers %in% extreme_out_sample)]	
		}
		non_outliers <- colnames(dat_filter_zscore)[which(dat_filter_zscore > -1 & dat_filter_zscore < 1)]

		## Subset variants to current MAF threshold
		var_list_subset <- unique(as.character(unlist(lapply(var_list, function(x) subset(x, maf>maf_min & maf<=maf_max & cadd_phred>=cadd_thresh)$gene_id))))

		## Collect results
		collect <- data.frame(names(dat_filter_zscore), stringsAsFactors=F)
		colnames(collect)[1] <- "gene_id"
		collect[, c("out_var", "out_all", "non_out_var", "non_out_all")] <- 0

		outliers_var <- outliers_final[which(outliers_final %in% var_list_subset)]
		collect$out_var[which(collect$gene_id %in% outliers_var)] <- 1
		collect$out_all[which(collect$gene_id %in% outliers_final)] <- 1 

		non_outliers_var <- non_outliers[which(non_outliers %in% var_list_subset)]
		collect$non_out_var[which(collect$gene_id %in% non_outliers_var)] <- 1
		collect$non_out_all[which(collect$gene_id %in% non_outliers)] <- 1

		## Write data
		write.table(collect, file=paste0("/out_var_gnomad_", var_class, "_maf-", round(maf_max, 6), "_z", z_thresh, "_cadd", cadd_thresh, "_", sample_id, ".txt"), col.names=T, row.names=F, sep="\t")
	}
}
