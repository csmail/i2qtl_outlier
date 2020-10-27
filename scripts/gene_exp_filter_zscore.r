#!/bin/R

## Read uncorrected RNA-seq counts, filter genes to protein-coding and LINC only
## Run PEER, remove global expression outlier samples

## Load packages
library(data.table)
library(annotables)
library(ggplot2)
library(peer)

##------------------- FUNCTIONS

## Get all outliers 
pick_global_outliers <- function(dat_zscore, z) {
	outliers <- apply(t(dat_zscore), 1, function(x) list(names(which(abs(x)>z))))
    extreme_out <- apply(t(dat_zscore), 1, function(x) names(which.max(abs(x)))) # which sample is most-extreme outlier
    outliers_extreme_check <- lapply(1:length(extreme_out), function(i) extreme_out[[i]] %in% as.character(unlist(outliers[[i]])))
    outliers_final <- data.frame(table(unlist(extreme_out[which(as.logical(unlist(outliers_extreme_check)))])), stringsAsFactors=F)
    colnames(outliers_final) <- c("sample_id", "freq")
    return(outliers_final)
}

##------------------- MAIN

## Read corrected gene expression data
dat <- fread("[RNAseq_featureCount].txt")
col.names <- as.character(unlist(dat$V1))
dat <- dat[,V1 := NULL]
dat_t <- t(dat)
colnames(dat_t) <- col.names
rownames(dat_t) <- gsub(".hs37d5", "", as.character(unlist(rownames(dat_t)))) # remove hipsci suffix
colnames(dat_t) <- as.character(unlist(sapply(colnames(dat_t), function(x) strsplit(x, '[_]')[[1]][1]))) 

## Subset expression matrix to protein-coding on LINC genes only
gene_list <- fread("[list_of_protein_coding_or_linc_genes].txt", header=F)
dat_filter_gene_list <- dat_t[, which(colnames(dat_t) %in% as.character(unlist(gene_list)))]
dat_filter_gene_list <- dat_filter_gene_list[, -which(colnames(dat_filter_gene_list) %in% subset(grch37, chr %in% c("X","Y"))$ensgene)] # remove X,Y chr

## Subset expression matrix to minimally expressed genes
dat_filter_study <- lapply(unique(fread("[study_key].txt")), function(x) {
		index <- grep(x, fread("[study_key].txt"))
		tmp_dat <- dat_filter_strand[index, ]
		tmp_exp <- names(which(apply(tmp_dat, 2, function(x) sum(x>0)>(nrow(tmp_dat)*0.5))))
	})

dat_filter_genes_intersect <- Reduce(intersect, dat_filter_study)
dat_filter_gene_exp <- dat_filter_strand[, which(colnames(dat_filter_strand) %in% dat_filter_genes_intersect)]

## Check for and remove any genes with zero variance
dat_var <- dat_filter_gene_exp[, apply(dat_filter_gene_exp, 2, function(x) !var(x)==0)]

## PEER factor correction
model = PEER()
PEER_setPhenoMean(model,as.matrix(dat_var))
dim(PEER_getPhenoMean(model))
PEER_setNk(model,200)
PEER_getNk(model)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)
peer_factors <- PEER_getX(model)

## Regress out PEER factors (first 50)
dat_peer_removed <- apply(dat_var, 2, function(x) lm(x~peer_factors[, 50])$residuals)

## Transform to Z-scores
dat_zscore <- scale(dat_peer_removed, center=T, scale=T)

## Remove global expression outliers
global_out_thresh <- 100 # define threshold for global outlier
global_z_thresh <- 2
outlier_count <- pick_global_outliers(dat_zscore, global_z_thresh)
global_out <- as.character(unlist(outlier_count[which(as.numeric(outlier_count[, 2])>global_out_thresh), 1]))
dat_peer_removed_global_remove <- dat_zscore[-which(rownames(dat_zscore) %in% global_out), ]

## Subset to samples with WGS
wgs_id <- fread("[list_of_samples_with_WGS].txt", header=T)
dat_filter_sample <- dat_peer_removed_global_remove[which(rownames(dat_peer_removed_global_remove) %in% as.character(unlist(wgs_id$RNA))), ]
dat_filter_sample <- scale(dat_filter_sample) # re-compute z-scores

## Convert to data table
dat_data_table <- data.table(rownames(dat_filter_sample), dat_filter_sample)
colnames(dat_data_table)[1] <- "sample_id"

## Write data
fwrite(dat_data_table, file="[iPSC_corrected_zscore].txt", sep="\t", col.names=T)

