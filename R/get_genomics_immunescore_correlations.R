library(tidyverse)

### function to run dna rna and immunescore correlations
run_corr <- function(dnadf, rnadf, isdf, cnv_thresh=0.2) {
  
  savecorr <- data.frame()
  
  for (i in 1:ncol(dnadf)) { # iterate through all genes, get corr
    
    if (i == ceiling(ncol(dnadf)/2) ) { print("half genes done") }
    
    df <- cbind( dnadf[,i, drop=F], rnadf[,i, drop=F], isdf[,"Ranked_Sum",drop=F]) # get all omics in one df
    colnames(df)[1:2] <- c(paste0(colnames(df)[1], "_DNA"), paste0(colnames(df)[2], "_RNA"))
    
    delfreq <- nrow(df[df[,1] < -cnv_thresh,]) / nrow(df) # deletion frequency of gene
    ampfreq <- nrow(df[df[,1] > cnv_thresh,]) / nrow(df) # amplification frequency of gene
    
    ### all samples correlations
    dna_rna_corr <- cor.test(df[,1], df[,2],  method = "spearman", exact = F) # DNA vs RNA
    dna_is_corr <- cor.test(df[,1], df[,3],  method = "spearman", exact = F) # DNA vs IS
    rna_is_corr <- cor.test(df[,2], df[,3],  method = "spearman", exact = F) # RNA vs IS
    
    ### deletions only
    df_del <- df[df[,1] < cnv_thresh, ] # less than 0.2 (neutral or deletions)
    dna_rna_del_corr <- cor.test(df_del[,1], df_del[,2],  method = "spearman", exact = F) # DNA vs RNA
    dna_is_del_corr <- cor.test(df_del[,1], df_del[,3],  method = "spearman", exact = F) # DNA vs IS
    rna_is_del_corr <- cor.test(df_del[,2], df_del[,3],  method = "spearman", exact = F) # RNA vs IS
    
    ### amplifications only
    df_amp <- df[df[,1] > -cnv_thresh, ] # greater than -0.2 (neutral or amplifications)
    dna_rna_amp_corr <- cor.test(df_amp[,1], df_amp[,2],  method = "spearman", exact = F) # DNA vs RNA
    dna_is_amp_corr <- cor.test(df_amp[,1], df_amp[,3],  method = "spearman", exact = F) # DNA vs IS
    rna_is_amp_corr <- cor.test(df_amp[,2], df_amp[,3],  method = "spearman", exact = F) # RNA vs IS
    
    ### save
    savecorr <- rbind(savecorr, data.frame(Gene = colnames(dnadf)[i],
                                           freq_del_0.2 = delfreq,
                                           freq_amp_0.2 = ampfreq,
                                           
                                           corr_DNA_RNA = unname(dna_rna_corr$estimate), 
                                           corr_DNA_RNA_Pval = unname(dna_rna_corr$p.value),
                                           corr_DNA_IS = unname(dna_is_corr$estimate), 
                                           corr_DNA_IS_Pval = unname(dna_is_corr$p.value),
                                           corr_RNA_IS = unname(rna_is_corr$estimate), 
                                           corr_RNA_IS_Pval = unname(rna_is_corr$p.value),
                                           
                                           corr_DNA_RNA_del = unname(dna_rna_del_corr$estimate), 
                                           corr_DNA_RNA_del_Pval = unname(dna_rna_del_corr$p.value),
                                           corr_DNA_IS_del = unname(dna_is_del_corr$estimate), 
                                           corr_DNA_IS_del_Pval = unname(dna_is_del_corr$p.value),
                                           corr_RNA_IS_del = unname(rna_is_del_corr$estimate), 
                                           corr_RNA_IS_del_Pval = unname(rna_is_del_corr$p.value),
                                           
                                           corr_DNA_RNA_amp = unname(dna_rna_amp_corr$estimate), 
                                           corr_DNA_RNA_amp_Pval = unname(dna_rna_amp_corr$p.value),
                                           corr_DNA_IS_amp = unname(dna_is_amp_corr$estimate), 
                                           corr_DNA_IS_amp_Pval = unname(dna_is_amp_corr$p.value),
                                           corr_RNA_IS_amp = unname(rna_is_amp_corr$estimate), 
                                           corr_RNA_IS_amp_Pval = unname(rna_is_amp_corr$p.value) ))
  }
  return(savecorr)
}
