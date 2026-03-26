source("R/immunescore_cnv_regression.R") # function to run regression models
source("R/plot_area_immunescore_cnv_regression.R") # function to plot outputs


############################### save xcell and cibersort tables ################################
# working_tables <- readRDS(file = "data/TCGA_working_tables.rds")
# subtypes <- unique(gsub("_.*", "", names(working_tables)))
# 
# tcgasubtypes <- read.csv("../TCGA_Cox/ABSOLUTE.purity.txt", sep="\t") # table with cancer names and subtypes - all samples in each group
# tcgasubtypes <- tcgasubtypes[!tcgasubtypes$cancer %in% "GBMLGG",] # remove gbm lgg
# 
# # xcell files
# xcellpath <- "../TCGA_Deconvolution/xcell_TCGA/"
# xcellfiles <- list.files(xcellpath)
# cancers <- gsub("_xcell.tsv", "", xcellfiles)
# 
# xcelltables <- list()
# for (ct in cancers) {
#    xcelldf <- read.csv(paste0(xcellpath, ct, "_xcell.tsv"), sep = "\t", check.names = F)
#    colnames(xcelldf) <- gsub(" ", "_", colnames(xcelldf))
#    colnames(xcelldf) <- gsub("\\(|\\)|\\+|-", "", colnames(xcelldf))
#    colnames(xcelldf)[1] <- "Sample"
#    xcelldf$Sample <- substr(xcelldf$Sample, 1, nchar(xcelldf$Sample)-3)
#    xcelldf <- xcelldf %>% group_by(Sample) %>% summarize(across(everything(), mean)) # duplicated samples - take mean across duplicates
#    xcelldf <- column_to_rownames(xcelldf, var = "Sample")
# 
#    xcelltables[[ct]] <- xcelldf
# 
#    possible_subtypes <- unique(grep(paste0(ct, "."), subtypes, value = T, fixed = T)) # get subtypes
#    if (ct=="READ") { possible_subtypes <- possible_subtypes[!grepl("COAD", possible_subtypes)] }
#    if ( length(possible_subtypes) > 1 ) {
#     for (sub in possible_subtypes) {
#       xcelltables[[sub]] <- xcelltables[[ct]][rownames(xcelltables[[ct]]) %in% tcgasubtypes[tcgasubtypes$cancer == sub, "sample"],]
#     }
#   }
# }
# # setdiff(subtypes, names(xcelltables)) ### add coadread (combine coad and read), remove KIPAN, LAML, and STES from downstream analysis
# xcelltables[["COADREAD"]] <- rbind(xcelltables$COAD, xcelltables$READ)
# xcelltables[["COADREAD.MSI"]] <- rbind(xcelltables$COAD.MSI, xcelltables$READ.MSI)
# xcelltables[["COADREAD.MSS"]] <- rbind(xcelltables$COAD.MSS, xcelltables$READ.MSS)
# xcelltables <- xcelltables[sort(names(xcelltables))]
# saveRDS(xcelltables, "data/TCGA_xCell_Tables.RDS")
# 
# 
# 
# # Cibersort files
# cibersortpath <- "../TCGA_Deconvolution/CIBERSORT_TCGA/"
# cibersortfiles <- list.files(cibersortpath)
# cancers <- gsub("_CIBERSORT.tsv", "", cibersortfiles)
# 
# cibersorttables <- list()
# for (ct in cancers) {
#   cibersortdf <- read.csv(paste0(cibersortpath, ct, "_CIBERSORT.tsv"), sep = "\t", check.names = F)
#   cibersortdf$`P-value` <- NULL; cibersortdf$Correlation <- NULL; cibersortdf$RMSE <- NULL
#   colnames(cibersortdf) <- gsub(" ", "_", colnames(cibersortdf))
#   colnames(cibersortdf) <- gsub("\\(|\\)", "", colnames(cibersortdf))
#   colnames(cibersortdf)[1] <- "Sample"
#   cibersortdf$Sample <- substr(cibersortdf$Sample, 1, nchar(cibersortdf$Sample)-3)
#   cibersortdf <- cibersortdf %>% group_by(Sample) %>% summarize(across(everything(), mean))
#   cibersortdf <- column_to_rownames(cibersortdf, var = "Sample")
# 
#   cibersorttables[[ct]] <- cibersortdf
# 
#   possible_subtypes <- unique(grep(paste0(ct, "."), subtypes, value = T, fixed = T)) # get subtypes
#   if (ct=="READ") { possible_subtypes <- possible_subtypes[!grepl("COAD", possible_subtypes)] }
#   if ( length(possible_subtypes) > 1 ) {
#     for (sub in possible_subtypes) {
#       cibersorttables[[sub]] <- cibersorttables[[ct]][rownames(cibersorttables[[ct]]) %in% tcgasubtypes[tcgasubtypes$cancer == sub, "sample"],]
#     }
#   }
# }
# 
# # setdiff(subtypes, names(cibersorttables)) ### add coadread (combine coad and read), remove KIPAN, LAML, and STES from downstream analysis
# cibersorttables[["COADREAD"]] <- rbind(cibersorttables$COAD, cibersorttables$READ)
# cibersorttables[["COADREAD.MSI"]] <- rbind(cibersorttables$COAD.MSI, cibersorttables$READ.MSI)
# cibersorttables[["COADREAD.MSS"]] <- rbind(cibersorttables$COAD.MSS, cibersorttables$READ.MSS)
# cibersorttables <- cibersorttables[sort(names(cibersorttables))]
# saveRDS(cibersorttables, "data/TCGA_cibersort_Tables.RDS")



################################ Read in data ################################
working_tables <- readRDS(file = "data/TCGA_working_tables.rds")
cancers <- unique(gsub("_.*", "", names(working_tables)))

xcelltables <- readRDS("data/TCGA_xCell_Tables.RDS")
# cibersorttables <- readRDS("data/TCGA_cibersort_Tables.RDS")
# 
cancers <- cancers[!cancers %in% setdiff(cancers, names(xcelltables))] #remove KIPAN, LAML, and STES from downstream analysis
# all(cancers %in% names(xcelltables))
# all(cancers %in% names(cibersorttables))
# 
# ################################ Run Regressions ################################
# outpath = "Deconvolutions_CNV_Regressions/TCGA/"
# 
# xcell_outputs_bi <- xcell_outputs_con <- cibersort_outputs_bi <- cibersort_outputs_con <- list() # empty lists
# 
# for (ct in cancers) {
#   print(ct)
# 
#   cnv_df = working_tables[[paste0(ct, "_CytoCNV")]] # cnv table
#   cnv_df = cbind(AS = rowSums(abs(cnv_df)), cnv_df) # add Aneuploidy Score (AS)
# 
#   xcelltable = xcelltables[[ct]] # table with all xcell scores
#   cibersorttable = cibersorttables[[ct]] # table with all cibersort scores
# 
#   # empty lists per cancer type
#   xcell_outputs_bi[[ct]] <- xcell_outputs_con[[ct]] <- cibersort_outputs_bi[[ct]] <- cibersort_outputs_con[[ct]] <- list()
# 
#   ### iterate thru all xcell scores
#   for (i in 1:ncol(xcelltable)) {
#     if ( sum(xcelltable[,i])==0 ) { next } # if no values in column, move on to next iter
# 
#     xcellname=names(xcelltable)[i] # name of score
# 
#     ### cnv binary
#     xcell_outputs_bi[[(ct)]][[xcellname]] <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = xcelltable[,i,drop=F],
#                                                                         use_AS = T, cnvlevel = "Cytoband",
#                                                                         cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = T, IS_binary = T,
#                                                                         model = "log", savepath = paste0(outpath, "xCell_outputs/", ct, "_", xcellname, "_xcell0.3_cytocnv_binary0.2_regression.csv"))
# 
#     ### cnv continuous
#     xcell_outputs_con[[ct]][[xcellname]] <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = xcelltable[,i,drop=F],
#                                                                        use_AS = T, cnvlevel = "Cytoband",
#                                                                        cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = F, IS_binary = T,
#                                                                        model = "log", savepath = paste0(outpath, "xCell_outputs/", ct, "_", xcellname, "_xcell0.3_cytocnv_continuous0.2_regression.csv"))
#   }
#   print("done xcell")
# 
#   ### iterate thru all cibersort scores
#   for (i in 1:ncol(cibersorttable)) {
#     if ( sum(cibersorttable[,i])==0 ) { next } # if no values in column, move on to next iter
# 
#     cibname=names(cibersorttable)[i] # name of score
# 
#     ### cnv binary
#     cibersort_outputs_bi[[(ct)]][[cibname]] <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = cibersorttable[,i,drop=F],
#                                                                           use_AS = T, cnvlevel = "Cytoband",
#                                                                           cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = T, IS_binary = T,
#                                                                           model = "log", savepath = paste0(outpath, "CIBERSORT_outputs/", ct, "_", cibname, "_cibersort0.3_cytocnv_binary0.2_regression.csv"))
# 
#     ### cnv continuous
#     cibersort_outputs_con[[ct]][[cibname]] <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = cibersorttable[,i,drop=F],
#                                                                          use_AS = T, cnvlevel = "Cytoband",
#                                                                          cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = F, IS_binary = T,
#                                                                          model = "log", savepath = paste0(outpath, "CIBERSORT_outputs/", ct, "_", cibname, "_cibersort0.3_cytocnv_continuous0.2_regression.csv"))
#   }
#   print("done cibersort")
# }
# 
# ### save object
# saveRDS(xcell_outputs_bi, file = paste0(outpath, "Regression_Outputs_xcell0.3_cytobandcnv_binary0.2.rds"))
# saveRDS(xcell_outputs_con, file = paste0(outpath, "Regression_Outputs_xcell0.3_cytobandcnv_continuous0.2.rds"))
# 
# saveRDS(cibersort_outputs_bi, file = paste0(outpath, "Regression_Outputs_cibersort0.3_cytobandcnv_binary0.2.rds"))
# saveRDS(cibersort_outputs_con, file = paste0(outpath, "Regression_Outputs_cibersort0.3_cytobandcnv_continuous0.2.rds"))
# 
# print("saved rds files")


############################################################ PLOTTING ############################################################


############################## plot xcell binary outputs ##############################
xcell_outputs_bi <- readRDS("Deconvolutions_CNV_Regressions/TCGA/Regression_Outputs_xcell0.3_cytobandcnv_binary0.2.rds")
biplots_zval <- biplots_pval <- list() # empty lists for holding all cancer types
for (ct in names(xcell_outputs_bi)) {
  print(ct)

  biplots_zval[[ct]] <- biplots_pval[[ct]] <- list() # empty lists to hold all score plots per cancer type
  for (score in names(xcell_outputs_bi[[ct]])) { # iterate thru all xcell scores per cancer type
    print(paste0("xCell ", score))

    # binary cnv, z val
    p1 <- plot_cnv_IS_regression(df = xcell_outputs_bi[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "Z_val",
                                 cnvtype = "bi", signif = 2, plottitle = score)
    # binary cnv, p val
    p2 <- plot_cnv_IS_regression(df = xcell_outputs_bi[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "P_val",
                                 cnvtype = "bi", signif = 2, plottitle = score)
    biplots_zval[[ct]][[score]] <- p1
    biplots_pval[[ct]][[score]] <- p2
  }

  ### save multiple plots per cancer type to one pdf
  fig1 <- ggarrange(plotlist = biplots_zval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/xCell_plots/", ct, "_xCell0.3_cytobandcnv_binary0.2_reg_plot_Zval.pdf"), width = 24, height = 12)
  print(annotate_figure(fig1, top = text_grob("xCell(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()

  fig2 <- ggarrange(plotlist = biplots_pval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/xCell_plots/", ct, "_xCell0.3_cytobandcnv_binary0.2_reg_plot_logP.pdf"), width = 24, height = 12)
  print(annotate_figure(fig2, top = text_grob("xCell(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
}

print("Done xcell binary plots")


############################## plot xcell continuous outputs ##############################
xcell_outputs_con <- readRDS("Deconvolutions_CNV_Regressions/TCGA/Regression_Outputs_xcell0.3_cytobandcnv_continuous0.2.rds")
conplots_zval <- conplots_pval <- list() # empty lists for holding all cancer types
for (ct in names(xcell_outputs_con)) {
  print(ct)

  conplots_zval[[ct]] <- conplots_pval[[ct]] <- list() # empty lists to hold all score plots per cancer type
  for (score in names(xcell_outputs_con[[ct]])) { # iterate thru all xcell scores per cancer type
    print(paste0("xCell ", score))

    # continuous cnv, z val
    p1 <- plot_cnv_IS_regression(df = xcell_outputs_con[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "Z_val",
                                 cnvtype = "con", signif = 2, plottitle = score)
    # continuous cnv, p val
    p2 <- plot_cnv_IS_regression(df = xcell_outputs_con[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "P_val",
                                 cnvtype = "con", signif = 2, plottitle = score)
    conplots_zval[[ct]][[score]] <- p1
    conplots_pval[[ct]][[score]] <- p2
  }

  ### save multiple plots per cancer type to one pdf
  fig1 <- ggarrange(plotlist = conplots_zval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/xCell_plots/", ct, "_xCell0.3_cytobandcnv_continuous0.2_reg_plot_Zval.pdf"), width = 24, height = 12)
  print(annotate_figure(fig1, top = text_grob("xCell(con0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()

  fig2 <- ggarrange(plotlist = conplots_pval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/xCell_plots/", ct, "_xCell0.3_cytobandcnv_continuous0.2_reg_plot_logP.pdf"), width = 24, height = 12)
  print(annotate_figure(fig2, top = text_grob("xCell(con0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
}

print("Done xcell continuous plots")


############################## plot cibersort binary outputs ##############################
cibersort_outputs_bi <- readRDS("Deconvolutions_CNV_Regressions/TCGA/Regression_Outputs_cibersort0.3_cytobandcnv_binary0.2.rds")
biplots_zval <- biplots_pval <- list() # empty lists for holding all cancer types
for (ct in names(cibersort_outputs_bi)) {
  print(ct)
  
  biplots_zval[[ct]] <- biplots_pval[[ct]] <- list() # empty lists to hold all score plots per cancer type
  for (score in names(cibersort_outputs_bi[[ct]])) { # iterate thru all cibersort scores per cancer type
    print(paste0("cibersort ", score))
    
    # binary cnv, z val
    p1 <- plot_cnv_IS_regression(df = cibersort_outputs_bi[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "Z_val",
                                 cnvtype = "bi", signif = 2, plottitle = score)
    # binary cnv, p val
    p2 <- plot_cnv_IS_regression(df = cibersort_outputs_bi[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "P_val",
                                 cnvtype = "bi", signif = 2, plottitle = score)
    biplots_zval[[ct]][[score]] <- p1
    biplots_pval[[ct]][[score]] <- p2
  }
  
  ### save multiple plots per cancer type to one pdf
  fig1 <- ggarrange(plotlist = biplots_zval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/CIBERSORT_plots/", ct, "_cibersort0.3_cytobandcnv_binary0.2_reg_plot_Zval.pdf"), width = 24, height = 12)
  print(annotate_figure(fig1, top = text_grob("CIBERSORT(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
  
  fig2 <- ggarrange(plotlist = biplots_pval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/CIBERSORT_plots/", ct, "_cibersort0.3_cytobandcnv_binary0.2_reg_plot_logP.pdf"), width = 24, height = 12)
  print(annotate_figure(fig2, top = text_grob("CIBERSORT(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
}

print("Done cibersort binary plots")

### plot cibersort continuous outputs
cibersort_outputs_con <- readRDS("Deconvolutions_CNV_Regressions/TCGA/Regression_Outputs_cibersort0.3_cytobandcnv_continuous0.2.rds")
conplots_zval <- conplots_pval <- list() # empty lists for holding all cancer types
for (ct in names(cibersort_outputs_con)) {
  print(ct)
  
  conplots_zval[[ct]] <- conplots_pval[[ct]] <- list() # empty lists to hold all score plots per cancer type
  for (score in names(cibersort_outputs_con[[ct]])) { # iterate thru all cibersort scores per cancer type
    print(paste0("cibersort ", score))
    
    # continuous cnv, z val
    p1 <- plot_cnv_IS_regression(df = cibersort_outputs_con[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "Z_val",
                                 cnvtype = "con", signif = 2, plottitle = score)
    # continuous cnv, p val
    p2 <- plot_cnv_IS_regression(df = cibersort_outputs_con[[ct]][[score]], cnvlevel = "Cytoband", scorename = score, stat = "P_val",
                                 cnvtype = "con", signif = 2, plottitle = score)
    conplots_zval[[ct]][[score]] <- p1
    conplots_pval[[ct]][[score]] <- p2
  }
  
  ### save multiple plots per cancer type to one pdf
  fig1 <- ggarrange(plotlist = conplots_zval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/CIBERSORT_plots/", ct, "_cibersort0.3_cytobandcnv_continuous0.2_reg_plot_Zval.pdf"), width = 24, height = 12)
  print(annotate_figure(fig1, top = text_grob("CIBERSORT(bi0.3) ~ Cytoband_CNV(con0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
  
  fig2 <- ggarrange(plotlist = conplots_pval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/CIBERSORT_plots/", ct, "_cibersort0.3_cytobandcnv_continuous0.2_reg_plot_logP.pdf"), width = 24, height = 12)
  print(annotate_figure(fig2, top = text_grob("CIBERSORT(bi0.3) ~ Cytoband_CNV(con0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
}
print("Done cibersort continuous plots")







