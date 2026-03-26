source("R/immunescore_cnv_regression.R") # function to run regression models
source("R/plot_area_immunescore_cnv_regression.R") # function to plot outputs


################################ Read in data ################################
working_tables <- readRDS(file = "data/TCGA_working_tables.rds")
cancers <- unique(gsub("_.*", "", names(working_tables)))

xcelltables <- readRDS("data/TCGA_xCell_Tables.RDS")
cancers <- cancers[!cancers %in% setdiff(cancers, names(xcelltables))] #remove KIPAN, LAML, and STES from downstream analysis
# all(cancers %in% names(xcelltables))


outpath = "Deconvolutions_CNV_Regressions/TCGA/outputs/"
xcell_outputs_con <- list()

for (ct in cancers) {
  print(ct)
  
  cnv_df = working_tables[[paste0(ct, "_ArmCNV")]] # cnv table
  cnv_df = cbind(AS = rowSums(abs(cnv_df)), cnv_df) # add Aneuploidy Score (AS)

  xcelltable = xcelltables[[ct]] # table with all xcell scores

  # empty list per cancer type
  xcell_outputs_con[[ct]] <- list()

  ### iterate thru all xcell scores
  for (i in 1:ncol(xcelltable)) {
    if ( sum(xcelltable[,i])==0 ) { next } # if no values in column, move on to next iter

    xcellname=names(xcelltable)[i] # name of score
    print(xcellname)

    ### cnv continuous
    xcell_outputs_con[[ct]][[xcellname]] <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = xcelltable[,i,drop=F], 
                                                                       use_AS = T, cnvlevel = "Arm", 
                                                                       cnv_thresh = 0.2, is_thresh = 0.3, 
                                                                       CNV_binary = F, IS_binary = T, model = "log")
  }
  
  print(paste0("done xcell for ", ct) )
}

saveRDS(xcell_outputs_con, file = paste0(outpath, "Regression_Outputs_xcell0.3_ARMcnv_continuous0.2.rds"))





############################################################ PLOTTING ############################################################

############################## plot xcell continuous outputs ##############################
# xcell_outputs_con <- readRDS(paste0(outpath, "Regression_Outputs_xcell0.3_ARMcnv_continuous0.2.rds"))
conplots_pval <- list() # empty lists for holding all cancer types
for (ct in names(xcell_outputs_con)) {
  print(ct)
  
  conplots_pval[[ct]] <- list() # empty lists to hold all score plots per cancer type
  for (score in names(xcell_outputs_con[[ct]])) { # iterate thru all xcell scores per cancer type
    print(paste0("xCell ", score))
    
    # continuous cnv, p val
    conplots_pval[[ct]][[score]] <- plot_cnv_IS_regression(df = xcell_outputs_con[[ct]][[score]], 
                                                           cnvlevel = "Arm", scorename = score, stat = "P_val",
                                                           cnvtype = "con", signif = 2, plottitle = score)
  }
  
  ### save multiple plots per cancer type to one pdf
  fig <- ggarrange(plotlist = conplots_pval[[ct]], common.legend = T, legend="bottom")
  pdf(paste0("Deconvolutions_CNV_Regressions/TCGA/xCell_Continuous_Arm_plots/", ct, "_xCell0.3_ARMcnv_continuous0.2_reg_plot_logP.pdf"), width = 24, height = 12)
  print(annotate_figure(fig, top = text_grob("xCell(con0.3) ~ Arm_CNV(bi0.2) + AS", color = "black", size = 12)))
  junk <- dev.off()
}

print("Done xcell continuous plots")











