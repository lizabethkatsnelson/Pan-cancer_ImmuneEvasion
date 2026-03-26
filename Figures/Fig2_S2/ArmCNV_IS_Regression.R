source("R/immunescore_cnv_regression.R") # function to run regression models
source("R/plot_area_immunescore_cnv_regression.R") # function to plot outputs

working_tables <- readRDS(file = "data/TCGA_working_tables.rds")
cancers <- unique(gsub("_.*", "", names(working_tables)))
# View(working_tables$ACC_ArmCNV)
# View(working_tables$ACC_IS)

outpath = "IS_CNV_Regression/TCGA/Arm/"

all_outputs_bi <- list()
all_outputs_con <- list()

for (ct in cancers) {
  print(ct)
  
  cnv_df = working_tables[[paste0(ct, "_ArmCNV")]]
  cnv_df = cbind(AS = rowSums(abs(cnv_df)), cnv_df) # add Aneuploidy Score (AS)
  is_df = working_tables[[paste0(ct, "_IS")]][,"Ranked_Sum",drop=F]

  ### cnv binary
  reg_out_bi <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = is_df, use_AS = T, cnvlevel = "Arm",
                                           cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = T, IS_binary = T, model = "log", 
                                           savepath = paste0(outpath, "Outputs/", ct, "_is0.3_ArmCNV_binary0.2_regression.csv"))
  all_outputs_bi[[ct]] <- reg_out_bi
  print("done with first regression")
  
  ### cnv continuous
  reg_out_con <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = is_df, use_AS = T, cnvlevel = "Arm",
                                            cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = F, IS_binary = T, model = "log", 
                                            savepath = paste0(outpath, "Outputs/", ct, "_is0.3_ArmCNV_continuous0.2_regression.csv"))
  all_outputs_con[[ct]] <- reg_out_con
  print("done with second regression")
}
### save object
saveRDS(all_outputs_bi, file = "IS_CNV_Regression/TCGA/Arm/Regression_Outputs_is0.3_ArmCNV_binary0.2.rds")
saveRDS(all_outputs_con, file = "IS_CNV_Regression/TCGA/Arm/Regression_Outputs_is0.3_ArmCNV_continuous0.2.rds")
# print("saved rds files")




############################################################ PLOTTING ############################################################
all_outputs_bi <- readRDS("IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_binary0.2.rds")
all_outputs_con <- readRDS("IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_continuous0.2.rds")

### plot binary outputs
for (ct in names(all_outputs_bi)) {
  print(ct)
  
  # # binary cnv, z val
  plot_cnv_IS_regression(df = all_outputs_bi[[ct]], cnvlevel = "Arm", scorename = "IS", stat = "Z_val", cnvtype = "bi",
                         signif = 2, plottitle = "IS(bi0.3) ~ Arm_CNV(bi0.2) + AS",
                         savepath = paste0("IS_CNV_Regression/TCGA/Arm/Plots/", ct, "_is0.3_armcnv_binary0.2_reg_plot_Zval.pdf"))

  # binary cnv, p val
  plot_cnv_IS_regression(df = all_outputs_bi[[ct]], cnvlevel = "Arm", scorename = "IS", stat = "P_val", cnvtype = "bi", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Arm_CNV(bi0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Arm/Plots/", ct, "_is0.3_armcnv_binary0.2_reg_plot_logP.pdf"))
}


### plot continuous outputs
for (ct in names(all_outputs_con)) {
  print(ct)
  
  # continuous cnv, z val
  plot_cnv_IS_regression(df = all_outputs_con[[ct]], cnvlevel = "Arm", scorename = "IS", stat = "Z_val", cnvtype = "con", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Arm_CNV(con0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Arm/Plots/", ct, "_is0.3_armcnv_continuous0.2_reg_plot_Zval.pdf"))
  # continuous cnv, p val
  plot_cnv_IS_regression(df = all_outputs_con[[ct]], cnvlevel = "Arm", scorename = "IS", stat = "P_val", cnvtype = "con", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Arm_CNV(con0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Arm/Plots/", ct, "_is0.3_armcnv_continuous0.2_reg_plot_logP.pdf"))
}