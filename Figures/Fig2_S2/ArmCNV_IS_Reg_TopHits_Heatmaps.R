######################################################## HEATMAP TOP ARM HITS ########################################################
source("R/get_top_regression_hits.R")

reg_out_con <- readRDS("../IS_CNV_Regression/TCGA/Arm/Regression_Outputs_is0.3_ArmCNV_continuous0.2.rds") ## continuous results
remove = c("BRCA","CHOL","COAD", "COAD.MSI", "COAD.MSS", "COADREAD", "DLBC", "ESCA", "HNSC", 
           "HNSC.TCGAdef.HPVneg", "HNSC.TCGAdef.HPVpos", "HNSC.HPVother", "KIPAN", "READ", "READ.MSI", "READ.MSS", "STES", "UCEC")

# continuous results
hm1 <- get_top_reg(reg_out = reg_out_con, remove = NULL, cnvlevel = "Arm", cluster_rows = F, cluster_columns = F, continuous = T,
                   show_column_names = F, show_row_names = T, savepath = "IS_CNV_Regression/TCGA/Arm/CON_alltumortypes_top_Cytobands")

hm2 <- get_top_reg(reg_out = reg_out_con, remove = remove, cnvlevel = "Arm", cluster_rows = F, cluster_columns = F, continuous = T,
                   show_column_names = F, show_row_names = T, savepath = "IS_CNV_Regression/TCGA/Arm/CON_filttumors_top_Cytobands")
write.csv(hm2$allregoutdf, "IS_CNV_Regression/TCGA/Arm/continuous_filttumors_tophits.csv", row.names = F)

hm3 <- get_top_reg(reg_out = reg_out_con, remove = remove, cnvlevel = "Arm", cluster_rows = T, cluster_columns = F, continuous = T,
                   show_column_names = F, show_row_names = T, savepath = "IS_CNV_Regression/TCGA/Arm/CON_filttumors_top_Cytobands_clust")