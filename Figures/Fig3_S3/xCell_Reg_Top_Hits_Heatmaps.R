######################################################## HEATMAP TOP HITS - XCELL ########################################################
source("R/get_top_regression_hits.R")


xcellreg_cyto <- readRDS("Deconvolutions_CNV_Regressions/TCGA/outputs/Regression_Outputs_xcell0.3_cytobandcnv_continuous0.2.rds")
xcellreg_arm <- readRDS("Deconvolutions_CNV_Regressions/TCGA/outputs/Regression_Outputs_xcell0.3_ARMcnv_continuous0.2.rds")
# workingtables <- readRDS(file = "data/TCGA_working_tables.rds")
# remove = c("BRCA","CHOL","COAD", "COAD.MSI", "COAD.MSS", "COADREAD", "DLBC", "ESCA", "HNSC", 
#            "HNSC.TCGAdef.HPVneg", "HNSC.TCGAdef.HPVpos", "HNSC.HPVother", "KIPAN", "READ", "READ.MSI", "READ.MSS", "STES", "UCEC")
# tumors <- unique(gsub("_.*", "", names(workingtables)))
# tumors <- tumors[!tumors %in% remove]


tumor_order = c("UCEC.MSI", "UCEC.MSS", "OV", "CESC", "LUSC", "HNSC.HPVneg", "ESCA.SC","BLCA",
                "BRCA.pos", "BRCA.neg", "PRAD", "LUAD", "LIHC", "COADREAD.MSI", "COADREAD.MSS", "ESCA.AD", "STAD",
                "PAAD", "KICH", "KIRC", "KIRP", "LGG", "GBM", "SKCM", "UVM", "ACC", "MESO", "PCPG", "SARC", "TGCT", "UCS")
tumor_groups = data.frame(row.names=tumor_order, 
                          Group=c("Gyn", "Gyn", "Gyn", "Squamous", "Squamous", "Squamous", "Squamous", "Squamous", 
                                  "Adeno", "Adeno", "Adeno", "Adeno", "Adeno", "GI", "GI", "GI", "GI", "GI",
                                  "Kidney", "Kidney", "Kidney", "NC-Derived", "NC-Derived", "NC-Derived", "NC-Derived", 
                                  "Other", "Other","Other","Other","Other","Other"))
tumor_group_cols <- c("Gyn" = "#6C2DC7", "Squamous" = "#028A0F", "Adeno" = "#6E0B14", 
                      "GI" = "#FFC30B", "Kidney" = "#2B65EC", "NC-Derived" = "#9B111E", "Other" = "#008080")

# cancers = unique(gsub("_.*", "", names(workingtables))) # all tumor types
cancers=names(xcellreg_cyto)
remove = setdiff(cancers, tumor_order) # tumors to remove from analysis

tumors = tumor_order

tcells <- c("CD8_T_cells", "CD4_T_cells", "Tregs")
myeloid <- c("Dendritic_cells", "Macrophage_M1", "Macrophage_M2")
other <- c("NK_cells", "Neutrophils", "B_cells")
nonimmune <- c("Endothelium", "Fibroblasts")


tcell_outputs_cyto <- list()
tcell_outputs_arm <- list()
for (tt in tumors) {
  for (celltype in tcells) {
    tcell_outputs_cyto[[celltype]][[tt]] <- xcellreg_cyto[[tt]][[celltype]]
    tcell_outputs_arm[[celltype]][[tt]] <- xcellreg_arm[[tt]][[celltype]]
  }
}

myeloid_outputs_cyto <- list()
myeloid_outputs_arm <- list()
for (tt in tumors) {
  for (celltype in myeloid) {
    myeloid_outputs_cyto[[celltype]][[tt]] <- xcellreg_cyto[[tt]][[celltype]]
    myeloid_outputs_arm[[celltype]][[tt]] <- xcellreg_arm[[tt]][[celltype]]
  }
}

other_outputs_cyto <- list()
other_outputs_arm <- list()
for (tt in tumors) {
  for (celltype in other) {
    other_outputs_cyto[[celltype]][[tt]] <- xcellreg_cyto[[tt]][[celltype]]
    other_outputs_arm[[celltype]][[tt]] <- xcellreg_arm[[tt]][[celltype]]
  }
}

nonimmune_outputs_cyto <- list()
nonimmune_outputs_arm <- list()
for (tt in tumors) {
  for (celltype in nonimmune) {
    nonimmune_outputs_cyto[[celltype]][[tt]] <- xcellreg_cyto[[tt]][[celltype]]
    nonimmune_outputs_arm[[celltype]][[tt]] <- xcellreg_arm[[tt]][[celltype]]
  }
}


###### tcells
## cd8
# arm
get_top_reg(reg_out = tcell_outputs_arm$CD8_T_cells, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/CD8_tcells_")
# cyto
get_top_reg(reg_out = tcell_outputs_cyto$CD8_T_cells, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/CD8_tcells_")

## cd4
# arm
get_top_reg(reg_out = tcell_outputs_arm$CD4_T_cells, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols,
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/CD4_tcells_")
# cyto
get_top_reg(reg_out = tcell_outputs_cyto$CD4_T_cells, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/CD4_tcells_")

## Tregs
# arm
get_top_reg(reg_out = tcell_outputs_arm$Tregs, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols,
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Tregs_")
# cyto
get_top_reg(reg_out = tcell_outputs_cyto$Tregs, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Tregs_")



###### myeloid cells
## Dendritic cells
# arm
get_top_reg(reg_out = myeloid_outputs_arm$Dendritic_cells, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Dendritic_cells_")
# cyto
get_top_reg(reg_out = myeloid_outputs_cyto$Dendritic_cells, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Dendritic_cells_")

## Macrophage_M1
# arm
get_top_reg(reg_out = myeloid_outputs_arm$Macrophage_M1, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Macrophage_M1_")
# cyto
get_top_reg(reg_out = myeloid_outputs_cyto$Macrophage_M1, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Macrophage_M1_")

## Macrophage_M2
# arm
get_top_reg(reg_out = myeloid_outputs_arm$Macrophage_M2, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Macrophage_M2_")
# cyto
get_top_reg(reg_out = myeloid_outputs_cyto$Macrophage_M2, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Macrophage_M2_")





###### other cells
## NK_cells
# arm
get_top_reg(reg_out = other_outputs_arm$NK_cells, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/NK_cells_")
# cyto
get_top_reg(reg_out = other_outputs_cyto$NK_cells, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/NK_cells_")

## Neutrophils
# arm
get_top_reg(reg_out = other_outputs_arm$Neutrophils, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Neutrophils_")
# cyto
get_top_reg(reg_out = other_outputs_cyto$Neutrophils, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Neutrophils_")

## B_cells
# arm
get_top_reg(reg_out = other_outputs_arm$B_cells, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/B_cells_")
# cyto
get_top_reg(reg_out = other_outputs_cyto$B_cells, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/B_cells_")


###### non immune cells
## Endothelium
# arm
get_top_reg(reg_out = nonimmune_outputs_arm$Endothelium, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Endothelium_")
# cyto
get_top_reg(reg_out = nonimmune_outputs_cyto$Endothelium, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Endothelium_")

## Fibroblasts
# arm
get_top_reg(reg_out = nonimmune_outputs_arm$Fibroblasts, remove = remove, cnvlevel = "Arm", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/Fibroblasts_")
# cyto
get_top_reg(reg_out = nonimmune_outputs_cyto$Fibroblasts, remove = remove, cnvlevel = "Cytoband", 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, 
            continuous = T, signif_thresh = 0.05, plotvals = "logP", savetable = T,
            savepath = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Cytoband_Hits/Fibroblasts_")





### arms
toparms_files <- list.files(path = "Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/", pattern = "table.csv")
celltypes <- gsub("_table.csv", "", toparms_files)

top_arms_cells_out <- data.frame()
for (celltype in celltypes) {
  df <- read.csv(paste0("Deconvolutions_CNV_Regressions/TCGA/Top_xCell_Continuous_Arm_Hits/", celltype, "_table.csv"))
  df_filt <- cbind(CellType=rep(celltype), na.omit(df[,c("TumorType", "Arm", "AllHits")]))
  top_arms_cells_out <- rbind(top_arms_cells_out, df_filt )
}

arms <- unique(df$Arm)


top_arms_cells_out$Event <- ifelse(grepl("Loss", top_arms_cells_out$AllHits), "Loss", "Gain")
top_arms_cells_out$Hit <- ifelse(grepl("Cold", top_arms_cells_out$AllHits), "Cold", "Hot")
top_arms_cells_out$Arm <- factor(top_arms_cells_out$Arm, levels = arms)


pdf("Deconvolutions_CNV_Regressions/TCGA/Top_Arms_Count_all_celltypes_barplots.pdf", width = 16, height = 10)
ggplot(top_arms_cells_out, aes(x=Arm, fill=Event)) + geom_bar() + facet_grid(rows = vars(CellType), cols = vars(Hit)) + theme_bw()
dev.off()





