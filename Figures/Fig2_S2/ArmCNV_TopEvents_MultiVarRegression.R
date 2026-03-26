library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

workingtables <- readRDS("data/TCGA_working_tables.rds")

############### top events per tumor type regressions - using continuous results
reg_con <- read.csv("IS_CNV_Regression/TCGA/Arm/Top_Regression_Hits/Continuous_Top_Arms_table.csv") %>% filter(ColdHits %in% c(-1, 1)) # only top hits

# get regression formulas per tumor type
reg_con_hits_tt <- as.data.frame(reg_con %>% group_by(TumorType) %>% 
                                   summarize(TopHits_Formula = paste(paste0("Arm", Arm, "_", Event), collapse = " + ")) %>%
                                   mutate(NumHits = str_count(TopHits_Formula, "\\+") + 1))

cancers <- reg_con_hits_tt[reg_con_hits_tt$NumHits > 1, "TumorType"] # tumor types with more than one top hit (use for multi regression)


##### RUN MULTI VAR REG

multi_reg_out <- data.frame()
for (tt in cancers) {
  print(tt)
  
  armdf <- workingtables[[paste0(tt, "_ArmCNV")]] # arm cnv df
  isdf <- workingtables[[paste0(tt, "_IS")]] # immune score df
  
  events <- strsplit(reg_con_hits_tt[reg_con_hits_tt$TumorType == tt, "TopHits_Formula"], " \\+ ")[[1]] # get string top events per tt
  eventsdf <- as.data.frame(t(as.data.frame(strsplit(events, "_"))),row.names = NA) # df of top arms and events per tt
  
  ### convert arm df to binary for models
  regdf <- armdf[,eventsdf$V1] # filter arm df for arms of interest
  for (i in 1:ncol(regdf)) { # convert arms of interest to binary depending on gain or loss
    if (eventsdf[i,2] == "Loss") { # if top event is loss for this arm, change binary for loss 
      regdf[,i] <- ifelse(regdf[,i] < -0.2, 1, 0)
    } else { # if top event is gain for this arm, change binary for gain
      regdf[,i] <- ifelse(regdf[,i] > 0.2, 1, 0)
    }
  }
  colnames(regdf) <- paste0(colnames(regdf), "_", eventsdf$V2) # fix names
  
  print("Arm event names correct: ")
  print(all(colnames(regdf) == events))
  
  ### check if events co-occur ######################## COME BACK ########################
  # comedf <- arm_co_me_filt[[tt]]
  # comedf_filt <- comedf[comedf$Event %in% c("Co_Occurence"),]
  # cooc_events <- data.frame()
  # for (i in events) {
  #   for (j in events) {
  #     if (nrow(comedf_filt[comedf_filt$Alt1 %in% i & comedf_filt$Alt2 %in% j,])>0) {
  #       print(paste0(i, " and ", j, " ", comedf_filt[comedf_filt$Alt1 %in% i & comedf_filt$Alt2 %in% j,"Event"]))
  #       cooc_events <- rbind(cooc_events, data.frame(Alt1 = i, Alt2 = j))
  #     }
  #   }
  # }
  # ### if cooc, add to regdf
  # if (nrow(cooc_events)>0) {
  #   for (i in 1:nrow(cooc_events)) {
  #     regdf[[ paste0(cooc_events[i,"Alt1"], "_and_", cooc_events[i,"Alt2"]) ]] <- 
  #       ifelse(regdf[,cooc_events[i,"Alt1"]] == 1 & regdf[,cooc_events[i,"Alt2"]] == 1, 1, 0 )
  #   }
  #   regdf[,unique(c(cooc_events$Alt1, cooc_events$Alt2))] <- NULL ######## remove orig columns
  # }
  
  
  ### add IS to model df
  regdf <- merge(isdf[,"Ranked_Sum", drop=F], regdf, by=0) %>% column_to_rownames("Row.names") # add immune score
  regdf$IS <- c(scale(regdf$Ranked_Sum)) # scale IS
  regdf$IS <- ifelse(regdf$IS >= unname(quantile(regdf$IS, probs = 0.7)), 1,  # change to binary
                     ifelse(regdf$IS <= unname(quantile(regdf$IS, probs = 0.3)), 0, NA))
  
  ### run multi reg model
  # regform <- as.formula(paste0("IS ~ ", reg_con_hits_tt[reg_con_hits_tt$TumorType == tt, "TopHits_Formula"])) # regression formula
  regform <- as.formula(paste0("IS ~ ", paste(grep("Arm", names(regdf), value = T), collapse = "+")))
  print(regform)
  
  glm_out <- glm(regform, data = regdf, family = binomial(link = "logit")) # run logistic model
  saveglm <- as.data.frame(summary(glm_out)$coefficients)[-1,]
  if (any(is.na(coef(glm_out)))) { # if any coef NA, add back to table
    saveglm <- merge(saveglm, data.frame(Estimate = coef(glm_out)[-1]), by=0, all=T)[,-6] %>%
      column_to_rownames("Row.names") %>% rename("Estimate" = "Estimate.x")
    saveglm <- saveglm[grep("Arm", names(regdf), value = T),] # order arms
  } 
  saveglm <- cbind(TumorType=rep(tt), saveglm) %>% rownames_to_column("Event")
  multi_reg_out <- rbind(multi_reg_out, saveglm)
}


colnames(multi_reg_out) <- c("Event", "TumorType", "Coefficient", "StdError", "Zval", "Pval")
write.csv(multi_reg_out, "IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/continuous_filttumors_tophits_MULTI_REGRESSION.csv")
# write.csv(multi_reg_out, "IS_CNV_Regression/TCGA/Arm/continuous_filttumors_tophits_CoOccur_MULTI_REGRESSION.csv")


##### plot results
multi_reg_out <- read.csv("IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/continuous_filttumors_tophits_MULTI_REGRESSION.csv", row.names = 1)
multi_reg_out <- multi_reg_out[!multi_reg_out$TumorType %in% "HNSC.HPVother",]

heatmaplist <- list()
# tt=unique(multi_reg_out$TumorType)[2]
for (tt in unique(multi_reg_out$TumorType)) {
  test <- multi_reg_out[multi_reg_out$TumorType %in% tt,]
  test_mat <- as.matrix(data.frame(row.names = gsub("_", " ", gsub("Arm", "", test$Event)), Z=test$Zval))
  
  hm <- Heatmap(test_mat, col = colorRamp2(c(-4, 0, 2), c("#521a63", "white", "#61631a")), 
                show_heatmap_legend = FALSE, rect_gp = gpar(col = "black", lwd = 1), 
                show_column_names = F, show_row_names = T, cluster_rows = F, cluster_columns = F, 
                row_names_side = "left", column_title = tt, 
                heatmap_width = unit(3.5, "cm"), heatmap_height = unit(2*nrow(test_mat), "cm"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1f", test_mat[i, j]), x, y, gp = gpar(fontsize = 10, col="white"))
                })
  heatmaplist[[tt]] <- grid.grabExpr(draw(hm))
}

pdf("IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/continuous_filttumors_tophits_MULTI_REGRESSION_Heatmaps.pdf", width = 14, height = 14)
grid.arrange(heatmaplist$ACC, heatmaplist$BLCA, heatmaplist$BRCA.pos, heatmaplist$CESC, heatmaplist$COADREAD.MSI, 
             heatmaplist$HNSC.HPVneg, heatmaplist$KIRC, heatmaplist$KIRP, heatmaplist$LGG, heatmaplist$LIHC, 
             heatmaplist$LUAD, heatmaplist$LUSC, heatmaplist$OV, heatmaplist$SARC, heatmaplist$SKCM, heatmaplist$STAD, 
             heatmaplist$TGCT, heatmaplist$UVM, ncol = 9, nrow=2) 
dev.off()







