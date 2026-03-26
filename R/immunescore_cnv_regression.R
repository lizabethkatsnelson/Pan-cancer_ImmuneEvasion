library(dplyr)
library(tidyverse)

### gene locations
gene_locations <- read.csv("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/scripts/Biomart_Gene_Locations.csv")

### function to run immune score and cnv regressions
immunescore_cnv_regression <- function(cnv_df, is_df, use_AS = T, cnv_thresh=0.2, is_thresh=0.3, 
                                       CNV_binary=T, IS_binary=T, model=c("log", "lin"), savepath=NULL, 
                                       cnvlevel=c("Gene", "Cytoband", "Arm")) {
  ### preprocess
  isname = names(is_df) # get IS name
  inputdf <- merge(is_df, cnv_df, by=0) %>% remove_rownames() %>% column_to_rownames(var="Row.names") # add IS to cnv table
  inputdf[,1] <- scale(inputdf[,1]) # scale IS
  
  ### set formula(s)
  form_gain <- as.formula(paste0(isname, " ~ Gain"))
  form_loss <- as.formula(paste0(isname, " ~ Loss"))
  if (use_AS) {
    inputdf[,2] <- scale(inputdf[,2]) # scale AS
    asname = names(inputdf)[2] # AS column name
    form_gain <- update(form_gain, paste(". ~ . +", asname))
    form_loss <- update(form_loss, paste(". ~ . +", asname))
  }
  
  ### convert IS to binary or keep as continuous
  if (IS_binary) {
    inputdf[,1] <- ifelse(inputdf[,1] >= unname(quantile(inputdf[,1], probs = 1-is_thresh)), 1,
                          ifelse(inputdf[,1] <= unname(quantile(inputdf[,1], probs = is_thresh)), 0, NA))
  }
  
  ### get model family
  if (model == "log") { family = binomial(link = "logit") } # logistic
  if (model == "lin") { family = gaussian } # linear
  
  ### iterate through all genes or cytobands and run models
  out <- data.frame() # output df
  colstart <- ifelse(use_AS, 3, 2) # which column do gene or cytoband cnv values start
  
  for (col in colstart:ncol(inputdf)) { # iterate thru all genes/cytobands/arms to run model
    colname=names(inputdf)[col] # gene or cytoband name
    
    ### get df for running model
    if (use_AS) { 
      df <- inputdf[,c(1:2,col)] # col1 = IS, col2 = AS, col3 = Gene/Cyto/Arm
    } else { df <- inputdf[,c(1,col)] } # col1 = IS, col2 = Cyto
    
    ### convert CNV to binary or keep as continuous
    if (CNV_binary) {
      df$Gain <- ifelse(df[,colname] > cnv_thresh, 1, 0) # gain = 1, all else = 0
      df$Loss <- ifelse(df[,colname] < -cnv_thresh, 1, 0) # loss = 1, all else = 0
    } else {
      df$Gain <- ifelse(df[,colname] > -cnv_thresh, df[,colname], NA) # if cnv greater than lower threshold, keep
      df$Loss <- ifelse(df[,colname] < cnv_thresh, df[,colname], NA) # if cnv lower than higher threshold, keep
    }
    
    ### MODEL FOR GAINS
    if ( sum(abs(df$Gain), na.rm = T) > 3 ) { # at least 4 samples with gain
      glm_gain <- glm(form_gain, data = df, family = family) # run gain model
      if ( any(is.na(unname(glm_gain$coefficients))) == F ) { # if there are no NAs in the output:
        glm_gain_res <- data.frame(summary(glm_gain)$coefficients['Gain',, drop=F]) # get results
        rownames(glm_gain_res) <- paste0(colname, "_Gain") # save in table
        out <- rbind(out, glm_gain_res)
      }
    }
    
    ### MODEL FOR LOSSES
    if ( sum(abs(df$Loss), na.rm = T) > 3 ) { # at least 4 samples with loss
      glm_loss <- glm(form_loss, data = df, family = family) # run loss model
      if ( any(is.na(unname(glm_loss$coefficients))) == F ) {
        glm_loss_res <- data.frame(summary(glm_loss)$coefficients['Loss',, drop=F])  # get results
        rownames(glm_loss_res) <- paste0(colname, "_Loss") # save in table
        out <- rbind(out, glm_loss_res)
      }
    }
  }
  
  ### process outputs
  if (model == "log") { colnames(out) <- c("Coefficients", "Std_Error", "Z_val", "P_val") } # fix column names
  if (model == "lin") { colnames(out) <- c("Coefficients", "Std_Error", "T_val", "P_val") }
  out$Name <- gsub("_.*", "", rownames(out)) # get gene or cyto name
  out$Event <- gsub(".*_", "", rownames(out)) # get event
  
  ### add "zero" rows where cnv didnt exist (fixes plot issues)
  gains <- out[out$Event == "Gain", "Name"]; losses <- out[out$Event == "Loss", "Name"] # get all genes for gain and loss
  addtogain <- setdiff(losses, gains); addtoloss <- setdiff(gains, losses) # any that don't intersect? aka missing because model couldnt run
  
  if ( length(addtogain) > 0 ) {
    addtogain_df <- data.frame(rep(NA, length(addtogain)), rep(NA, length(addtogain)), rep(0, length(addtogain)),
                               rep(1, length(addtogain)), addtogain, rep("Gain", length(addtogain)))
    colnames(addtogain_df) <- colnames(out)
    out <- rbind(out, addtogain_df)
  }
  
  if ( length(addtoloss) > 0 ) {
    addtoloss_df <- data.frame(rep(NA, length(addtoloss)), rep(NA, length(addtoloss)), rep(0, length(addtoloss)),
                               rep(1, length(addtoloss)), addtoloss, rep("Loss", length(addtoloss)))
    colnames(addtoloss_df) <- colnames(out)
    out <- rbind(out, addtoloss_df)
  }
  
  ### add chromosome locations
  colnames(out)[5] <- cnvlevel
  if (cnvlevel == "Gene") {
    out_loc <- unique(merge(out, gene_locations, by="Gene", all.x=T))
    out_loc$Chromosome <- factor(out_loc$Chromosome, levels = paste0("Chr", c(1:22, "X", "Y"))) # set chromosome order
    out_loc <- out_loc[order(out_loc$Chromosome, out_loc$Start),]
  }
  if (cnvlevel == "Cytoband") {
    out_loc <- out
    out_loc$Chromosome <- factor(paste0("Chr", gsub("p.*|q.*", "", out_loc$Cytoband)), levels = paste0("Chr", c(1:22, "X", "Y")))
    out_loc$Arm <- ifelse(grepl("p", out_loc$Cytoband), "p", "q")
    out_loc$ChrArm <- paste0(out_loc$Chromosome, out_loc$Arm)
    out_loc$Cytoband <- factor(out_loc$Cytoband, levels = unique(out_loc$Cytoband))
    out_loc <- out_loc[order(out_loc$Cytoband),]
  }
  if (cnvlevel == "Arm") {
    out_loc <- out
    armfactors <- c(paste0(rep(1:12, each=2), c("p","q")), paste0(13:15, "q"), paste0(rep(16:20, each=2), c("p","q")), paste0(21:22, "q"), "Xp", "Xq")
    out_loc$Arm <- factor(gsub("Arm", "", out_loc$Arm), levels = armfactors)
    
    if (any(!armfactors %in% unique(out_loc$Arm))) { ### any arms missing? add zeros to tables
      armstoadd <- armfactors[!armfactors %in% unique(out_loc$Arm)]
      dftoadd <- data.frame(Coefficients=rep(NA), Std_Error=rep(NA), Z_val=rep(0), P_val=rep(1), 
                          Arm = sort(c(armstoadd, armstoadd)), Event=rep(c("Gain", "Loss")))
      out_loc <- rbind(out_loc, dftoadd)
    }
    out_loc$Chromosome <- factor(paste0("Chr", gsub("p.*|q.*", "", out_loc$Arm)), levels = paste0("Chr", c(1:22, "X")))
    out_loc$ChrArm <- ifelse(grepl("p", out_loc$Arm), paste0(out_loc$Chromosome, "p"), paste0(out_loc$Chromosome, "q"))
    out_loc <- out_loc[order(out_loc$Arm, out_loc$Event),]
  } 
  
  ### save and return
  if (!is.null(savepath)) {
    write.csv(out_loc, savepath, row.names = F)
  }
  
  return(out_loc)
}



