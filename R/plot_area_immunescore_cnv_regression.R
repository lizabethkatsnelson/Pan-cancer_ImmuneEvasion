library(dplyr)
library(ggplot2)
library(ggpubr)


### Functino ot plot IS CNV regression model output
plot_cnv_IS_regression <- function(df, cnvlevel=c("Gene", "Cytoband", "Arm"), scorename="IS", 
                                   stat=c("P_val", "Z_val", "T_val"), cnvtype=c("bi", "con"),
                                   signif=2, plottitle="", savepath=NULL) {
  
  ### process table
  df <- df[!is.na(df$Chromosome),] # remove rows with missing chrom info
  df$Chromosome <- factor(gsub("Chr", "", df$Chromosome), levels = c(1:22, "X", "Y"))
  df$ChrArm <- factor(gsub("Chr", "", df$ChrArm), levels = na.omit(unique(gsub("Chr", "", df$ChrArm))))
  
  if (cnvlevel == "Gene") {
    df <- df[order(df$Chromosome, df$Start),] # order by chr and start site
  }
  if (cnvlevel == "Cytoband") {
    # df <- df[order(df$Cytoband),] # order by cytoband
    df$Cytoband <- factor(df$Cytoband, levels = unique(df$Cytoband))
  }
  if (cnvlevel == "Arm") {
    df <- df[order(df$ChrArm),] # order by arm
  }
  
  # plotting with p val and binary cnv
  if (stat == "P_val" & cnvtype == "bi") { 
    df$signedlogP <- ifelse(df$Coefficients > 0, -log10(df$P_val), -log10(df$P_val)*-1 ) # if coef pos, logP is pos, if neg logP is neg
    df$signedlogP[is.na(df$signedlogP)] <- 0 # change NA to 0
    stat = "signedlogP" # change which var to use for plot
  }
  
  # plotting with pval and continuous cnv
  if (stat == "P_val" & cnvtype == "con") { 
    df$signedlogP <- -log10(df$P_val)
    df$signedlogP <- ifelse(df$Event == "Gain" & df$Coefficients < 0, df$signedlogP *-1, # if gain and (-) weight, logP (-)
                            ifelse(df$Event == "Loss" & df$Coefficients > 0, df$signedlogP *-1, df$signedlogP)) # if loss and (+) weight, logP (-)
    df$signedlogP[is.na(df$signedlogP)] <- 0 # change NA to 0
    stat = "signedlogP" # change which var to use for plot
  }
  
  # plotting with z val (log) or t val (lin) and binary cnv - don't need to change anything, statistic already correct direction
  
  # plotting with z or t(lin) val and continuous cnv
  if (stat %in% c("Z_val", "T_val") & cnvtype == "con") { ### because stat already signed, just need to change direction for losses
    df[,stat] <- ifelse(df$Event == "Loss", df[,stat]*-1, df[,stat]) # if loss, change dir of stat (neg becomes pos, pos becomes neg)
  }

  ### get x axis cutoffs for genome
  chr.cut = na.omit(df %>% group_by(Ref = Chromosome) %>% summarise(vals = last(get(cnvlevel)))) # end of chr

  arm.cut <- df %>% group_by(Ref = ChrArm) %>% summarise(vals = last(get(cnvlevel))) # end of arms
  arm.cut$vals <- ifelse(grepl("q", arm.cut$Ref), NA, as.character(arm.cut$vals)) # remove q arm (don't add line to end of q)
  arm.cut$Ref <- gsub("p|q", "", arm.cut$Ref)
  arm.cut <- arm.cut[!duplicated(arm.cut$Ref),]
  
  chr.mid <- df %>% group_by(Ref = Chromosome) %>% filter(row_number()==ceiling(n()/2)) %>% dplyr::select(Ref, vals=cnvlevel) # middle of each chr
  
  y = sym(stat) ## turn string into symbol 
  x = sym(cnvlevel)
  
  ymax = ceiling(max(df[,stat]))
  ymin = floor(min(df[,stat]))
  
  ymax = max(abs(ymin), ymax) 
  if (ymax < 4) { ymax = 4 }
  
  if (cnvlevel %in% c("Gene", "Cytoband")) {  
    p <- ggplot(df) +
      geom_area(aes(x = !!x, y = !!y, group = Event, fill = Event), position = "identity") + 
      scale_fill_manual(values = alpha(c("#ad0224","#020dad"),0.6)) +
      labs(title = plottitle) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            title = element_text(size=10), axis.text.y = element_text(size=9)) +
      scale_x_discrete(limits = unique(df[,cnvlevel])) + 
      ylim(-ymax, ymax) +
      geom_vline(xintercept = arm.cut$vals, linetype="dotted", color = "#8F8F8F", linewidth=0.15) +
      geom_vline(xintercept = chr.cut$vals, linetype="longdash", color = "#8F8F8F", linewidth=0.3) +
      geom_hline(yintercept = c(-signif, signif), linetype="dashed", color = "black", linewidth=0.3) +
      geom_hline(yintercept = c(-1, 1), linetype="dashed", color = "black", linewidth=0.15) +
      geom_hline(yintercept = 0, color = "black", linewidth=0.15) +
      geom_text(data = chr.mid, mapping = aes(x = vals, y = -ymax, label = Ref, hjust = 0, vjust = 0), angle=90, size=2)
  }
  
  if (cnvlevel == "Arm") {
    p <- ggplot() +
      geom_bar(data = df[df$Event == "Loss",], aes(x = !!x, y = !!y, fill = Event), stat = "identity", width = 1) + 
      geom_bar(data = df[df$Event == "Gain",], aes(x = !!x, y = !!y, fill = Event), stat = "identity", width = 1) + 
      scale_fill_manual(values = alpha(c("#ad0224","#020dad"),0.6)) +
      labs(title = plottitle) + ylim(-ymax, ymax) + theme_bw()
  }
  
  ### save
  if (is.null(savepath)) { return(p) } else { # if no path given to save pdf, return the plot as output for function
    pdf(savepath, width = 10, height = 5)
    print(p)
    junk <- dev.off()
  }
}