library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

cytoband_order <- read.csv("data/cytobandorder.csv") # order of cytobands

plot_cnv_frequency_area <- function(df, thresh=0.2, plottitle=NULL, savepath=NULL, 
                                    ylabel=NULL, output=c("plot", "table", "both")) {
  
  totaln <- nrow(df) # num samples
  df_cat <- as.data.frame(apply(df, 2, function(x) ifelse(x < -thresh, "Loss", ifelse( x > thresh, "Gain", "None")))) # convert to binary
  df_freq <- as.data.frame(t(sapply(df_cat, function(x) table(factor(x, levels = c("Loss", "None", "Gain")))))) # get frequencies
  
  colnames(df_freq) <- paste0(colnames(df_freq), "_Count") # colnames
  df_freq <- df_freq %>% mutate(Gain_Freq = Gain_Count/totaln, Loss_Freq = Loss_Count/totaln) # gain and loss frequencies
  df_freq$Cytoband <- factor(rownames(df_freq), levels = cytoband_order$Cytoband) # cytobands
  df_freq$Chromosome <- factor(gsub("p.*|q.*", "", df_freq$Cytoband), levels = c(1:22, "X", "Y")) # add chr columns
  df_freq <- df_freq[order(df_freq$Cytoband),] # order cytobands
  
  ### df for plotting
  df_plot <- pivot_longer(df_freq[,c(6,7,4,5)], cols = c(Gain_Freq, Loss_Freq), names_to = "Event", values_to = "Frequency") 
  df_plot$Event <- gsub("_Freq", "", df_plot$Event)
  df_plot$ChrArm <- paste0(df_plot$Chromosome, ifelse(grepl("p", df_plot$Cytoband), "p", "q"))
  df_plot$ChrArm <- factor(df_plot$ChrArm, levels = unique(df_plot$ChrArm))
  df_plot$Frequency <- ifelse(df_plot$Event == "Loss", df_plot$Frequency*-1, df_plot$Frequency)
  
  ### for x axis
  chr.cut = na.omit(df_plot %>% group_by(Ref = Chromosome) %>% summarise(vals = last(Cytoband))) # end of chr
  
  arm.cut <- df_plot %>% group_by(Ref = ChrArm) %>% summarise(vals = last(Cytoband)) # end of arms
  arm.cut$vals <- ifelse(grepl("q", arm.cut$Ref), NA, as.character(arm.cut$vals)) # remove q arm (don't add line to end of q)
  arm.cut$Ref <- gsub("p|q", "", arm.cut$Ref)
  arm.cut <- arm.cut[!duplicated(arm.cut$Ref),]
  
  chr.mid <- df_plot[,c("Chromosome", "Cytoband")] %>% 
    group_by(Chromosome) %>% slice(n() %/% 2 + 1) %>% ungroup() %>% rename(Ref=Chromosome, vals=Cytoband)
  
  p <- ggplot(df_plot) +
    geom_area(aes(x = Cytoband, y = Frequency, group = Event, col = Event, fill=Event), position = "identity" ) +
    scale_fill_manual(values = alpha(c("#ad0224","#020dad"),0.2)) +
    scale_color_manual(values = alpha(c("#ad0224","#020dad"))) +
    labs(title = plottitle, y=ylabel) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = margin(0, 0, 0, 0, "cm")) +
    scale_x_discrete(limits = unique(df_plot$Cytoband)) + 
    ylim(-1, 1) +
    geom_vline(xintercept = arm.cut$vals, linetype="dotted", color = "#8F8F8F", linewidth=0.2) +
    geom_vline(xintercept = chr.cut$vals, linetype="longdash", color = "#8F8F8F", linewidth=0.3) +
    geom_hline(yintercept = c(-1, 0,1), color = "black", linewidth=0.2) +
    geom_hline(yintercept = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), linetype="dashed", color = "#8F8F8F", linewidth=0.2) +
    geom_text(data = chr.mid, mapping = aes(x = vals, y = -0.9, label = Ref, hjust = 0.5, vjust = 0), angle=0, size=3)
  
  
  if (!is.null(savepath)) {
    ### save pdf
    pdf(paste0(savepath, "_cytoband_freq.pdf"), width = 10, height = 5)
    print(p)
    junk <- dev.off()
    
    ### save frequency table
    write.csv(df_freq, paste0(savepath, "_cytoband_freq.csv"))
  }
  
  if (output == "plot") { return(p) } 
  if (output == "table") { return(df_plot) }
  if (output == "both") { return(list(plot=p, table=df_plot)) }
}
