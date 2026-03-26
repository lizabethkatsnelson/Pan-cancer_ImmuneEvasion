library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

### function to filter summary tables based on dels or amps
custom_filter <- function(df, filter_del=list(), filter_amp=list(), plottitle="", savepath=NULL, 
                          topN = NULL, rankorder = list(del=c(), amp=c()) ) {
  totalgenes <- nrow(df) # original number genes
  
  ### deletions
  df_del <- df # new variable that will be used for just deletion filtering
  filtered_deletions <- data.frame(Deletions_Step = "Full Table", Deletions_NumberGenes = totalgenes)
  for (step in filter_del) { # iterate through all criteria for deletions, save how many rows (genes)
    df_del <- df_del %>% filter(eval(parse(text = step))) # take filter criteria for this step and filter table
    filtered_deletions <- rbind(filtered_deletions, data.frame(Deletions_Step = step, Deletions_NumberGenes = nrow(df_del))) # save number of rows per iteration
  }
  
  ### amplifications
  df_amp <- df # new variable that will be used for just amp filtering
  filtered_amplifications <- data.frame(Amplifications_Step = "Full Table", Amplifications_NumberGenes = totalgenes)
  for (step in filter_amp) { # iterate through all criteria for amplifications, save how many rows (genes)
    df_amp <- df_amp %>% filter(eval(parse(text = step))) # take filter criteria for this step and filter table
    filtered_amplifications <- rbind(filtered_amplifications, data.frame(Amplifications_Step = step, Amplifications_NumberGenes = nrow(df_amp))) # save number of rows per iteration
  }
  
  ## combine stats output
  gene_counts_output <- cbind(filtered_deletions, filtered_amplifications)
  
  ## genes in common and unique between del and amp after filtering
  incommon <- intersect(df_del$Gene, df_amp$Gene)
  only_deletions <- setdiff(df_del$Gene, df_amp$Gene)
  only_amplifications <- setdiff(df_amp$Gene, df_del$Gene)
  filtered_genes <- list(Both_Amp_Del = incommon, Only_Deletions = only_deletions, Only_Amplifications = only_amplifications)
  
  ## add only amp and only del gene counts to filtered counts table
  gene_counts_output <- rbind(gene_counts_output, data.frame(Deletions_Step="Deletions only", Deletions_NumberGenes=length(only_deletions), 
                                                             Amplifications_Step="Amplifications only", Amplifications_NumberGenes=length(only_amplifications)))
  
  ## topN (if ranking liptak)
  if (!is.null(topN)) {
    # del
    df_del <- df_del %>% arrange(eval(parse(text = rankorder$del))) %>% slice(1:topN) # order by given variables, take top N
    # amp
    df_amp <- df_amp %>% arrange(eval(parse(text = rankorder$amp))) %>% slice(1:topN)
  }
  
  ## PLOTTING
  # deletions
  plot_del <- data.frame(Gene = df_del$Gene, Chromosome = factor(gsub("p.*|q.*", "", df_del$Cytoband), levels = c(1:22, "X", "Y")),
                         Arm = factor(ifelse(grepl("p", df_del$Cytoband), "p", "q"), levels = c("p", "q")))
  plot_del <- plot_del[plot_del$Gene %in% only_deletions,] ########### filter for only deletion hits
  plot_del <- plot_del[order(plot_del$Chromosome, plot_del$Arm),] %>% 
    group_by(Chromosome, Arm) %>% summarise(Count = n()) # group by chromosome arm and count genes
  plot_del$Arm <- paste0(plot_del$Arm, "_deletions") # add "deletions" to arm labels for plot colors
  plot_del$Count <- plot_del$Count * -1 # deletion counts will be negative in plot
  
  # amplifications
  plot_amp <- data.frame(Gene = df_amp$Gene, Chromosome = factor(gsub("p.*|q.*", "", df_amp$Cytoband), levels = c(1:22, "X", "Y")),
                         Arm = factor(ifelse(grepl("p", df_amp$Cytoband), "p", "q"), levels = c("p", "q")))  
  plot_amp <- plot_amp[plot_amp$Gene %in% only_amplifications,] ########### filter for only amps hits
  plot_amp <- plot_amp[order(plot_amp$Chromosome, plot_amp$Arm),] %>%
    group_by(Chromosome, Arm) %>%
    summarise(Count = n()) # counts for amp will be positive
  plot_amp$Arm <- paste0(plot_amp$Arm, "_amplifications")
  
  # combine
  plot_all <- rbind(plot_amp, plot_del) # combine del and amp tables
  plot_all$Arm <- factor(plot_all$Arm, levels = c("p_amplifications", "q_amplifications", "p_deletions", "q_deletions"))
  
  ## plot and save
  plotcolors <- c("p_amplifications"="#B74E4E", "q_amplifications"="#4a0101", "p_deletions"="#4EB7A5", "q_deletions"="#13473e")
  p = ggbarplot(plot_all, x = "Chromosome", y = "Count", fill = "Arm", 
                xlab = "Chromosome", ylab = "Gene Count", 
                lab.size = 11, width = 0.85, title = plottitle,
                color = NA, palette = plotcolors) + theme_bw()
  if (!is.null(savepath)) {
    
    write.csv(gene_counts_output, paste0(savepath, "_filtered_gene_counts.csv"), row.names = F) # save gene counts table
    
    pdf(paste0(savepath, "_filtered_genes_distribution_plot.pdf"), width = 6, height = 4) # save plot
    print(p)
    junk <- dev.off()
    
    if (nrow(df_del) > 0 & nrow(df_amp) > 0) { # if del and amp hits
      save_gene_hits_summary <- rbind(cbind(Filtered_Hit = rep("Deletion"), df_del),
                                      cbind(Filtered_Hit = rep("Amplification"), df_amp))
    } else if (nrow(df_del) > 0 & nrow(df_amp) == 0) { # if only del hits (no amp genes found)
      save_gene_hits_summary <- cbind(Filtered_Hit = rep("Deletion"), df_del)
    } else if (nrow(df_del) == 0 & nrow(df_amp) > 0) { # if only amp hits (no del genes found)
      save_gene_hits_summary <- cbind(Filtered_Hit = rep("Amplification"), df_amp)
    } else {
      save_gene_hits_summary <- NULL
    }
    
    if (!is.null(save_gene_hits_summary)) { # save the table if it exists (otherwise there were no hits found)
      write.csv(save_gene_hits_summary, paste0(savepath, "_filtered_summary_table.csv"), row.names = F)
    }
    
  }
  
  return(list(Filtered_Deletion_Table = df_del, 
              Filtered_Amplification_Table = df_amp, 
              Gene_Counts = gene_counts_output, 
              Gene_Lists = filtered_genes,
              Distribution_Plot = p))
}
