########################################################################
#clean.workplace
########################################################################
rm(list = ls())
dev.off()
cat("\014")
########################################################################
#Library and Function
########################################################################

library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(readxl)
library(parallel)
library(forcats) 
library(Hmisc)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(cowplot)

`%ni%` <- Negate(`%in%`)
set.seed(1234567)

########################################################################
#Obtained health score and its related phenotype 
#from /data_preparation/03-health-score
#use the result directly
########################################################################
phenotypic_combine      <- openxlsx::read.xlsx("./Data_upload/Phenotypic data.xlsx")

rownames(phenotypic_combine) <- phenotypic_combine$Strain_id
Strain_id               <- phenotypic_combine[,c("Strain_id", "strain", "diet")]

data.plot               <- phenotypic_combine[, c("Strain_id", "strain", "diet", "MHS")]
data.plot               <- data.plot[!is.na(data.plot$MHS), ]
data.plot.tt            <- split(data.plot, data.plot$strain)

strain_orders           <- names(data.plot.tt)
BXD_strains             <- strain_orders[grepl("BXD", strain_orders)]
BXD_strains             <- BXD_strains[order(as.numeric(gsub("BXD|a|b", "", BXD_strains)), decreasing = F)]
BXD_strains             <- c(BXD_strains, strain_orders[!grepl("BXD", strain_orders)])

P.all                   <- lapply(BXD_strains, function(BXD_strain){ #tt= data.plot.tt[[3]]
  
  tt                    <- data.plot.tt[[BXD_strain]]
  print(BXD_strain)
  
  p <- ggplot(tt, aes(x=diet, y=MHS, color=diet, fill = diet)) +
    geom_boxplot(alpha = 0.5, outlier.color = "white") +
    #facet_grid(cols = vars(name), scales = "free", space = "free") +
    geom_jitter(aes(colour = diet)) +
    stat_boxplot(geom= 'errorbar' , width = 0.3, position = position_dodge(width = 0.75) ) +
    stat_compare_means(comparisons = list(c("CD", "HFD")), label = "p.signif",position = "identity", method = "t.test") +
    theme_bw()+
    theme(axis.text=element_text(size=8,color="black"),
          plot.title=element_text(face = "bold",size=11,hjust = 0.5,vjust = 0, margin=margin(0,0,5,0)),
          plot.subtitle = element_text(size=10,hjust = 0.5, vjust = 0, margin=margin(0,0,5,0)),
          axis.title=element_text(size=5,hjust = 0.5),
          legend.title = element_blank(),
          legend.text=element_text(size=11),
          strip.text = element_text(face = "bold",size=10),
          strip.background = element_blank(),
          #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
          element_line(color = "gray25", size = 0.5),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          panel.spacing = unit(0, "lines"),
          panel.grid =element_blank(),legend.position = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
    scale_colour_manual(values = c("CD" = "#FF9300", "HFD" = "#107F40")) +
    scale_fill_manual(values = c("CD" = "#FF9300", "HFD" = "#107F40")) +
    labs(x="",
         y="") + ggtitle(unique(tt$strain)) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
})

p.all.plot             <- plot_grid(plotlist = P.all, align = "hv", nrow = 7)
