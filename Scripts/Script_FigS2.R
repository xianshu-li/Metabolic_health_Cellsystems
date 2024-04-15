####################################
#clean.workplace
####################################

rm(list = ls())

####################################
#Library
####################################

library(tidyverse)
library(data.table)
library(ggpubr)

`%ni%` <- Negate(`%in%`)


######Read phenotypic data from a CC founders######

phenotypic_combine

#select the components of MHS
risk_factors            <- c("Fat_mass_perc", "ogtt_glycemia_t0", "ogtt_insulin_t0", "plasma_Cholesterol", "plasma_TGL")
pheno                   <- phenotypic_combine[, risk_factors]
rownames(pheno)         <- phenotypic_combine$ID
#only use the mouse with all components of MHS
pheno.filter            <- pheno[rowSums(is.na(pheno)) < 1, ]

#calculate MHS
Zscore_for              <- function(x){ 
  z <- (x - mean(x)) / sd(x)
  return(z)
}
counttemp2              <- data.frame(apply(pheno.filter, 2,Zscore_for))
counttemp               <- data.frame(apply(counttemp2, 1, sum))
colnames(counttemp)     <- "MHS"
counttemp$MHS           <- counttemp$MHS * (-1)
counttemp$Strain_id     <- rownames(counttemp)

#merge MHS into the previous phenotype data table
phenotypic_combine      <- merge(phenotypic_combine,counttemp, by.x = "ID", by.y = "Strain_id")
phenotypic_combine$strain_diet <- paste0(phenotypic_combine$StrainShort, "_", phenotypic_combine$Diet)

#use the mean value of MHS in WD condition to order strains

WD_phenotype           <- phenotypic_combine[phenotypic_combine$Diet == "WD", c("ID", "Strain", "MHS")]

mean_value             <- do.call(rbind, lapply(split(WD_phenotype, WD_phenotype$Strain), function(x){
  tmp                  <- data.table(Strain = unique(x$Strain), mean = mean(x$MHS, na.rm = T))
}))


order                   <- mean_value$Strain[order(mean_value$mean, decreasing = T)]
phenotypic_combine$Strain <- factor(phenotypic_combine$Strain, levels = order)

###FigS2A
p <- ggplot(phenotypic_combine, aes(x=Diet, y=MHS, color= Diet)) +
  geom_boxplot(outlier.color = "white", fill = NA) +
  geom_jitter() +
  scale_colour_manual(values = c("CD" = "#FF9300", "WD" = "#107F40")) +
  stat_boxplot(geom= 'errorbar' , width = 0.3, position = position_dodge(width = 0.75) ) +
  facet_grid(cols = vars(Strain), scales = "free", space = "free") +
  stat_compare_means(comparisons = list(c("CD", "WD")), label = "p.signif",position = "identity",method = "t.test") +
  theme_bw()+
  theme(axis.text.y =element_text(size=8,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_text(face = "bold",size=11,hjust = 0.5,vjust = 0, margin=margin(0,0,5,0)),
        plot.subtitle = element_text(size=10,hjust = 0.5, vjust = 0, margin=margin(0,0,5,0)),
        axis.title=element_text(size=5,hjust = 0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=11),
        strip.text = element_text(face = "bold",size=10),
        strip.background = element_blank(),
        #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
        element_line(color = "gray25", size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.spacing = unit(0, "lines")) +
  labs(x="",
       y="MHS") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))



###FigS2B
#select the NASH-related parameters
traits <- c("Steatosis_score", "Lobular_inflammation_score", "NAS_score")

p_all_traits <- lapply(traits, function(trait){ #trait = traits[1]
  print(trait)
  plot.tmp   <- phenotypic_combine[, c("Diet", trait, "Strain", "StrainColors")]
  colnames(plot.tmp) <- c("Diet", "trait", "Strain", "StrainColors")
  #keep the same order as the MHS
  plot.tmp$Strain <- factor(plot.tmp$Strain, levels = order)
  
p <- ggplot(plot.tmp, aes(x=Diet, y=trait, color= Diet)) +
  geom_boxplot(outlier.color = "white", fill = NA) +
  geom_jitter() +
  scale_colour_manual(values = c("CD" = "#FF9300", "WD" = "#107F40")) +
  stat_boxplot(geom= 'errorbar' , width = 0.3, position = position_dodge(width = 0.75) ) +
  facet_grid(cols = vars(Strain), scales = "free", space = "free") +
  stat_compare_means(comparisons = list(c("CD", "WD")), label = "p.signif",position = "identity",method = "t.test") +
  theme_bw()+
  theme(axis.text.y =element_text(size=8,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_text(face = "bold",size=11,hjust = 0.5,vjust = 0, margin=margin(0,0,5,0)),
        plot.subtitle = element_text(size=10,hjust = 0.5, vjust = 0, margin=margin(0,0,5,0)),
        axis.title=element_text(size=5,hjust = 0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=11),
        strip.text = element_text(face = "bold",size=10),
        strip.background = element_blank(),
        #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
        element_line(color = "gray25", size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.spacing = unit(0, "lines")) +
  labs(x="",
       y= trait) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))


})


