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
library(parallel)
library(Hmisc)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(cowplot)
library(FactoMineR)

`%ni%` <- Negate(`%in%`)

########################################################################
#Obtained health score and its related phenotype 
########################################################################
phenotypic_combine           <- openxlsx::read.xlsx("./Data_upload/Phenotypic data.xlsx")
rownames(phenotypic_combine) <- phenotypic_combine$Strain_id

which.col               <- c("MHS", "body_fat", "Blood_Glucose_.mmol.L.", "Resting_insulin", "Blood_Triglycerides_.mmol.L.", "Blood_TotalCholesterol_.mmol.L.")

result                  <- phenotypic_combine[, which.col]
result                  <- result[!is.na(result$MHS),]
result                  <- result[!is.na(result$Resting_insulin),]
result                  <- result[!is.na(result$Blood_Glucose_.mmol.L.),]
result                  <- result[!is.na(result$Blood_Triglycerides_.mmol.L.),]
result                  <- result[!is.na(result$body_fat),]
result                  <- result[!is.na(result$Blood_TotalCholesterol_.mmol.L.),]
result$Strain_id        <- rownames(result)
result_melt             <- reshape2::melt(result, id.var = "Strain_id")
result_melt             <- result_melt[order(result_melt$value, decreasing = T), ]
result_melt             <- merge(Strain_id, result_melt, by = "Strain_id")
result_melt$split       <- paste0(result_melt$strain, "_", result_melt$diet, "_", result_melt$variable)

data_lt                 <- split(result_melt, result_melt$split, drop = T)

data_ave                <- do.call(rbind, 
                                   lapply(data_lt, function(y){ #y= data_lt[["BXD39_HFD_Resting_insulin"]]
                                     ave      <- as.data.frame(mean(y$value))
                                     ave$sd   <- sd(y$value)
                                     colnames(ave) <- c("ave", "sd")
                                     ave$strain    <- unique(y$strain)
                                     ave$pheno     <- unique(y$variable)
                                     ave$diet      <- unique(y$diet)
                                     return(ave)
                                   }))


data     <- data_ave
data$id  <- paste0(data$strain, "_",data$diet)
data     <- reshape2::dcast(data, id ~ pheno, value.var = "ave"  )
data$diet <- gsub(".*_", "", data$id)
data$strain <- gsub("_.*", "", data$id)

####Figure1 B
#Z score
Zscore_for      <- function(x){
  z <- (x - mean(x)) / sd(x)
  return(z)
}
#need to do Z-score
data_ave$merge          <- paste0(data_ave$strain, "_", data_ave$diet)
data_ave                <- data_ave[order(data_ave$ave, decreasing = F),]
tt                      <- split(data_ave, data_ave$pheno)
data_ave                <- do.call(rbind, lapply(tt, function(lt){ #lt = tt[[1]]
  lt$zscore             <- Zscore_for(lt$ave)
  lt
}))

data_ave$label          <- paste0(gsub("_", " (", data_ave$merge), ")")

# Hierarchical clustering on phenotypes1
data_dcast <- reshape2::dcast(data_ave, merge~pheno, value.var="zscore")
data_dcast[is.na(data_dcast)] <- 0 # For the purpose of clustering, we consider "NA" equivalent to no correlation
rownames(data_dcast) <- data_dcast$merge
hc <- hclust(dist(t(data_dcast[,-1])))
rowInd <- hclust(dist(data_dcast[,-1]))$order
colInd <- hclust(dist(t(data_dcast[,-1])))$order

orders                  <- data_ave[data_ave$pheno == "MHS", ]
orders                  <- orders$label[order(orders$zscore, decreasing = T)]
data_ave$label          <- factor(data_ave$label, levels = orders)
data_ave$pheno          <- factor(data_ave$pheno, levels = c( "Blood_Triglycerides_.mmol.L.", "Blood_Glucose_.mmol.L.", "Blood_TotalCholesterol_.mmol.L.", "Resting_insulin", "body_fat", "MHS" ), 
                                  labels = c("Triglycerides", "Glucose", "Total Cholesterol", "Insulin", "Body fat", "Metabolic health score"))

Heatmap_palette         <- c(rev(brewer.pal(7,"Blues")), "white", brewer.pal(7,"Reds"))
g <- ggplot(data_ave, aes(x = label, y = pheno, col = zscore, fill = zscore)) + 
  geom_tile() +
  #facet_grid(cols = vars(diet) , space="free", scales="free") +
  scale_fill_gradientn(colours= Heatmap_palette , limits=c(-2,2) ,na.value="gray87", oob = squish) +
  scale_color_gradientn(colours= Heatmap_palette , limits=c(-2,2) ,na.value="gray87", oob = squish) +
  scale_y_discrete(expand=c(0, 0))+
  #coord_fixed(2) +
  theme_bw(base_size = 8) +
  ylab(NULL) +
  xlab(NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title = element_text(size = 8),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, color = "black", vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 8),
    legend.position = "right"
    #panel.border = element_blank()
  )

g2 <- ggplot(data_ave, aes(x = label, y = "diet", col = diet, fill = diet)) + 
  geom_tile() +
  scale_color_manual(values = c("CD" = "#FF9300", "HFD" = "#107F40")) + 
  scale_fill_manual(values = c("CD" = "#FF9300", "HFD" = "#107F40")) + 
  scale_y_discrete(expand=c(0, 0))+
  #coord_fixed(2) +
  theme_bw(base_size = 8) +
  ylab(NULL) +
  xlab(NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title = element_text(size = 8),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, color = "black", vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 8),
    legend.position = "right"
    #panel.border = element_blank()
  )

g_all   <- plot_grid(g,g2, align = "hv", nrow = 2, ncol =1, rel_heights = c(2.8,1.2))

########################################################################
#Use scatter plot to show Metabolic health score (Figure 1C)
########################################################################

Strain_id               <- phenotypic_combine[,c("Strain_id", "strain", "diet")]
which.col               <- c("Strain_id","strain","diet","MHS")
result                  <- phenotypic_combine[, which.col]
result                  <- result[!is.na(result$MHS),]
result$split            <- paste0(result$strain, "_", result$diet)

data_lt                 <- split(result, result$split, drop = T)

data_ave                <- do.call(rbind, 
                                   lapply(data_lt, function(y){ #y= data_lt[[1]]
                                     ave      <- as.data.frame(mean(y$MHS))
                                     ave$sd   <- sd(y$MHS)
                                     colnames(ave) <- c("ave", "sd")
                                     ave$strain    <- unique(y$strain)
                                     ave$diet      <- unique(y$diet)
                                     return(ave)
                                   }))

data_df                  <- reshape2::dcast(data_ave, strain~diet, value.var = "ave")

data.plot               <- phenotypic_combine[, c("Strain_id", "strain", "diet", "MHS")]
data.plot               <- data.plot[!is.na(data.plot$MHS), ]
data.plot.tt            <- split(data.plot, data.plot$strain)

Pvalue                  <- do.call(rbind, lapply(data.plot.tt, function(tt){ #tt= data.plot.tt[[3]]
  print(unique(tt$strain))
  if (length(unique(tt$diet)) > 1) {
    
    HFD                  <- tt[tt$diet == "HFD", ]
    CD                   <- tt[tt$diet == "CD", ]
    
    if(nrow(HFD) > 1 & nrow(CD) >1){
      P <- t.test(HFD$MHS, CD$MHS)
      out <- data.table(strain = unique(tt$strain),
                        P      = P$p.value)
    }else{
      out <- NULL
    }
  }else{
    out <- NULL
  }
  out
}))
tmp          <- Pvalue[Pvalue$P > 0.05, ]
#data_df$is_highlight     <- ifelse(abs(data_df$CD - data_df$HFD) < 1, "Yes", "No")
data_df$is_highlight     <- ifelse((data_df$strain %in% tmp$strain) & abs(data_df$CD - data_df$HFD) < 1, "Yes", "No")
data_df$label            <- data_df$strain
#data_df$label[data_df$is_highlight == "No"] <- NA

g <- ggplot(data = data_df, aes(x= CD, y = HFD, label = label, color = is_highlight)) + 
  geom_point() +
  scale_color_manual(values = c("Yes" = "orange", "No" = "#969696")) + 
  geom_text_repel(size=3, parse = F) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  theme_bw() + 
  theme(axis.text.x = element_text(size =10, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11),
        strip.background = element_rect( fill = "white"),
        panel.grid = element_blank(),
        #panel.border = element_blank()
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "right"
  )


########################################################################
#the correlation with other phenotype (Figure 1D)
########################################################################

which.col  <- c("MHS", "HOMA_IR", "OGTT_AUC_Glucose", "BW_gained", "scWAT.", "Liver.", "Blood_ALAT_.u.L.", "diet")
#conditions
Diets      <- c("HFD", "CD")
all(which.col %in% colnames(phenotypic_combine) )
data       <- phenotypic_combine[, colnames(phenotypic_combine) %in% which.col]
data       <- data[!is.na(data$MHS),]


lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = method,  se=FALSE) +
    stat_cor(
      aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
      method = "pearson"
    )
  p
}

g <- ggpairs(data, 
             mapping = aes(color = diet, group = diet),
             lower = list(continuous = wrap(lowerFn, method = "lm")),
             upper  = list(continuous = "blankDiag"),
             diag  = list(continuous = "blankDiag"),
             #insulin resistance
             columns = c("health_score_real","HOMA_IR", "OGTT_AUC_Glucose"),
             columnLabels = c("Metabolic health score", "HOMA (IR)", "Glucose (AUC)")) +
  
  # columns = c("health_score_real", "BW_gained", "scWAT."),
  # columnLabels = c("Metabolic health score", "BW (gained)", "scWAT (%)")) +
  # columns = c("health_score_real", "Liver.", "Blood_ALAT_.u.L."),
  # columnLabels = c("Metabolic health score", "Liver weight (%)", "ALT")) +
  scale_colour_manual("",values = c("CD" = "#FF9300", "HFD" = "#107F40")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
    #panel.border = element_blank()
  ) 
