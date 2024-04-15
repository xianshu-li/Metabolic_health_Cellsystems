
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
library(WGCNA)
library(FactoMineR)

`%ni%` <- Negate(`%in%`)


########################################################################
#WGCNA
########################################################################

clusterWGCNA <- function(dataset, powers, analysis, organ){ #dataset = datExpr analysis = diet
  print("running pickSoftThreshold")
  set.seed(20000001)
  sft <- pickSoftThreshold(t(dataset), 
                           powerVector = powers, 
                           verbose = 5, 
                           networkType = "signed hybrid", 
                           blockSize = 25000, 
                           corFnc = "bicor",corOptions = list(quick = 0.5, use = 'pairwise.complete.obs'),
                           moreNetworkConcepts = T)
  
  pdf(paste0(file_out, organ, "_", analysis, "_09_module_generation_pickSoftThreshold.pdf"), width = 10, height = 10)
  cex1 = 0.9;
  
  par(mfrow = c(1,2));
  plot(sft$fitIndices[,1], sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  dev.off()    
  
  fit.sequence <- sft$fitIndices[,"SFT.R.sq"]
  names(fit.sequence) <- sft$fitIndices[,"Power"]
  fit.sequence[1] <- 0 #never take the first power
  bestpower <-  as.numeric(ifelse(max(fit.sequence)> 0.85, names(fit.sequence)[min(which(fit.sequence > 0.85))], names(fit.sequence)[which.max(fit.sequence)]))
  minmodsize <- 30
  
  set.seed(20000001)
  net <- blockwiseModules(t(dataset), power = bestpower, networkType = "signed hybrid",
                          TOMType = "signed", minModuleSize = minmodsize,
                          reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                          numericLabels = TRUE, pamStage = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, minKMEtoStay = 0.2, minCoreKME = 0.3, minCoreKMESize = 5,
                          saveTOMFileBase = paste0(file_out), useBranchEigennodeDissim = FALSE,
                          verbose = 1, maxBlockSize = 25000, corType = "bicor", quickCor = 1,
                          nThreads = 0)
  
  mergedColors = labels2colors(net$colors)
  
  
  rownames(net$MEs) <- colnames(dataset)
  ME <- net$MEs
  ME.colors <- labels2colors(as.numeric(gsub("ME", "", colnames(ME))))
  METree <- hclust(dist(t(ME)), method = "average");
  
  net$METree <- METree
  net$ME.colors <- ME.colors
  net$sft <- sft
  clusteringResult <- data.frame(cluster = as.numeric(net$colors), gene = rownames(dataset))
  net$clusteringResult <- clusteringResult
  net$analysis <- analysis
  net$dataset <- dataset
  net$powers  <- powers
  net$bestpower <- bestpower
  saveRDS(net,paste0(file_out, organ, "_", analysis, "_net.RDS"))
}

#dataset is the gene expression data
#analysis is an indication
#organ is liver
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
clusterWGCNA(dataset = dataset, powers = powers, analysis = "all", organ = "liver")

########################################################################
#Figure S3
########################################################################

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

########################################################################
#check the MHS-associated modules
########################################################################

MEs2       <- as.data.frame(net$MEs)
#module 0 is not a meaningful module
res_lt <- do.call(rbind, mclapply(colnames(MEs2)[colnames(MEs2) %ni% c("id", "ME0")], function(probe){ #probe= colnames(expression)[2]
  expression_tmp <- as.data.frame(MEs2[, probe])
  rownames(expression_tmp) <- rownames(MEs2)
  colnames(expression_tmp) <- "module"
  expression_tmp$sample    <- rownames(expression_tmp)
    phenotype_tmp           <- as.data.frame(pheno[,"MHS"])
    rownames(phenotype_tmp) <- rownames(pheno)
    colnames(phenotype_tmp) <- "phenotype"
    phenotype_tmp$sample    <- rownames(phenotype_tmp)
    tmp                     <- merge(expression_tmp, phenotype_tmp, by = "sample")
    tmp$diet                <- gsub(".*_", "", tmp$sample)
    tmp$strain              <- gsub("_.*", "", tmp$sample)
    tmp                     <- tmp[!is.na(tmp$phenotype), ]
    set.seed(12072021)
    lm1                     <- lm(formula =  phenotype ~ module + diet, data = tmp)
    coe                     <- summary(lm1)[["coefficients"]]
    cor                     <- data.table(
      slope   = coe[2,1],
      p.slope = coe[2,4],
      module  = probe
    )
}, mc.cores = 20))

res_lt$adjp.slope <- p.adjust(res_lt$p.slope, method = "BH")


########################################################################
#get the enriched genesets of the MHS-related modules
########################################################################

res             <- clusterProfiler::enricher(gene = genes.names.tmp,
                                             TERM2GENE = Genesets[[geneset]][, c("gs_name", "gene_symbol")],
                                             universe = universe,
                                             minGSSize = 10,
                                             maxGSSize = 1500)



################Correlation analyses##################################

data_cor <- do.call(rbind, lapply(colnames(lipid_ave)[colnames(lipid_ave) %ni% c("id")], function(lipid){ #lipid = colnames(lipid_ave)[colnames(lipid_ave) %ni% c("id")][1]
  #print(lipid)
  phenotype_tmp           <- as.data.frame(lipid_ave[,c(lipid, "id")])
  colnames(phenotype_tmp) <- c("lipid", "id")
  tmp                     <- merge(pheno_data, phenotype_tmp, by = "id")
  tmp                     <- tmp[!is.na(tmp$MHS), ]
  tmp                     <- tmp[!is.na(tmp$lipid), ]
  cor_test                <- cor.test(tmp$MHS, tmp$lipid, method = "pearson")
  cor                     <- data.table(r = cor_test$estimate, p = cor_test$p.value, pheno = "MHS", lipid = lipid)
  return(cor)
}))




