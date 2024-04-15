########################################################################
#clean.workplace
########################################################################

rm(list = ls())

########################################################################
#library
########################################################################

library(qtl2)
library(ggplot2)
require(plyr)
library(dplyr)
library(data.table)
library(ggrepel)
library(qtl2convert)
library(readxl)
require(emma)
require(mlmm)
library(MASS)
library(parallel)
library(ggpubr)
library(effsize)
library(RColorBrewer)
library(scales)

########################################################################
#QTL mapping
#same code for BXD and CC mice
########################################################################

# Operator for "Not in"
`%ni%` = Negate(`%in%`)

#pheno.norm#
pheno.norm <- function(pheno.data){ # pheno.data <- data_tmp[, 20]

  # test if there is any zero or negative values in the phenotype data
  if(length(which(pheno.data <= 0)) != 0){
    # if non-positive value exists, use quantile transformation
    message(paste0("Phenotype has non-positive values, quantile transformation will be applied!"))
    quantNorm = function(x){
      qnorm(rank(x, na.last = "keep",ties.method = "average")/(length(x)+1))
    }
    pheno.data.norm <- quantNorm(pheno.data)
  }else{
    # if all values are positive, boxcox transformation could be used
    #message(paste0("Phenotype has only positive values, Boxcox transformation will be applied!"))
    # scale all data into data that centers around 1, by dividing the average
    if (length(which(pheno.data >= 0)) != 0) {
      pheno.center <- pheno.data / mean(pheno.data, na.rm = TRUE)
      # run the box-cox transformation
      bc <- boxcox(pheno.center ~ 1, lambda = seq(-2, 2, 0.25), plotit=FALSE)
      trans <- bc$x[which.max(bc$y)]
      if(trans == 0){
        pheno.data.norm <- log(pheno.center)
      }else{
        pheno.data.norm <- (pheno.center^trans - 1)/trans
      }
    }else{
      pheno.data.norm <- NA
    }
  }
  return(pheno.data.norm)
}


########################################################################
#Data input and output paths
########################################################################


if (!dir.exists(file_out)) {
  dir.create(file_out, recursive = T)
}


#obtained from genenetwork
geno
colnames(geno)               <- gsub("C57BL.6J", "C57BL/6J", gsub("DBA.2J", "DBA/2J", colnames(geno)))
#because the order of these two locus in the physical and genetic map does not match, remove them
geno                         <- geno[geno$Locus %ni% c("UNC5348732", "rs6377403"),]

phenotypic_combine           <- openxlsx::read.xlsx("./Data_upload/Phenotypic data.xlsx")
rownames(phenotypic_combine) <- phenotypic_combine$Strain_id

Diets                        <- c("CD","HFD")

lapply(Diets, function(diet){ #diet = Diets[1]
  print(diet)
  
  if (!dir.exists(paste0(file_out, diet, "/rawdata"))) {
    dir.create(paste0(file_out, diet, "/rawdata"), recursive = T)
  }
  
  print("read phenotypic data")
  
  pheno                      <- phenotypic_combine[phenotypic_combine$diet == diet,]
  pheno                      <- pheno[pheno$strain %in% colnames(geno),]
  
  #Perform QTL mapping for the MHS and its components
  phenotype_names            <- colnames(phenotypic_combine)[colnames(phenotypic_combine) %ni% c("Strain_id", "strain", "diet", "id")]
  
  tmp_res                    <- do.call(rbind, lapply(phenotype_names, function(phenotype_name){ #phenotype_name = phenotype_names[9]
    
    print(phenotype_name)
    pheno_tmp                <- pheno[, c("Strain_id", phenotype_name)]
    pheno_tmp                <- pheno_tmp[!is.na(pheno_tmp[,2]), ]
    pheno_norm               <- as.data.frame(pheno.norm(pheno_tmp[,2]))
    rownames(pheno_norm)     <- pheno_tmp$Strain_id
    colnames(pheno_norm)     <- phenotype_name
  
  write2csv(pheno_norm, paste0(file_out, diet, "/rawdata/bxd_pheno_", phenotype_name, ".csv"), overwrite=TRUE, row.names = "Strain",
            "BXD phenotype data")
  
  # read the genotypes (and maps)
  geno_id <- lapply(pheno_tmp$Strain_id, function(x){ # x <- pheno$id[1]
    #print(x)
    y <- unique(pheno[pheno$Strain_id == x,]$strain)
    z <- data.frame(geno[, which(colnames(geno) == y)])
    colnames(z) <- x
    return(z)
  })
  geno_id <- do.call(cbind, geno_id)
  g <- cbind(geno[, c("Chr", "Locus", "Mb_mm9", "Mb_mm10", "cM_BXD")], geno_id)
  #g <- cbind(geno[, c("Chr", "Locus", "cM", "Mb")], geno_id)
  g <- g[-which(g$Chr=="Y"),]
  g <- g[-which(g$Chr=="M"),]
  
  # unique(g$Chr)
  # grab genetic and physical maps
  gmap <- data.frame(g[,c("Locus","Chr","cM_BXD")])
  colnames(gmap) <- c("marker", "chr", "pos")
  rownames(gmap) <- gmap[,"marker"]
  
  pmap <- data.frame(g[,c("Locus","Chr","Mb_mm10")])
  colnames(pmap) <- c("marker", "chr", "pos")
  rownames(pmap) <- pmap[,"marker"]
  
  # reorder genetic and physical map by physical order
  pmap <- pmap[order(factor(pmap$chr, c(1:19,"X")), as.numeric(pmap$pos)),]
  gmap <- gmap[rownames(pmap),]
  # write genetic and physical maps
  write2csv(gmap, paste0(file_out, diet, "/rawdata/bxd_gmap_", phenotype_name, ".csv"), comment="Genetic map for BXD data", overwrite=TRUE)
  write2csv(pmap, paste0(file_out, diet, "/rawdata/bxd_pmap_", phenotype_name, ".csv"), comment="Physical map (mm10 Mbp) for BXD data", overwrite=TRUE)
  # grab genotypes
  g <- g[, !(colnames(g) %in% c("Chr", "Mb_mm9", "Mb_mm10", "cM_BXD"))]
  colnames(g)[colnames(g)=="Locus"] <- "marker"
  rownames(g) <- g[,"marker"]
  # reorder markers as in the maps
  g <- g[rownames(pmap),]
  # write the genotypes
  write2csv(g, paste0(file_out, diet, "/rawdata/bxd_geno_", phenotype_name, ".csv"), comment="Genotypes for BXD data", overwrite=TRUE)
  # write cross info (BxD for all strains)
  crossinfo <- data.frame(id=names(g)[-1], cross_direction=rep("BxD", ncol(g)-1))
  write2csv(crossinfo, paste0(file_out, diet, "/rawdata/bxd_crossinfo_", phenotype_name, ".csv"), overwrite=TRUE,
            comment=paste0("Cross info for BXD data\n#",
                           "(all lines formed from cross between female B and male D)"))
  
  # write json file
  write_control_file(paste0(file_out, diet,"/rawdata/bxd_", phenotype_name, ".json"),
                     description="BXD mouse data from the fasted study",
                     crosstype="risib",
                     geno_file=paste0("bxd_geno_", phenotype_name, ".csv"),
                     geno_transposed=TRUE,
                     geno_codes=list(B=1, D=2),
                     xchr="X",
                     gmap_file=paste0("bxd_gmap_", phenotype_name, ".csv"),
                     pmap_file=paste0("bxd_pmap_", phenotype_name, ".csv"),
                     pheno_file=paste0("bxd_pheno_", phenotype_name, ".csv"),
                     crossinfo_file = paste0("bxd_crossinfo_", phenotype_name, ".csv"),
                     crossinfo_codes = c("BxD"=0),
                     alleles=c("B", "D"),
                     na.strings=c("-", "H","NA"),
                     overwrite=TRUE)
  ###########read.data################
  print("read.data")
  BXD <- read_cross2(paste0(file_out, diet,"/rawdata/bxd_", phenotype_name, ".json"))
  summary(BXD)
  ###########Calculating genotype probabilities#######               
  map <- insert_pseudomarkers(BXD$gmap, step=1)
  ############calculating the QTL genotype probabilities
  pr <- calc_genoprob(BXD, map, error_prob=0.002)
  apr <- genoprob_to_alleleprob(pr)
  ###########Calculating a kinship matrix##########
  print("Calculating a kinship matrix")
  grid <- calc_grid(map = map, step=1)
  pr_grid <- probs_to_grid(pr, grid)
  kinship_grid <- calc_kinship(pr_grid, "loco")
  ##########Performing a genome scan########
  print("Performing a genome scan")
  set.seed(20210609)
  out_kinship <- scan1(pr_grid, BXD$pheno, kinship_grid, cores = 45)
  dir.create(paste0(file_out, diet,"/result"))
  save(out_kinship, file = paste0(file_out, diet,"/result/out_kinship_", phenotype_name, ".RData"))
  ##############Performing a permutation test########
  print("Performing a permutation test")
  set.seed(20210609)
  operm <- scan1perm(pr_grid, BXD$pheno, kinship_grid, n_perm=10000, cores = 35)
  save(operm, file = paste0(file_out, diet,"/result/permutation_test_", phenotype_name, ".RData"))
  # load(paste0(file_out, diet,"/result/permutation_test.RData"))
  # load(paste0(file_out, diet,"/result/out_kinship.RData"))
  #use P value = 0.05 and P value = 0.1 to find QTL peaks
  res3        <- find_peaks(out_kinship, map, threshold= summary(operm, 0.1), prob=0.95, peakdrop=5)
  
  res3
  
  }))
  write.table(tmp_res, paste0(file_out, diet,"/result/significant_qtl2_gene_position_0.1.txt"), quote=F, sep='\t', row.names=F)
  
})


########################################################################
#Manhattan plot (Fig 4A)
########################################################################

theme_graphs <- theme_bw(base_size = 12) + 
  theme(axis.text=element_text(size=8,color="black"),
        plot.title=element_text(face = "bold",size=7,hjust = 0.5,vjust = 0.5),
        plot.subtitle = element_text(size=7,hjust = 0.5, vjust = 0.5),
        axis.title=element_text(size=12,hjust = 0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        strip.text = element_text(size=10),
        strip.background = element_blank(),
        legend.position = "right",
        #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
        element_line(color = "gray25", size = 0.5),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        panel.spacing = unit(0, "lines"))


#get the physical and genetic position of each marker
#diet = "HFD"
BXD                     <- read_cross2(paste0(file_out, diet,"/rawdata/bxd_MHS.json"))
map                     <- insert_pseudomarkers(BXD$gmap, step=1)
map_mb                  <- interp_map(map, BXD$gmap, BXD$pmap)

position                <- do.call(rbind, lapply(c(1:19, "X"), function(x){ #x= "X"
  position_cm           <- as.data.frame(map[[x]])
  colnames(position_cm) <- "position_cm"
  position_cm$locus     <- rownames(position_cm)
  position_mb           <- as.data.frame(map_mb[[x]])
  colnames(position_mb) <- "position_mb"
  position_mb$locus     <- rownames(position_mb)
  position_mer          <- merge(position_mb, position_cm, by = "locus")
  position_mer$chr      <- x
  position_mer
}))


#geno preparation
position$chr            <- as.numeric(gsub("X", "20", position$chr))
Geno                    <- as.data.frame(position[,c("chr","locus","position_mb")])
colnames(Geno)          <- c("Chr","Locus","Mb_mm10")
Geno$BP                 <- Geno$Mb_mm10*10^6


Diets_lt                <- lapply(Diets, function(diet){ #diet = Diets[1]
      print(diet)
      #Find gene region of the peak
      print("read out_kinship")
      list_result_file      <- list.files(path = paste0(file_out, "/", diet,"/result"), full.names = T)
      list_out_kinship      <- list_result_file[grepl("_kinship_MHS",list_result_file)]
      load(list_out_kinship, verbose = T)
      out_kinship           <- as.data.frame(out_kinship)
      out_kinship$locus     <- rownames(out_kinship)
      out_kinship           <- reshape2::melt(out_kinship[, c("locus", "MHS")], id.vars = "locus")
      out_kinship$diets     <- diet

     colnames(out_kinship)  <- c("Locus", "phenotypes", "LOD", "diets")
     data_1                 <- merge(out_kinship, Geno, by = "Locus")
     data_1                 <- data_1[,c("Locus", "phenotypes", "LOD", "Chr","BP")]
     colnames(data_1)       <-c("SNP", "phenotype", "LOD", "CHR", "BP")
     data_1$CHR              <- factor(data_1$CHR, levels = c(1:20), labels = c(1:19, "X"))

     # Compute chromosome size
     a <- data_1 %>% group_by(CHR) %>% 
     summarise(chr_len=max(BP)) %>%
     # Calculate cumulative position of each chromosome
     mutate(tot=as.numeric(cumsum(chr_len))- as.numeric(chr_len)) %>%
     #select(-chr_len) %>%
     # Add this info to the initial dataset
     left_join(data_1, ., by=c("CHR"="CHR")) %>%
     # Add a cumulative position of each SNP
     arrange(CHR, BP) %>%
     mutate(BPcum=as.numeric(BP)+as.numeric(tot))

     axisdf = a %>% group_by(CHR) %>% dplyr::summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

     if(diet == "HFD"){
        color = "#1b9e77"
    }else{
        color = "#d95f02"
    }
     
     load(paste0(file_out, diet,"/result/permutation_test_MHS.RData"), verbose = T)
     thres    <- summary(operm, 0.05)[1]
     thres_s  <- summary(operm, 0.1)[1]
  p1          <- ggplot(a, aes(x=BPcum, y=LOD)) +
    geom_line(size = 0.2, color = color) +
    geom_hline(yintercept = thres, color = color, size = 0.2, linetype = c("solid")) +
    #geom_point(data=subset(x, signif=="Yes"), color = pals::brewer.dark2(3)[3], shape = 17, size =2, alpha= 0.8) +
    geom_hline(yintercept = thres_s, color = color, size = 0.2, linetype = c("twodash")) +
    scale_color_manual("",values = color) +
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0, 0.2)) + 
    theme_graphs + 
    xlab("Chromosomes") + ylab("LOD score")
  
  ggsave(filename = paste0("./Figures/11-manhattan_",diet, ".png"), p1, width = 5.5, height = 2.5)
  ggsave(filename = paste0("./Figures/11-manhattan_",diet, ".pdf"), p1, width = 5.5, height = 2.5, useDingbats = FALSE)
  
  return(p1)
})

########################################################################
#Boxplot plot (Fig 4B, 4C)
########################################################################
#get locus with the highest LOD under the significant QTL peak
locus                <- c("rs32906553", "rs31874398")
names(locus)         <- c("CD_Chr8", "HFD_Chr7")

mapply(locu = locus, name = names(locus), function(locu, name){
  print(name)
  #name = names(locus)[2]
  #locu = locus[2]
  pheno            <- phenotypic_combine
  geno             <- read.csv(paste0(file_out, gsub("_.*", "", name),"/rawdata/bxd_geno.csv"),comment.char = "#")
  geno.tmp         <- geno[geno$marker == locu, ]
  geno.tmp         <- reshape2::melt(geno.tmp, id.vars = "marker")
  geno.tmp$Strain  <- gsub("^X", "", gsub("\\.", "-", geno.tmp$variable))
  pheno.tmp        <- pheno[, c("Strain_id", "MHS")]
  
  tmp              <- merge(pheno.tmp, geno.tmp, by.x = "Strain_id", by.y = "Strain")
  tmp              <- tmp[!is.na(tmp$MHS), ]
  
  p <- ggplot(tmp, aes(x=value, y= MHS, color= value, fill = value)) +
    geom_boxplot(alpha = 0.5, outlier.color = "white") +
    geom_jitter(aes(colour = value)) +
    stat_boxplot(geom= 'errorbar' , width = 0.3, position = position_dodge(width = 0.75) ) +
    stat_compare_means(comparisons = list(c("B", "D")), label = "p.signif",position = "identity",method = "t.test") +
    theme_graphs + theme(panel.grid =element_blank(),legend.position = "right") +
    scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
    scale_colour_manual(values = c("B" = "#CB181D", "D" = "#2171B5")) +
    scale_fill_manual(values = c("B" = "#CB181D", "D" = "#2171B5")) +
    labs(x="",
         y="") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
  
})

########################################################################
#Dot plot (Fig 4D) 
#Obtained QTL peaks of the components of MHS
#then plot the QTL peaks
########################################################################

phenotype_name             <- data.frame(phenotype = colnames(phenotypic_combine)[colnames(phenotypic_combine) %ni% c("Strain_id", "strain", "diet")],
                                         name = c("MHS", "BW (gained)", "Fat (%)", "Glucose (AUC)", "Insulin (AUC)", "Glucose (OGTT)", "Insulin (OGTT)", "HOMA (IR)", "eWAT (%)",
                                                  "Liver (%)", "pWAT (%)", "scWAT (%)", "ALAT", "Glucose", "HDL", "LDL", "Total Cholesterol", "Triglycerides"))

lapply(Diets, function(diet){
  
list_result_file         <- list.files(path = paste0(file_out, diet,"/result"), full.names = T)
list_out_kinship         <- list_result_file[grepl("_kinship",list_result_file)]

res  <- do.call(cbind, lapply(list_out_kinship, function(i){ # i= list_out_kinship[1]
  
  load(i, verbose = T)
  
  return(out_kinship)
  
}))


res                       <- as.data.frame(res[, phenotype_name$phenotype])

res$Locus                  <- rownames(res)

data_melt                  <- reshape2::melt(res, id.vars = "Locus")
colnames(data_melt)        <- c("Locus", "phenotypes", "LOD")

data_melt_table            <- as.data.table(data_melt)
data_1                     <- merge(data_melt_table, position, by.x = "Locus", by.y = "locus")
colnames(data_1)           <- c("SNP", "phenotype", "LOD","BP","CM", "CHR")
data_1$CHR                 <- factor(data_1$CHR, levels = c(1:20), labels = c(1:19, "X"))
data_1$BP                  <- data_1$BP * 10^6
# Compute chromosome size
a <- data_1 %>% group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=as.numeric(cumsum(chr_len))- as.numeric(chr_len)) %>%
  #select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(data_1, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=as.numeric(BP)+as.numeric(tot))
axisdf = a %>% group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

accm  <- unique(a[, c("CHR", "tot")])


#Find significant QTLs
print("read significant QTLs")
list_significant_QTLs    <- list_result_file[grepl("gene_position_0.1.txt",list_result_file)]
significant_QTLs         <- read.table(list_significant_QTLs, header = T)

QTLs                     <- do.call(rbind, lapply(1:nrow(significant_QTLs), function(i){ #i =1
  x                          <- significant_QTLs[i,]
  
  start                      <- position[position$position_cm == x$ci_lo & position$chr == x$chr, "position_mb"][1] * 10^6
  end                        <- position[position$position_cm == x$ci_hi & position$chr == x$chr, "position_mb"][1] * 10^6
  pos_mb                     <- position[position$position_cm == x$pos & position$chr == x$chr, "position_mb"][1] * 10^6
  
  start_accm                 <- start + accm$tot[accm$CHR == x$chr ]
  end_accm                   <- end   + accm$tot[accm$CHR == x$chr ]
  pos_accm                   <- pos_mb + accm$tot[accm$CHR == x$chr ]
  tmp                        <- data.frame(phenotype = x$lodcolumn, chr = x$chr, pos = x$pos, pos_mb = pos_mb, LOD = x$lod, 
                                           start_cm = x$ci_lo, end_cm = x$ci_hi, start_mb = start, end_mb = end, 
                                           start_accm = start_accm, end_accm = end_accm, pos_accm = pos_accm)
  tmp
}))


QTLs                    <- merge(QTLs, phenotype_name, by = "phenotype")

if(diet == "HFD"){
  color = "#1b9e77"
  #QTLs$name            <- factor(QTLs$name, levels = rev(c("MHS", "Glucose", "HOMA (B)", "BW (gained)", "Fat (%)", "Lean (%)", "Blood pressure (Systolic)","Liver (%)", "eWAT (%)","Gastrocnemius (%)","Kidney (%)", "ALPL")))
  
}else{
  color = "#d95f02"
  #QTLs$name            <- factor(QTLs$name, levels = rev(c("MHS", "Glucose", "Glucose (AUC)", "Insulin (AUC)", "BW (gained)", "Lean (%)", "Heart (%)", "Blood pressure (Diastolic)", "ASAT", "ALAT", "ALPL", "RER (Night)")))
}


p1 <- ggplot(a, aes(x=BPcum)) +
  geom_segment(data = QTLs, aes(x=start_accm, xend= end_accm, y=name, yend=name), color="black") +
  geom_point(data = QTLs, aes(x=pos_accm, y=name), color=color, size = 3) +
  geom_vline(xintercept = accm$tot[-1], size = 0.2, linetype = c("twodash")) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.01, 0.01), limits= c(0, max(a$BPcum))) +
  theme_bw(base_size = 14) + 
  theme(axis.text=element_text(size=8,color="black"),
        plot.title=element_text(size=8,hjust = 0.5,vjust = 0.5),
        plot.subtitle = element_text(size=8,hjust = 0.5, vjust = 0.5),
        axis.title=element_text(size=8,hjust = 0.5),
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        legend.position = "none",
        strip.text = element_text(size=8),
        strip.background = element_blank(),
        #strip.background = element_rect(fill = "gray88" ,colour="black", size = 0.3),
        element_line(color = "gray25", size = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  xlab("Chromosomes") + ylab("LOD score")

})

########################################################################
#Dot plot (Fig 4E) 
#Effect size
########################################################################

Diets_lt <- do.call(rbind, lapply(Diets, function(diet){ #diet = Diets[1]
  print(diet)
  print("read phenotypic data")
  
  phenotypic_data       <- phenotypic_combine[phenotypic_combine$diet == diet, ]
  phenotypic_data       <- phenotypic_data[phenotypic_data$strain %in% colnames(geno), ]
  phenotypic_data
  # read the genotypes (and maps)
  geno_id <- lapply(phenotypic_data$Strain_id, function(x){ #
    #print(x)
    y <- phenotypic_data[phenotypic_data$Strain_id==x,]$strain
    z <- data.frame(geno[, which(colnames(geno) == y)])
    colnames(z) <- x
    return(z)
  })
  geno_id <- do.call(cbind, geno_id)
  g <- cbind(geno[, c("Chr", "Locus", "Mb_mm9", "Mb_mm10", "cM_BXD")], geno_id)
  g <- g[!grepl("Y|M", g$Chr),]



  pheno_selected     <- data.table(

    phenotype = c( "scWAT.", "pWAT.", "body_fat", "Blood_Glucose_.mmol.L.","Resting_insulin", "OGTT_AUC_Insulin_.AUC.", "HOMA_IR","Liver.", "Blood_ALAT_.u.L.", "Blood_TotalCholesterol_.mmol.L."),
    labels    = c( "scWAT (%)", "pWAT (%)","Body fat (%)", "Glucose", "Insulin","insulin (AUC)", "HOMA (IR)","Liver weight (%)", "ALT", "Total cholesterol"),
    category  = c( "Obesity","Obesity","Obesity", "Insulin resistance", "Insulin resistance", "Insulin resistance", "Insulin resistance","Liver damage", "Liver damage", "Liver damage")

  )
  
  result       <- do.call(rbind, lapply(locus, function(x){ #x = locus[1]
    geno_data  <- as.data.frame(g[g$Locus == x,])
    geno_data  <- geno_data[, colnames(geno_data)[!grepl("Chr|Mb_mm9|Mb_mm10|cM_BXD", colnames(geno_data))]]
    geno_data  <- reshape2::melt(geno_data, id.vars = "Locus")
    colnames(geno_data) <- c("Locus", "id", "genotype")
    
    pheno.data  <- phenotypic_data[, pheno_selected$phenotype]
    
    pheno.data_norm <- do.call(cbind, lapply(pheno_selected$phenotype, function(i){ # i = phenotype_name$phenotype[1]
      
      
      tmp       <- pheno.data[, i]
      names(tmp)<- rownames(pheno.data)
      tmp       <- tmp[!is.na(tmp)]
      tmp       <- pheno.norm(tmp)
      tmp_return<- as.data.frame(tmp[match(rownames(pheno.data), names(tmp))])
      colnames(tmp_return) <- i
      rownames(tmp_return) <- rownames(pheno.data)
      tmp_return
      
    }))
    
    pheno.data_norm$Strain_id <- rownames(pheno.data_norm)
    
    pheno.data_norm           <- reshape2::melt(pheno.data_norm, id.vars = "Strain_id")
    colnames(pheno.data_norm) <- c("id", "phenotype", "value")
    pheno.geno.data <- merge(pheno.data_norm, geno_data, by = "id")
    pheno.geno_lt   <- split(pheno.geno.data, pheno.geno.data$phenotype)
    
    result_2        <- do.call(rbind, lapply(pheno.geno_lt, function(y){ #y = pheno.geno_lt[[1]]
      
      print(as.character(unique(y$phenotype)))
      y          <- y[y$value < 1e80, ]
      y          <- y[!is.na(y$value), ]
      if (length(y$value) >20) {
        bb            <- y[y$genotype == "B", "value"]
        dd            <- y[y$genotype == "D", "value"]
        p             <- t.test(bb,dd)$p.value
        effect_size   <- cohen.d(bb[!is.na(bb)], dd[!is.na(dd)])$estimate
        effect_factor <- as.character(cohen.d(bb[!is.na(bb)], dd[!is.na(dd)])$magnitude)
        
        p             <- data.frame(phenotype = as.character(unique(y$phenotype)),
                                    locus     = x,
                                    P_value   = p,
                                    effect_size = effect_size,
                                    effect_factor = effect_factor,
                                    stringsAsFactors = F)
      }else{
        p <- NULL
      }
      return(p)
    }))
  }))
  
  result$adjustP <- p.adjust(result$P_value, method = "BH")
  result$log10P  <- -log10(result$adjustP)
  
  result$diet    <- diet
  return(result)
}))

result_data        <- Diets_lt[(Diets_lt$locus == "rs32906553" & Diets_lt$diet == "CD")|(Diets_lt$locus == "rs31874398" & Diets_lt$diet == "HFD"), ]
result_data        <- merge(pheno_selected, result_data, by = "phenotype")

result_data       <- result_data[order(result_data$effect_size, decreasing = T),]

# Hierarchical clustering on phenotypes1
data_dcast <- reshape2::dcast(result_data, labels~locus, value.var="effect_size")
data_dcast[is.na(data_dcast)] <- 0 # For the purpose of clustering, we consider "NA" equivalent to no correlation
rownames(data_dcast) <- data_dcast$labels
hc <- hclust(dist(t(data_dcast[,-1])))
rowInd <- hclust(dist(data_dcast[,-1]))$order
colInd <- hclust(dist(t(data_dcast[,-1])))$order

result_data$locus <- factor(result_data$locus, levels = rev(colnames(data_dcast)[-1][colInd]))
result_data$labels <- factor(result_data$labels, levels = data_dcast$labels[rowInd])

Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))

result_data$stars <- cut(result_data$adjustP, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

result_data$locus <- factor(result_data$locus, levels = c("rs32906553", "rs31874398"))

g <- ggplot(result_data, aes(x = locus, y = labels, fill = effect_size)) +
  geom_point(aes(size = log10P, fill = effect_size, color = effect_size)) +
  scale_size(range = c(3,8)) +
  scale_fill_gradientn(colours= Heatmap_palette , limits=c(-1,1) ,na.value="gray87", oob=squish) +
  scale_color_gradientn(colours= Heatmap_palette , limits=c(-1,1) ,na.value="gray87", oob=squish) +
  facet_grid(rows = vars(category), space="free", scales="free") +
  geom_text(aes(label=stars), color="black", size=4, vjust = 0.8, hjust = 0.5) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    #axis.title = element_text(size = 13),
    axis.text.x = element_text(angle= 90, hjust= 1, vjust=0.5,size = 11, color = "black"),
    #axis.text.x = element_text( hjust= 0.5, vjust=0.5,size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_blank(),
    #axis.ticks = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 11)
    #panel.border = element_blank()
  ) 



