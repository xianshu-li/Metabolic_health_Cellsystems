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
#GWAS analyses for 43 inbred strains
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

####phewad.clac######
phewas.calc <- function(pheno, geno, kmatrix, pmatrix=NULL){
  pheno.num <- ncol(pheno)
  for(i in 1:pheno.num){ #i=1
    print(i)
    pheno.i <- pheno[, i]
    names(pheno.i) <- rownames(pheno)
    pheno.i.not_na <- !is.na(pheno.i)
    pheno.i <- pheno.i[pheno.i.not_na]
    pheno.length <- length(pheno.i)
    geno.at_index  <- t(geno)
    geno.at_index <- geno.at_index[rownames(geno.at_index)%in%names(pheno.i),]
    if(pheno.length < 15){
      pheno.name <- colnames(pheno)[i]
      message(paste0("phenotype ", pheno.name, " has ", pheno.length, " strains, not enough for analysis!"))
      next
    }else{
      pheno.name <- colnames(pheno)[i]
      geno_imp.at_index <- as.matrix(geno.at_index)
      kmatrix.i <- kmatrix[rownames(kmatrix)%in%names(pheno.i), colnames(kmatrix)%in%names(pheno.i)]
      names(pheno.i)
      mygwas <- mlmm(Y = pheno.i, X = geno_imp.at_index, K = kmatrix.i, nbchunks = 2, maxsteps = 3)
      pval <- mygwas$pval_step[[1]]$out
      colnames(pval)[2] <- pheno.name
      if (length(pmatrix) == 0){
        pmatrix <- pval
      }else{
        pmatrix <- join(pmatrix, pval, by = "SNP", type = "left", match = "all")
      }
    }
  }
  return(pmatrix)
}


########################################################################
# function to calculate kinship matrix from genotype data
# input:
# geno: a p by m numeric genotype matrix containing p markers in rows and m samples in columns
#     colnames(genotype) should be the sample IDs
#     NOTE: the input genotype matrix should be coded in a minor allele fashion,
#     i.e. minor allel = 0, major allele = 1, heterozygous = 0.5
########################################################################
kinship.emma <- function(geno) { #geno = geno_ori[, 1:2]
  K_mat <- emma.kinship(geno)            # Calculate kinship using the emma package
  colnames(K_mat) <- colnames(geno)      # emma does not set the rownames and column names, so we should do it
  rownames(K_mat) <- colnames(geno)
  return(K_mat)
}


