library(tidyverse)
library(ggrepel)
library(patchwork)


###############################################################################################

# Figure 3F:

#calculate the relationship between the proportions of copy number change before (MCN = 1) and after WGD (MCN>=2) with the gene density (including essential genes, tumor suppersor genes and oncogenes ) in each arm.

###############################################################################################

#Define functions.

#' get.seg.mat.arm
#' 
#' get the arm-level CNV from segment files.
#' 
#' 
#' 
get.seg.mat.arm <- function (seg.mat.copy){
  ### centromere positions
  data("centromere", package = "EstimateClonality")
  seg.mat.copy$PQ <- NA
  
  seg.mat.copy.arm  <- c()
  for (chr in c(1:22)) {
    chr.mat.copy <- seg.mat.copy[seg.mat.copy$chr==chr,]
    
    #label those that finish before and after the centromere with obvious labeling
    chrPs <- chr.mat.copy$endpos<centromere[centromere[,1]==chr,3]
    if(TRUE%in%chrPs)
    {
      chr.mat.copy[chr.mat.copy$endpos<centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"p",sep="")
      
    }
    
    chrQs <- chr.mat.copy$startpos>centromere[centromere[,1]==chr,3]
    if(TRUE%in%chrQs)
    {
      chr.mat.copy[chr.mat.copy$startpos>centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"q",sep="")
    }
    
    chr.mat.copy.done <- chr.mat.copy[!is.na(chr.mat.copy$PQ),,drop=FALSE]
    # work out how many cases do we need to do more
    chr.mat.missing <- chr.mat.copy[is.na(chr.mat.copy$PQ),,drop=FALSE]
    
    if(nrow(chr.mat.missing)==0)
    {
      seg.mat.copy.arm <- rbind(seg.mat.copy.arm,chr.mat.copy.done)
      
    }
    
    if(nrow(chr.mat.missing)>=1)
    {
      chr.mat.missing.p <- chr.mat.missing.q <-chr.mat.missing 
      
      chr.mat.missing.p$endpos <- as.numeric(centromere[centromere[,1]==chr,3])-1
      chr.mat.missing.q$startpos <- as.numeric(centromere[centromere[,1]==chr,3])+1
      
      chr.mat.missing.p[chr.mat.missing.p$endpos<centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"p",sep="")
      chr.mat.missing.q[chr.mat.missing.q$startpos>centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"q",sep="")
      
      chr.mat.copy.compl <- rbind(chr.mat.copy.done,chr.mat.missing.p,chr.mat.missing.q)
      
      chr.mat.copy.compl <- chr.mat.copy.compl[order(chr.mat.copy.compl$sample,chr.mat.copy.compl$startpos),,drop=FALSE]
      
      seg.mat.copy.arm <- rbind(seg.mat.copy.arm,chr.mat.copy.compl)
      
    }
    
  }
  
  return(seg.mat.copy.arm)
}



#" Plot Scores.
funPlotDavoliScore <- function(DavoliScore,
                               GenomicMeasure,
                               cexSize=1.5,
                               y_lab="Essential Score",
                               x_lab="Haploid Mean",
                               chrNames= NULL,
                               label.x.cor.pos = 0.9, label.y.cor.pos = 0.9,
                               title = NULL){
  
  df = data.frame(
    x = as.numeric(GenomicMeasure),
    y = as.numeric(DavoliScore),
    label = chrNames
  )
  
  annotate_text = df %>%
    summarise(p.value = (cor.test(x, y)$p.value ),
              rho = (cor.test(x, y, method = "spearman")$estimate ),
              x.pos = min(x)+ ((max(x)-min(x)) *label.x.cor.pos ),
              y.pos = min(y)+ ((max(y)-min(y))*label.y.cor.pos )) %>%
    mutate(text = sprintf('R[sp]~"="~%.2f""~italic(P)~"="~%.1e', rho, p.value  ) ) %>%
    mutate(text1 = c(paste('r=',signif(rho,3),'\np=',signif(p.value,2)))  )
  
  
  df %>%
    ggplot(aes(x = x, y = y ) ) + ggpubr::theme_pubr() +
    geom_point(pch=21, fill="#3b518be5", size=2)+
    geom_text_repel(aes(label = chrNames), size = 3 ) +
    labs(x = x_lab, y = y_lab, title = title) +
    geom_smooth(method = "glm", formula = "y~x", se = T, col='blue', linetype = 2, size = 1.1) +
    geom_text(aes(x=x.pos, y= y.pos, label = text1) ,  col = "#258E4F", hjust = 0, vjust = 0, data= annotate_text, parse = F ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  
}


#' getEssential
#' get the correlation between  

getEssential = function(seg.mat.copy, cancer = "cancer"){
  
  seg.mat.copy = seg.mat.copy %>%
    rename(
      chr = Chr,
      startpos = Start,
      endpos = End,
      sample = SampleID,
      nAraw = nA,
      nBraw = nB
    ) %>%
    mutate(
      startpos = as.numeric(startpos),
      endpos = as.numeric(endpos)
    )
  
  seg.mat.copy.arm = get.seg.mat.arm(seg.mat.copy)
  
  tmp.seg <- get.seg.mat.arm(seg.mat.copy = seg.mat.copy)
  
  if( sum(tmp.seg$PQ=='21p') > 0){
    tmp.seg[tmp.seg$PQ=='21p',]$PQ <- '21q'
  }
  
  DataFolder = "geneSets/"
  
  #Get the essential scores for each arms.
  davoli <-   gdata::read.xls(paste( "resources/", DataFolder,"1-s2.0-S0092867413012877-mmc6.xlsx",sep=""),skip=1,stringsAsFactors=FALSE)
  davoli <- davoli[1:39,]
  davoli$Arm<-as.character(davoli$Arm)
  
  unique(tmp.seg$PQ)
  unique(davoli$Arm)
  
  
  # now let's count haploid LOH for each arm
  chrArmHaploidMedian <- c()
  chrArmHaploidMean <- c()
  chrArmHaploidSum <- c()
  chrArmHaploidSegs  <- c()
  #chrArmPostGDLoss   <- c()
  chrArmpostWGDMedian <- c()
  chrArmpostWGDMean <- c()
  chrArmpostWGDSum <- c()
  chrArmpostWGDSegs <- c()

  chrSize <- c()
  for (chrArm in unique(tmp.seg$PQ)){
    #message(chrArm)
    chrArmSeg  <- tmp.seg[tmp.seg$PQ==chrArm,]
    chrArmSize <- max(as.numeric(chrArmSeg$endpos))-min(as.numeric(chrArmSeg$startpos))
    chrSize<- c(chrSize,chrArmSize)
    haploidLOH_samples <- rep(0,length(unique(chrArmSeg$sample)))
    postGD_samples <- rep(0,length(unique(chrArmSeg$sample)))
    
    names(haploidLOH_samples) <- unique(chrArmSeg$sample)
    names(postGD_samples)     <- names(haploidLOH_samples)
    chrArmSegHap <- chrArmSeg[chrArmSeg$nAraw<=1&chrArmSeg$nBraw<=0.75,,drop=FALSE]
    chrArmSegWGDloss <- chrArmSeg[chrArmSeg$nBraw<1.5&chrArmSeg$nBraw>=0.75&chrArmSeg$nAraw>=1.5&(chrArmSeg$nAraw+chrArmSeg$nBraw)<=chrArmSeg$Ploidy,,drop=FALSE]
    
    for (sample in unique(chrArmSegHap$sample)) {
      specSampleHap <- chrArmSegHap[chrArmSegHap$sample%in%sample,,drop=FALSE]
      sum(as.numeric(specSampleHap$endpos)- as.numeric(specSampleHap$startpos))/chrArmSize
      
      haploidLOH_samples[sample] <- sum(as.numeric(specSampleHap$endpos)-as.numeric(specSampleHap$startpos))/chrArmSize
      
    }
    
    for (sample in unique(chrArmSegWGDloss$sample)) {
      specSampleWGDLoss <- chrArmSegWGDloss[chrArmSegWGDloss$sample%in%sample,,drop=FALSE]
      sum(specSampleWGDLoss$endpos-specSampleWGDLoss$startpos)/chrArmSize
      
      postGD_samples[sample] <- sum(as.numeric(specSampleWGDLoss$endpos)-as.numeric(specSampleWGDLoss$startpos))/chrArmSize
      
    }
    
    
    chrArmHaploidMedian <- c(chrArmHaploidMedian,median(haploidLOH_samples))
    chrArmHaploidMean <- c(chrArmHaploidMean,mean(haploidLOH_samples))
    chrArmHaploidSum <- c(chrArmHaploidSum,sum(haploidLOH_samples))
    chrArmHaploidSegs <- c(chrArmHaploidSegs,length(haploidLOH_samples[which(haploidLOH_samples>0)]))
    
    chrArmpostWGDMedian <- c(chrArmpostWGDMedian,median(postGD_samples))
    chrArmpostWGDMean <- c(chrArmpostWGDMean,mean(postGD_samples))
    chrArmpostWGDSum <- c(chrArmpostWGDSum,sum(postGD_samples))
    chrArmpostWGDSegs <- c(chrArmpostWGDSegs,length(postGD_samples[which(postGD_samples>0)]))
    
  }
  
  names(chrArmHaploidMean) <- unique(tmp.seg$PQ)
  names(chrArmHaploidSum) <- unique(tmp.seg$PQ)
  names(chrSize) <-  unique(tmp.seg$PQ)
  
  DavoliScore <- davoli$CharmEssential_score
  GenomicMeasure <- chrArmpostWGDMean
  
  davoli = davoli %>%
    mutate( cancer = cancer, 
      chrArmHaploidMedian = chrArmHaploidMedian,
           chrArmHaploidMean = chrArmHaploidMean,
           chrArmHaploidSum = chrArmHaploidSum,
           chrArmHaploidSegs  = chrArmHaploidSegs,
           chrArmpostWGDMedian = chrArmpostWGDMedian,
           chrArmpostWGDMean = chrArmpostWGDMean,
           chrArmpostWGDSum = chrArmpostWGDSum,
           chrArmpostWGDSegs = chrArmpostWGDSegs)
  
  p1 = funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score
                     ,GenomicMeasure = davoli$chrArmHaploidMean
                     ,y_lab="Essential Score"
                     ,x_lab = "Mean proportion haploid",
                     label.x.cor.pos = 0.5,
                     chrNames = davoli$Arm
                     )
  
  p2 = funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score
                     ,GenomicMeasure = davoli$chrArmpostWGDMean
                     ,x_lab = "Mean proportion post duplication loss",
                     label.x.cor.pos = 0.1,
                     chrNames = davoli$Arm
                     )
  
  p1 + p2
  
  list(
    davoli = davoli,
    plot = p1 +p2
  )
}

#########################################################################################################

plotSegs = function(davoli, title = NULL){
  
  plist = list()
  
  plist[[1]] = funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score,
                                  GenomicMeasure = davoli$chrArmHaploidMean,
                                  x_lab = "Chromosome arm proportion major allele <=1",
                                  y_lab = "Chromosome arm essential score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  
  plist[[2]] = funPlotDavoliScore(DavoliScore = davoli$CharmTSG_score,
                                  GenomicMeasure = davoli$chrArmHaploidMean,
                                  x_lab = "Chromosome arm proportion major allele <=1",
                                  y_lab = "Chromosome arm TSG score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  
  plist[[3]] = funPlotDavoliScore(DavoliScore = davoli$CharmOG_score,
                                  GenomicMeasure = davoli$chrArmHaploidMean,
                                  x_lab = "Chromosome arm proportion major allele <=1",
                                  y_lab = "Chromosome arm oncogenes score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  
  
  plist[[4]] = funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score,
                                  GenomicMeasure = davoli$chrArmpostWGDMean,
                                  x_lab = "Chromosome arm proportion loss post WGD",
                                  y_lab = "Chromosome arm essential score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  
  plist[[5]] = funPlotDavoliScore(DavoliScore = davoli$CharmTSG_score,
                                  GenomicMeasure = davoli$chrArmpostWGDMean,
                                  x_lab = "Chromosome arm proportion loss post WGD",
                                  y_lab = "Chromosome arm TSG score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  
  plist[[6]] = funPlotDavoliScore(DavoliScore = davoli$CharmOG_score,
                                  GenomicMeasure = davoli$chrArmpostWGDMean,
                                  x_lab = "Chromosome arm proportion loss post WGD",
                                  y_lab = "Chromosome arm oncogenes score",
                                  label.x.cor.pos = 0.7,
                                  chrNames = davoli$Arm,
                                  title = title
  )
  plist
}


#############################################################################################################

#Plot the main plot, with chromosomes arms before( haploid LOH ) and after WGD (neutral LOH)

load(file = "data/CNVinfo.rda")

load("data/manec.seg.sequenza.rda")

seg.mat.copy = seg.mat.copy.list$segments %>%
  filter(SampleID %in% 
           ( CNVinfo %>% filter( (Molecular_subtype %in% c("Nec", "Aca")) & GD.status == "GD") %>%
               pull(SampleID) %>%
               unique() )
  )

Manec.CINs = getEssential(seg.mat.copy)


plist = plotSegs(Manec.CINs$davoli)


pdf("figures/GD.Segs.Arms.Manec.pdf", width = 12, height = 8)

wrap_plots(plist, nrow = 2)

dev.off()


