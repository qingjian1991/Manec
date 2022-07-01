#see the segs status, LOH, 


###############################################################################################

# Figure 3C and Supplementary Figures:

# Make comparisons between MANEC-ACA, MANEC-NEC, STAD-CIN:GD and STAD-CIN:noGD. including:
#(1) proportional of copy neutral-LOH, haploid LOH
#(2) clonal and subclonal LOH
#(3) MATH scores (tumor heterogeneity scores)

###############################################################################################




library(ggpubr)
library(rstatix)


#' CNV_type
#' Test the CNV status according to the major and minor allele number.
#' @param nA: number of major allele copy number.
#' @param nB: number of minor allele copy number.

CNV_Type = function(nA, nB){
  
  TestCNV = function(nA, nB){
    
    if(is.na(nA) | is.na(nB)){
      return(NA)
    }
    
    if(nB == 0){
      if(nA == 0){
        "comp-Del"
      }else if(nA == 1){
        "Hap-LOH"
      }else{
        "neutral-LOH"
      }
    }else if(nA + nB > 2){
      "Amp"
    }else{
      "neutral"
    }
  }
  
  if( length(nA) != length(nB) ){
    stop( sprintf("#nA should be equal with #nB, #nA: %s, #nB: %s", length(nA), length(nB)) )
  }
  
  cnv = rep(NA, length(nA))
  
  for(i in 1:length(nA)){
    cnv[i] = TestCNV(nA[i], nB[i])
  }
  
  cnv
}


#' getCNV_summary
#' get summary of CNVs. Using seg.mat.copy.list$segments as the input
#' @param Cancer: names of the cancer
#' @param GD.status: dataframe of GD.status

getCNV_summary = function(segments, GD.status, Cancer = "Cancer" , chrSilent =  c("X","Y")){
  #load hg19 chromosomal information.
  
  #chr info for hg19
  chrSilent = c("X","Y")
  chrLens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
              158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
              114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
              63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  names(chrLens) = c(1:22, "X","Y")
  chrLens = chrLens[!names(chrLens) %in% chrSilent]
  
  if( ! "Molecular_subtype" %in% colnames(segments)){
    segments$Molecular_subtype = Cancer
  }
  
  segments = segments %>%
    mutate( CNV_Type =  CNV_Type(nA, nB)) %>%
    filter(
      !Chr %in% chrSilent,
      !is.na(CNV_Type)
    ) %>%
    mutate(
      width = End - Start
    ) %>%
    group_by(Molecular_subtype, SampleID, CNV_Type) %>%
    summarise(
      width = sum(width)
    ) %>%
    mutate(
      Cancer= Cancer
    )
  
  segments %>%
    group_by(Molecular_subtype, SampleID) %>%
    filter(CNV_Type != "neutral") %>%
    mutate( total_width = sum(width)  ) %>%
    mutate( prop_cnv = round(width/total_width, 4),
            prop_cnv_total = round(width/sum(chrLens), 4)
            ) %>%
    left_join(GD.status)
}


#####################################################################################
#load CNV Information.

CNVinfo = list()

#load MANEC
load("data/manec.seg.sequenza.rda")
load("data/mutdata_all.rda")


CNVinfo[[1]] = getCNV_summary(
              seg.mat.copy.list$segments %>% 
                 mutate(Molecular_subtype =  mapply(function(x) x[3], str_split(SampleID, "_")) ),
              GD.status = mutdata_all %>%
                select(patient, GD.status) %>%
                unique.data.frame() %>%
                setNames(c("SampleID","GD.status") ),
               Cancer = "MANEC"
               ) %>%
  filter(
    GD.status == "GD"
  )

#load STAD
load("data/STAD.ascat.rda")
#Get the GD status.
load("data/mutdata_stad.rda")

segmentsInfo = seg.mat.copy.list$segments

CNVinfo[[2]] = getCNV_summary(seg.mat.copy.list$segments,
                              GD.status = mutdata_stad %>%
                                select(patient, GD.status) %>%
                                unique.data.frame() %>%
                                setNames(c("SampleID","GD.status") ),
                              Cancer = "STAD"
               )
CNVinfo = purrr::reduce(CNVinfo, rbind)

CNVinfo = CNVinfo %>%
  filter(
    Molecular_subtype %in% c("Aca", "Nec", "CIN")
  ) %>%
  filter(!is.na(GD.status)) %>%
  mutate(CNV_Type = factor(CNV_Type, levels = c("neutral-LOH","Hap-LOH","comp-Del", "Amp") )) %>%
  mutate(subtype = str_c(Molecular_subtype, GD.status, sep = ":") ) %>%
  mutate(subtype = factor(subtype, levels = c("Aca:GD","Nec:GD","CIN:GD", "CIN:nGD") )  ) %>%
  filter(CNV_Type %in% c("neutral-LOH", "Hap-LOH") )

#non-params test: wilcox_test
wilcox.test.cnvinfo = CNVinfo  %>%
  group_by(CNV_Type) %>%
  wilcox_test(prop_cnv  ~ subtype) %>%
  add_xy_position(x = "subtype",  fun = "max" ) %>%
  #add_significance(p.col = "p") %>%
  filter(group2 == "CIN:nGD")

p1 = CNVinfo %>%
  ggplot(aes(x = subtype, y = prop_cnv) ) + theme_bw() +
  geom_boxplot(aes(fill = GD.status), size = 0.6 ) +
  facet_wrap( ~ CNV_Type) +
  labs(y = "Proportion", x = NULL) +
  scale_fill_manual( values = c("#5959AF","#878C8D") ) +
  theme(
    axis.text.x = element_text( size = 10, angle = 90, hjust = 1, vjust = 0.5 ),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 15),
    legend.position = "none"
  )+
  stat_pvalue_manual(wilcox.test.cnvinfo) 

p1




pdf("figures/GD.LOH_within_GD_boxplot-p.value.pdf", width = 5.47, height = 4 )
p1
dev.off()


##################################################################################

## see samples with and without WGD, their ITH.


head(CNVinfo)

load("data/heterogeneity.index.rda")

het.index


het = CNVinfo %>%
  mutate(Tumor_Sample_Barcode = ifelse( grepl("TCGA", SampleID), str_sub(SampleID, 1, 12), SampleID)
         ) %>%
  inner_join(
    het.index
  ) %>%
  ungroup() %>%
  select(Molecular_subtype, Tumor_Sample_Barcode, cancer, GD.status, MATH, shannon, simpson, Molecular.Subtype) %>%
  unique.data.frame()%>%
  filter( Molecular_subtype %in% c("Aca", "Nec", "CIN")) %>%
  mutate(subtype = str_c(Molecular_subtype, GD.status, sep = ":") ) %>%
  mutate(subtype = factor(subtype, levels = c("Aca:GD","Nec:GD","CIN:GD", "CIN:nGD") )  ) 


#non-params test: wilcox_test
wilcox.test.het = het  %>%
  wilcox_test(MATH  ~ subtype) %>%
  add_xy_position(x = "subtype",  fun = "max" ) %>%
  add_significance(p.col = "p") %>%
  filter(group2 == "CIN:nGD")

p.het = het %>%
  ggplot(aes(x = subtype, y = MATH ) ) + theme_classic() +
  geom_boxplot( size = 0.6, aes(fill = GD.status) ) +
  scale_fill_manual( values = c("#5959AF","#878C8D") ) +
  labs(x = NULL ) +
  theme(
    axis.text.x = element_text( size = 10, angle = 90, hjust = 1, vjust = 0.5 ),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 15),
    legend.position = "none"
  )+
  stat_pvalue_manual(wilcox.test.het, label = "p.signif") 



pdf("figures/GD.MATH_within_GD_boxplot-p.value.pdf", width = 4.4, height = 4.2)

p.het

dev.off()


############################################################################################

# see the proportion of clonal and subclonal LOH.

chrLens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
            158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
            114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
            63811651, 62435964, 46944323, 49691432, 154913754, 57772954)

#segStatus.rds, see SCNA_ITH.R
segStatus = readRDS("data/segStatus.rds")

clonalLOH = segStatus %>%
  ungroup() %>%
  filter(!seqnames %in% c("X","Y") ) %>%
  mutate(CNV_Type = CNV_Type(A, B)) %>%
  group_by(SampleID, CNV_Type, clonalStatus ) %>%
  summarise(
    width = sum(width)
  ) %>%
  mutate(prop_cnv_total = round(width/sum(chrLens), 4)) %>%
  mutate(Molecular_subtype =  mapply(function(x) x[3], str_split(SampleID, "_")))

clonalLOH = clonalLOH %>%
  filter(CNV_Type %in% c("neutral-LOH", "Hap-LOH") ) %>%
  filter(clonalStatus %in% c("Clonal","Subclonal") ) %>%
  filter(Molecular_subtype %in% c("Aca","Nec") ) %>%
  mutate(CNV_Type = factor(CNV_Type, levels = c("neutral-LOH","Hap-LOH"),
                           labels = c("Neutral-LOH","Hap-LOH")
                           ) ) 

#non-params test: wilcox_test
wilcox.test = clonalLOH  %>%
  group_by(CNV_Type, Molecular_subtype) %>%
  wilcox_test(prop_cnv_total ~ clonalStatus) %>%
  add_xy_position(x = "Molecular_subtype", fun = "max" ) %>%
  add_significance()

p2 = clonalLOH %>%
  ggboxplot(
    x = "Molecular_subtype", y = "prop_cnv_total", fill = "clonalStatus",
    facet.by = "CNV_Type", scales = "free",
    #add = "jitter", add.params = list(width = 0.1, size = 1, shape = 1) ,
    ggtheme = theme_bw() +
      theme(
        axis.text.x = element_text( size = 11, angle = 0 ),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 15)
      )
    
    ) +
  labs(x = NULL, y = "% genome LOH") +
  scale_fill_manual( values = c("#F4BD70","#878C8D") ) +
  stat_pvalue_manual(wilcox.test, label = "p.signif") +
  guides(fill = "none"  )

pdf("figures/GD.LOH_clonalStatus_boxplot.pdf", width = 5.47, height = 3.29 )
p2
dev.off()

#####################################################################################

save(CNVinfo, file = "data/CNVinfo.rda")


