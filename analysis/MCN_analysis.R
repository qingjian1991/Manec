#######################################################################################
#
# Fig 3B: The scatter plot of tumor ploidy and the fraction of autosomal tumor genomes with MCN >=2
#
#######################################################################################

#Purpose:
# According to Bielski et al.(2018), we could measure the fraction of patients with the Major Copy Number (MCN) in autosomal tumor genomes >=2

#References: Bielski, C. M., Zehir, A., Penson, A. V., Donoghue, M. T. A., Chatila, W., Armenia, J., . . . Taylor, B. S. (2018). Genome doubling shapes the evolution and prognosis of advanced cancers. Nat Genet, 50(8), 1189-1195. doi:10.1038/s41588-018-0165-1


library(ggsci)

#' getMCN
#' calculate the fraction of automsomal tumor genomes with an MCN of two or greater for all patients.
#' @param segments: segments file
#' @param GD.status: GD.status from method1
#' @param Cancer: the names of cancer subtypes.
#' @param chrSilent: Chromosomal to ignored.


getMCN = function(segments, GD.status,
                  Cancer = "Cancer", chrSilent =  c("X","Y")){
  #load hg19 chromosomal information.
  chrLens = data.frame(
    Chr = c(1:22, "X","Y"),
    ChrLen = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
               158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
               114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
               63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  )
  
  if( ! "Molecular_subtype" %in% colnames(segments)){
    segments$Molecular_subtype = Cancer
  }
  
  segments_chr = segments %>%
    filter(
      !Chr %in% chrSilent,
    ) %>%
    mutate(
      width = End - Start
    ) %>%
    mutate(
      MCN = ifelse(nA >=2, "Amp", "Neutral")
    ) %>%
    group_by(SampleID, Molecular_subtype, Chr, MCN, Ploidy) %>%
    summarise(
      width = sum(width)
    ) %>%
    left_join(
      chrLens
    ) 
  
  segments_genome = segments_chr %>%
    group_by(SampleID, Molecular_subtype, MCN, Ploidy) %>%
    summarise(
      width = sum(width)
    ) %>%
    pivot_wider(
      id_cols = c("SampleID","Molecular_subtype","Ploidy"),
      names_from = "MCN",
      values_from = "width",
      values_fill = 0
    ) %>%
    mutate(
      ratio = Amp/(Amp + Neutral)
    ) %>%
    full_join(
      GD.status
    )
  
  segments_genome %>%
    ggplot(aes(x= ratio, fill = Molecular_subtype) ) + theme_classic() +
    geom_histogram()
   
  return(
    segments_genome
  )
}


#######################################################################################

#load MANEC
load("data/manec.seg.sequenza.rda")
load("data/mutdata_all.rda")

#segments = seg.mat.copy.list$segments %>% 
#  mutate(Molecular_subtype =  mapply(function(x) x[3], str_split(SampleID, "_")) )

CNVinfo = list()

CNVinfo[[1]] = getMCN(
  seg.mat.copy.list$segments %>% 
    mutate(Molecular_subtype =  mapply(function(x) x[3], str_split(SampleID, "_")) ),
  GD.status = mutdata_all %>%
    select(patient, GD.status) %>%
    unique.data.frame() %>%
    setNames(c("SampleID","GD.status") )
) 

#load STAD
load("data/STAD.ascat.rda")
#Get the GD status.
load("data/mutdata_stad.rda")

CNVinfo[[2]] = getMCN(seg.mat.copy.list$segments,
                              GD.status = mutdata_stad %>%
                                select(patient, GD.status) %>%
                                unique.data.frame() %>%
                                setNames(c("SampleID","GD.status") )
)


CNVinfo = purrr::reduce(CNVinfo, rbind)

MCNInfo = CNVinfo %>%
  filter(
    Molecular_subtype %in% c("Aca", "Nec")
  )

p1 = MCNInfo %>%
  ggplot(aes(x= ratio, fill = Molecular_subtype) ) + theme_classic() +
  geom_histogram() +
  geom_vline(xintercept = 0.5, linetype = 2, size = 1, col = "red4") +
  geom_vline(xintercept = 0.4, linetype = 6, size = 0.6, col = "grey20") +
  scale_fill_d3() +
  labs(x = "% autosomal geneomes with MCN >=2", y = "Frequency") +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14)
  )


pdf("figures/GD.MCN.distribution.MANEC.pdf", width = 6, height = 4.2)

p1

dev.off()

MCNInfo1 = MCNInfo %>%
  filter(
    Molecular_subtype %in% c("Aca", "Nec")
  ) %>%
  mutate(
    GD2 = ifelse(ratio >= 0.5, "GD2", "nGD2")
  ) 

table(MCNInfo1$GD.status, MCNInfo1$GD2)
#     GD2  nGD2
#GD   91    4
#nGD   0    6


MCNInfo1 = MCNInfo %>%
  filter(
    Molecular_subtype %in% c("Aca", "Nec")
  ) %>%
  mutate(
    GD2 = ifelse(ratio >= 0.4, "GD3", "nGD3")
  ) 

table(MCNInfo1$GD.status, MCNInfo1$GD2)

#    GD3   nGD3
#GD   95    0
#nGD   1    5-



##############################################################################################

### Read the GD information directed from EstimatedClonality.

#purity info
purity_table = read.table("data/purityV1.csv", header = T, sep = ",", row.names = 1)
#rownames(purity_table) = purity_table$sample

purityInfo = purity_table %>%
  dplyr::select(sample, sequenza, ploidy.mean.cn, new_sample)  %>%
  dplyr::filter(!is.na(new_sample))


purityInfo # This is in GenomeDoubling.
purity_table


GD.infor = list()

for(i in purityInfo$new_sample){
  
  GD.infor[[i]] = read.table( sprintf("data/Manec/Manec/%s/%s.GD.tsv", i, i) , header = T, sep = "\t")

}

GD.infor = purrr::reduce(GD.infor, rbind)

GD.infor = GD.infor %>%
  mutate(Tumor = mapply(function(x) x[3], str_split(patient, "_") ) ) %>%
  filter(Tumor != "Squca") %>%
  left_join(purity_table %>% select(new_sample, ploidy.mean.cn, sample ), by = c("patient" = "new_sample"))


p1 = ggplot(GD.infor, aes(x = prop.major.even.obs, y = ploidy.mean.cn, col = GD.status, shape =  Tumor)) +
  geom_jitter(size =1.8) +
  theme_classic2()+
  geom_hline(yintercept = 2, linetype =2, size =1, col = "grey") +
  geom_vline(xintercept = 0.5, linetype =2, size =1, col = "red4") +
  labs(x = "% autosomal genomes with MCN >=2", y = "Tumor Ploidy") +
  theme(
    axis.title = element_text(size = 15)
  ) +
  scale_color_manual(values = c("#F8766D","#00BFC4") ) +
  scale_shape_manual(values= c(2,4) )

pdf("figures/GD.infor.pdf", width = 5.2, height = 3.6)
p1
dev.off()


### Read the GD timing detected from EstimatedClonality

GD.time = list()

for(i in purityInfo$new_sample){
  
  GD.time[[i]] = read.table( sprintf("data/Manec/Manec/%s/%s.GD.Timing.tsv", i, i) , header = T, sep = "\t")
  
}

GD.time = purrr::reduce(GD.time, rbind)

GD.time = GD.time %>%
  mutate(Tumor = mapply(function(x) x[3], str_split(sample, "_") ),
         simple = mapply(function(x) x[1], str_split(sample, "_"))
         ) %>%
  filter(Tumor != "Squca") %>%
  left_join(GD.infor %>% select(patient, GD.status, sample ), by = c("sample" = "patient")) %>%
  filter(GD.status != "nGD") %>%
  mutate(abs = before - after)

GD.time %>%
  pivot_wider(
    id_cols = c("sample", "Tumor"),
    names_from = source,
    values_from = abs
  ) %>%
  ggplot( aes(x = segments, y = arms, col = Tumor ) ) + theme_bw() +
  geom_abline(slope = 1, linetype =2, size = 1.2) +
  geom_point() 


GD.time = GD.time %>%
  filter(source == "segments") 

table(GD.time$abs >= 0 )
#FALSE  TRUE 
#61    34 

#Early Changes 61, Later Changes 34.


p1 = GD.time %>%
  mutate(sites = mapply(function(x) str_c(x[1], x[2], sep = "_"), str_split(sample, "_") )) %>%
  mutate(sites = ifelse(sites == "M34_1G","M34_1E", sites ) ) %>%
  pivot_wider(
    id_cols = sites,
    names_from = Tumor,
    values_from = abs
  ) %>%
  na.omit() %>%
  mutate(Type = ifelse(Aca <0 & Nec <0, "After", ifelse(Aca >0 & Nec >0, "Before", "Undetermined") ) ) %>%
  ggplot( aes(x = Aca, y = Nec, col =  Type) ) + theme_classic2() +
  geom_abline(slope = 1, linetype =2, size = 1.2) +
  geom_point(size = 2 ) +
  scale_color_manual(values = c("#E63915","#0083B6","#6E7B8B")) +
  labs(x = "Prop Loss before GD - Prop Loss after GD for ACA", y = "Prop Loss before GD - Prop Loss after GD for NEC" ) +
  geom_hline(yintercept = 0, linetype =2, size = 1.2, col = "grey") +
  geom_vline(xintercept = 0, linetype =2, size = 1.2, col = "grey") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13)
  )

pdf("figures/GD.timing.pdf", width = 5.2, height = 4.8)
p1
dev.off()
