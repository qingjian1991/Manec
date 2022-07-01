############################################################################################
## Estimate tumor whole genome doubling (WGD) status.
############################################################################################


#This code shows the analysis code of R package "EstimateClonality". 
#see code: github.com/qingjian1991/EstimateClonality

#This code has not been well checked. Please contact me if have any questions.

library(EstimateClonality)
library(tidyverse)


blac.mut = system.file("extdata/Input/BLCA.mutation.table.txt", package = "EstimateClonality")
blac.seg = system.file("extdata/Input/tcga.blca.seg.hg19.rdata", package = "EstimateClonality")


clonality.estimation(mutation.table.loc= blac.mut,
                     seg.mat.loc= blac.seg,
                     data.type='TCGA_BLCA',
                     TCGA.barcode="TCGA-BT-A42C",
                     ANALYSIS.DIR="example1/")

mutdata = read.table(blac.mut, header = T, sep = "\t")

load(blac.seg)

blac.seg = get(seg.mat.copy.list)

#seg.mat.copy <- load(seg.mat.loc)
#seg.mat.copy.list <- get(seg.mat.copy)
#seg.mat.copy <- seg.mat.copy.list$segments
#barcodes <- intersect(unique(seg.mat.copy$SampleID), unique(mutation.table[, 1]))

####################################################################################################################

#Prepare input.

#for mutations.

#all mutations, after filtering. 
#maf.all is the mutations after filering, avoiding false positive.
#maf.all includes all mutations, including Nec, Aca, Squaca, Lym and .
maf.all =  readRDS("../isma/maftools.data/manec.maf.rds")

maf.all = maf.all %>%
  dplyr::select(c('Tumor_Sample_Barcode' ,'Chromosome' ,'Start_Position' ,'End_Position' ,'Reference_Allele' ,'Tumor_Seq_Allele2' ,'Alt_allele_depth' ,'Ref_allele_depth' ,'Hugo_Symbol' ,'Variant_Classification'))


colnames(maf.all) = c('Patient' ,'Chr' ,'Start_position' ,'End_position' ,'Reference' ,'Alternate' ,'Variant_freq' ,'Ref_freq' ,'Hugo_Symbol' ,'Variant_Classification')

write.table(maf.all, file = "datasets/manec.maf.all.txt", sep = "\t", row.names = F, quote = F)


#for cnvs.

#format is sequenza

#purity info
purity_table = read.table("../ITH/purityV1.csv", header = T, sep = ",", row.names = 1)
#rownames(purity_table) = purity_table$sample

purityInfo = purity_table %>%
  dplyr::select(sample, sequenza, ploidy.mean.cn, new_sample)  %>%
  dplyr::filter(!is.na(new_sample))
 
clinical = read.table("../isma/clinical.20210508.txt", header = T, sep = "\t")

if(formate == "sequenza"){
  #for CMVs Path
  cnvPath="/data1/qingjian/Project/QiuMZ/sequenza/Seg" 
  cnvfiles = list.files(path = cnvPath)
  
  cnvlist = list()
  
  for(i in 1:nrow(purityInfo)){
    sampleNames = purityInfo$new_sample[i]
    cnv = read.table( sprintf("%s/%s", cnvPath, str_c(purityInfo$sample[i], "_segments.txt") ), header = T, sep = "\t" ) %>%
      cbind( purityInfo %>% filter(new_sample == sampleNames) ) %>%
      dplyr::select(new_sample, chromosome, start.pos,  end.pos, N.BAF, CNt, A, B, ploidy.mean.cn, sequenza) %>%
      setNames( nm = c("SampleID", "Chr", "Start", "End", "nProbes", "cn", "nA", "nB", "Ploidy", "Aberrant Cell Fraction") )
    
    cnvlist[[sampleNames]] = cnv
  }
  
  cnvlist = purrr::reduce(cnvlist, rbind) %>%
    mutate(Chr = str_remove(Chr, "chr"))
  
  #store in seg.mat.copy.list
  
  seg.mat.copy.list = list()
  seg.mat.copy.list$segments = cnvlist
  seg.mat.copy.list$sampleInfo = purityInfo
  
  save(seg.mat.copy.list, file = "datasets/manec.seg.sequenza.rda")
}

################################################################################

#run as examples.

#Running 
sc.runs = easypar::run(
  FUN = function(i){
    manec.mut = "datasets/manec.maf.all.txt"
    manec.seg = "datasets/manec.seg.sequenza.rda"
    message( purityInfo$new_sample[i] )
    clonality.estimation(mutation.table.loc= manec.mut,
                         seg.mat.loc= manec.seg,
                         data.type='Manec',
                         TCGA.barcode= purityInfo$new_sample[i],
                         ANALYSIS.DIR= "Manec/",
                         min.alt.reads = 4,	
                         min.depth = 20	
    )
  },
  PARAMS = lapply( 1:105 , list),
  parallel = TRUE,
  progress_bar = TRUE,
  packages = "EstimateClonality",
  export = c("purityInfo")
)

###################################################################################

#combined all samples.


mutdata = list()

for(i in purityInfo$new_sample){
  mutdata[[i]] = read.table( sprintf("Manec/Manec/%s/%s.earlylate.tsv", i, i), header = T, sep = "\t" )
}

mutdata = purrr::reduce(mutdata, rbind)

## get the GD status

mutdata_combined = maf.all %>%
  mutate(mutation_id = str_c(Patient, Chr, Start_position, Reference, sep = ":")) %>%
  right_join(
    mutdata
  ) %>%
  left_join(
    clinical %>% dplyr::select(Tumor_Sample_Simple, Simple, Tumor_ID), by = c("Patient" = "Tumor_Sample_Simple")
  )%>%
  #see clonal status from single samples.
  mutate(clonal.status =  ifelse(prob.clonal > prob.subclonal, "Clonal", "Subclonal")  )


mutdata_combined %>%
  select(Patient, GD.status, Simple, Tumor_ID) %>%
  unique.data.frame() %>%
  filter(
    GD.status != "GD"
  )


#   Patient GD.status Simple Tumor_ID
#1  M14_1D_Aca       nGD    M14      Aca
#2  M16_13H_Aca       nGD    M16      Aca
#3  M22_1F_Aca       nGD    M22      Aca
#4  M23_1B_Aca       nGD    M23      Aca
#5  M29_1P_Aca       nGD    M29      Aca
#6  M29_1P_Nec       nGD    M29      Nec



####################################################################################

#analysis.

#Plot Drivers according to the timing of drivers.

#see only driver genes.
drivers = c("TP53",
            "FAT4","FAT3","FAT2","FAT1", 
            "APC", "APC2", 
            "LRP1B","LRP1","LRP2",
            "CTNNB1",
            "KMT2D","KMT2A","KMT2B","KMT2C",
            "CDH11",
            "REBB2","ERBB3","ERBB4",
            "FBXW7","GNAS", "RB1", "AMER1",
            "PTPRT","SYNE1","BIRC6","TCF4", "RNF213","ZNF521","ABL1"
)

#combined mutations of all regions to call clonal and sub-clonal.

#get Number of Samples for each tumor.
nSamples = clinical %>%
  filter(Tumor_ID %in% c("Aca","Nec")) %>%
  group_by(Simple) %>%
  summarise(nSamples = n())

#get the clonal status by exam the overlapping mutations between regions.
mutClonalStatus =  mutdata_combined %>%
  filter(Tumor_ID != "Squca") %>%
  mutate(mutation_id = str_c( Simple, Chr, Start_position, Reference, sep = ":" ) ) %>%
  group_by(
    Simple, mutation_id
  ) %>%
  summarise(
    nMut = n()
  ) %>%
  left_join(nSamples) %>%
  mutate(clonal.status1 = ifelse(nMut >=7/8*nSamples, "Clonal","Subclonal"  ))
  
#combined data. 
mutdata_all = mutdata_combined %>%
  filter(Tumor_ID != "Squca") %>%
  mutate(mutation_id = str_c( Simple, Chr, Start_position, Reference, sep = ":" ) ) %>%
  left_join(mutClonalStatus %>% select(mutation_id, clonal.status1), 
            by = c("mutation_id" = "mutation_id")) %>%
  mutate(clonal.status2 = ifelse( (clonal.status == "Clonal" |  clonal.status1 == "Clonal"  ), "Clonal","Subclonal"  )  )
  
#There are more clonal mutations when by MRS.

#clonal status by single samples.
#table(mutdata_all$clonal.status)
#Clonal Subclonal 
#7602     25011 

#clonal status by exam overlapping mutations between regions.
# table(mutdata_all$clonal.status1)
#Clonal Subclonal 
#12677  21357

#table(mutdata_all$clonal.status2)
#Clonal Subclonal 
#12677     19936

#see status in drivers.


driverStatus = mutdata_all %>%
  #filter(Tumor_ID == "Nec") %>%
  filter(Hugo_Symbol %in% drivers) %>%
  #removed non-GD 
  filter( GD.status != "nGD") %>%
  filter( !(comb.timing == "Early" & clonal.status2 == "Subclonal") ) %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site") ) %>%
  group_by(Hugo_Symbol, comb.timing, clonal.status2) %>%
  summarise(
    num = n()
  ) %>%
  mutate(
    class = str_c(comb.timing,clonal.status2, sep = "_")
  )

genesOrder = driverStatus %>%
  reshape2::dcast(
    Hugo_Symbol  ~ class,
    value.var = "num",
    fill = 0
  ) %>%
  mutate(
    p_early_clonal = Early_Clonal/(Early_Clonal + Late_Clonal + Late_Subclonal),
    p_late_clonal = Late_Clonal/(Early_Clonal + Late_Clonal + Late_Subclonal),
    p_late_subclonal = Late_Subclonal/(Early_Clonal + Late_Clonal + Late_Subclonal),
  ) %>%
  arrange(desc(p_early_clonal), desc(p_late_clonal), desc(p_late_subclonal)) %>%
  #group drivers
  mutate(
    Genes_Group = ifelse(
      p_early_clonal > 0.4, "Pre-WGD",
      ifelse(p_early_clonal >0, "Post-WGD", "Later")
    )
    
  )


pdf("figures/Drivers_Clonal_Status.pdf", width = 9.3, height = 3.5)

driverStatus %>%
  left_join( genesOrder %>% select(Hugo_Symbol, Genes_Group)  ) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genesOrder$Hugo_Symbol)) %>%
  mutate(Genes_Group = factor(Genes_Group, levels = c("Pre-WGD","Post-WGD","Later") ) ) %>%
  ggplot(aes(x = Hugo_Symbol, y = num)) + theme_classic() +
  geom_col(aes(fill = class), position = position_fill(reverse = TRUE), width = 0.8 ) +
  geom_text(aes(label = num), position = position_fill(vjust = 1 ), size = 3  ) +
  facet_grid( ~ Genes_Group, scales = "free_x", space = "free" ) +
  ggpubr::theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14, hjust = 1, angle = 90, vjust = 0.5)
  ) +
  scale_fill_manual(values = c("#C7543C","#ACC9D8","#5B88C3") )

dev.off()


#Plot as whole and facet them later.

driverStatus_split = mutdata_all %>%
  #filter(Tumor_ID == "Nec") %>%
  filter(Hugo_Symbol %in% drivers) %>%
  #removed non-GD 
  filter( GD.status != "nGD") %>%
  filter( !(comb.timing == "Early" & clonal.status2 == "Subclonal") ) %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site") ) %>%
  group_by(Hugo_Symbol, comb.timing, Tumor_ID, clonal.status2) %>%
  summarise(
    num = n()
  ) %>%
  mutate(
    class = str_c(comb.timing,clonal.status2, sep = "_")
  )


pdf("figures/Drivers_Clonal_Status.Facte.pdf", width = 9.4, height = 4.7)


driverStatus_split %>%
  left_join( genesOrder %>% select(Hugo_Symbol, Genes_Group)  ) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genesOrder$Hugo_Symbol)) %>%
  mutate(Genes_Group = factor(Genes_Group, levels = c("Pre-WGD","Post-WGD","Later") ) ) %>%
  ggplot(aes(x = Hugo_Symbol, y = num)) + theme_classic() +
  geom_col(aes(fill = class), position = position_fill(reverse = TRUE), width = 0.8 ) +
  geom_text(aes(label = num), position = position_fill(vjust = 1 ), size = 3  ) +
  facet_grid( Tumor_ID ~ Genes_Group, scales = "free_x", space = "free" ) +
  ggpubr::theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14, hjust = 1, angle = 90, vjust = 0.5)
  ) +
  scale_fill_manual(values = c("#C7543C","#ACC9D8","#5B88C3") )


dev.off()




#mutdata_all: GD status of all mutations.
mutdata_all = mutdata_all %>%
  mutate(Simple.x = NULL) %>%
  dplyr::rename(Simple = Simple.y)

save(mutdata_all, file = "datasets/mutdata_all.rda")








