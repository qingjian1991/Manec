############################################################################################

#######  Running TCGA data for tumor whole genome doubling (WGD) status .

#This code shows the analysis code of R package "EstimateClonality". 
#see code: github.com/qingjian1991/EstimateClonality

#This code has not been well checked. Please contact me if have any questions.


library(EstimateClonality)


#(1) for mutation data

load( sprintf("/data1/qingjian/Rproject/AminoAcids/ProcessMutations/MC3/TCGA-%s.MC3.rda", "STAD") )
#loading variable is maf

head(maf)

maf.select = maf %>%
  dplyr::select(c('Tumor_Sample_Barcode' ,'Chromosome' ,'Start_Position' ,'End_Position' ,'Reference_Allele' ,'Tumor_Seq_Allele2' ,'t_alt_count' ,'t_ref_count' ,'Hugo_Symbol' ,'Variant_Classification', "HGVSp_Short")) %>%
  mutate(
    Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, 1, 15)
  )

colnames(maf.select) = c('Patient' ,'Chr' ,'Start_position' ,'End_position' ,'Reference' ,'Alternate' ,'Variant_freq' ,'Ref_freq' ,'Hugo_Symbol' ,'Variant_Classification', 'Protein_Change')

write.table(maf.select, file = "datasets/STAD.maf.txt", sep = "\t", row.names = F, quote = F)


##### Read CNAs of COSMIC
# 
# cnas = read.csv("datasets/Cosmic_STAD_V92_37_CNA.csv", header = T)
# 
# #get purity information from cbioportal
# tcga_samples = read.delim("/data1/database/cbioportal/stad_tcga_pub/data_clinical_sample.txt", header = T,
#                           comment.char = "#", stringsAsFactors = F
#                           )
# 
# 
# cnas = inner_join(
#   cnas %>% select(SAMPLE_NAME, CHROMOSOME_G_START_G_STOP, TOTAL_CN, MINOR_ALLELE, MINOR_ALLELE) %>% unique.data.frame() , 
#   tcga_samples %>% select(SAMPLE_ID, PATIENT_ID, MOLECULAR_SUBTYPE , ABSOLUTE_EXTRACT_PURITY, ABSOLUTE_EXTRACT_PLOIDY),
#   by = c("SAMPLE_NAME" ="SAMPLE_ID")
# ) %>%
#   filter(!is.na(ABSOLUTE_EXTRACT_PURITY))
# 
# cnas = str_split(cnas$CHROMOSOME_G_START_G_STOP, "[:|..]", simplify  = TRUE)[,c(1,2,4)] %>%
#   as.data.frame() %>%
#   setNames(nm = c("Chr", "Start", "End")) %>%
#   cbind(cnas) %>%
#   setNames(c("Chr", "Start", "End", "SampleID", "CHROMOSOME_G_START_G_STOP", "cn", "nB", "PATIENT_ID", "MOLECULAR_SUBTYPE", "Aberrant Cell Fraction","Ploidy") ) %>%
#   mutate(nA = cn - nB,
#          nProbes = NA,
#          ) %>%
#   select(c("SampleID", "Chr", "Start", "End", "nProbes", "cn", "nA", "nB", "Ploidy", "Aberrant Cell Fraction", "PATIENT_ID","MOLECULAR_SUBTYPE")
# )
# 
#   
# seg.mat.copy.list = list()
# seg.mat.copy.list$segments = cnas
# 
# save(seg.mat.copy.list, file = "datasets/STAD.seg.rda")

#Read CNVs of TCGA.

#


#Download ASCAT from TCGA, hg38

#query
query <- TCGAbiolinks::GDCquery(project = "TCGA-STAD", 
                                data.category = "Copy Number Variation",
                                data.type = "Allele-specific Copy Number Segment")
#download
TCGAbiolinks::GDCdownload(query)

#get results
stad_ascat <- TCGAbiolinks::GDCprepare(query)

#hg38 to hg19
stad_ascat_hg19 = liftOverData(data = stad_ascat,
                               chainFiles = "hg38ToHg19.over.chain",
                               chrVersion = "hg19")

save(stad_ascat, file = "stad_ascat_hg38.rda")
#There a lot of small regions in the 
save(stad_ascat_hg19, file = "stad_ascat_hg19.rda")

stad_ascat_hg19_1 = stad_ascat_hg19 %>% filter(width >=1e5)

#grs =  stad_ascat_hg19$gr_muts_ch

#grs1 = grs[grs$Sample == "TCGA-BR-8296-10A-01D-2338-01", ]
#grs2 =reduce(grs1)

cnas = stad_ascat_hg19_1 %>%
  mutate(SAMPLE_ID = str_sub(Sample, 1, 15),
         Chromosome = str_remove(Chromosome, "chr")
         ) %>%
  inner_join(
    tcga_samples %>% select(SAMPLE_ID, PATIENT_ID, MOLECULAR_SUBTYPE , ABSOLUTE_EXTRACT_PURITY, ABSOLUTE_EXTRACT_PLOIDY) %>%
      filter( ! is.na(ABSOLUTE_EXTRACT_PURITY )),
    by = c("SAMPLE_ID" ="SAMPLE_ID")
  ) %>%
  mutate(nProbes = NA) %>%
  select(
    SAMPLE_ID, Chromosome, Start, End, nProbes, Copy_Number, Major_Copy_Number, Minor_Copy_Number, ABSOLUTE_EXTRACT_PLOIDY, ABSOLUTE_EXTRACT_PURITY, MOLECULAR_SUBTYPE,
  ) %>%
  setNames(nm = c("SampleID", "Chr", "Start", "End", "nProbes", "cn", "nA", "nB", "Ploidy", "Aberrant Cell Fraction", "Molecular_subtype"))

seg.mat.copy.list = list()
seg.mat.copy.list$segments = cnas

save(seg.mat.copy.list, file = "datasets/STAD.ascat.rda")

#see overlap between mutations and cnvs

patients = intersect(cnas$SampleID, maf.select$Patient) 

save(patients, file = "datasets/stad.patients.rda")

#Running 
sc.runs = easypar::run(
  FUN = function(i){
    manec.mut = "datasets/STAD.maf.txt"
    manec.seg = "datasets/STAD.ascat.rda"
    message( patients[i] )
    EstimateClonality::clonality.estimation(mutation.table.loc= manec.mut,
                         seg.mat.loc= manec.seg,
                         data.type='STAD',
                         TCGA.barcode= patients[i],
                         ANALYSIS.DIR= "STAD/",
                         min.alt.reads = 4,	
                         min.depth = 20	
    )
  },
  PARAMS = lapply( 1:length(patients) , list),
  parallel = TRUE,
  progress_bar = TRUE,
  packages = "EstimateClonality",
  export = c("patients")
)

#############################################################################################

#Read datas and see the results.

mutdata_stad = list()

for(i in patients){
  mutdata_stad[[i]] = read.table( sprintf("STAD/STAD/%s/%s.earlylate.tsv", i, i), header = T, sep = "\t" )
}

mutdata_stad = purrr::reduce(mutdata_stad, rbind)

mutdata_stad = mutdata_stad %>%
  #join molecular subtypes.
  left_join(
    tcga_samples %>% select(SAMPLE_ID, MOLECULAR_SUBTYPE , ABSOLUTE_EXTRACT_PURITY, ABSOLUTE_EXTRACT_PLOIDY),
    by = c("patient" = "SAMPLE_ID")
  )%>%
  #see clonal status from single samples.
  mutate(clonal.status =  ifelse(prob.clonal > prob.subclonal, "Clonal", "Subclonal")  ) %>%
  #join mutations information
  left_join(
    maf %>%
      mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, 1, 15)) %>%
      mutate(mutation_id = str_c( Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, sep = ":" )) %>%
      select(mutation_id, Hugo_Symbol, Variant_Classification, Variant_Type, HGVSp_Short, SIFT, PolyPhen, IMPACT, VARIANT_CLASS) ,
    by = c("mutation_id" = "mutation_id")
    
  )

save(mutdata_stad, file = "datasets/mutdata_stad.rda")

#see some summary.

table( mutdata_stad$clonal.status)
#Clonal Subclonal 
#17294     48738

stad.GD = mutdata_stad %>% 
  select(patient, GD.status, MOLECULAR_SUBTYPE) %>%
  unique.data.frame() %>%
  group_by(GD.status, MOLECULAR_SUBTYPE) %>%
  summarise( num = n())

stad.GD
# A tibble: 8 x 3
# Groups:   GD.status [2]
#GD.status MOLECULAR_SUBTYPE   num
#<chr>     <chr>             <int>
#1 GD        CIN                  47
#2 GD        EBV                   5
#3 GD        GS                    6
#4 GD        MSI                   9
#5 nGD       CIN                  15
#6 nGD       EBV                   8
#7 nGD       GS                   23
#8 nGD       MSI                  20

#Summary Number of GD status in the STAD.

#Percentages of GD patients in the CIN sub groups.
#47/(15+47)
#0.7580645
GD.status.label = stad.GD %>%
  group_by(MOLECULAR_SUBTYPE) %>%
  mutate(total_num = sum(num),
         precentage = round(num/total_num, 2)
         ) 

pdf("figures/STAD.GD.barplot.pdf", width = 4, height = 3)

GD.status.label %>%
  ggplot(aes(x = MOLECULAR_SUBTYPE, y = num, fill = GD.status ) ) + theme_pubr() +
  geom_bar(stat = "identity", position = position_dodge2(), 
           width = 0.8) +
  labs(x = "Molecular Subtypes from STAD", y = "# of cases") +
  geom_text( aes(
    label = precentage),
    position = position_dodge2(0.8),
    vjust = 1.5,
    size = 3
    )

dev.off()

###############################################################################################

#focus on the CIN and GD tumors.

## see TP53 and RB1

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


driverStatus = mutdata_stad %>%
  #select STAD CIN
  filter(MOLECULAR_SUBTYPE == "CIN") %>%
  filter(Hugo_Symbol %in% drivers) %>%
  #removed non-GD 
  filter( GD.status != "nGD") %>%
  filter( !(comb.timing == "Early" & clonal.status == "Subclonal") ) %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site") ) %>%
  group_by(Hugo_Symbol, comb.timing, clonal.status) %>%
  summarise(
    num = n()
  ) %>%
  mutate(
    class = str_c(comb.timing, clonal.status, sep = "_")
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
  pull(Hugo_Symbol) %>%
  unique()

pdf("figures/Drivers_Clonal_Status.STAD_CIN.pdf", width = 9, height = 5.5)

driverStatus %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genesOrder)) %>%
  ggplot(aes(x = Hugo_Symbol, y = num)) + theme_classic() +
  geom_bar(aes(fill = class), stat = "identity", position =   position_fill(), width = 0.8 ) +
  ggpubr::theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14, hjust = 1, angle = 90, vjust = 0.5)
  )

dev.off()





  