#########################################################################


load("datasets/Proportional.rda")

#See the ITH of SNVs or CNVs between Nec and Aca

#snvStatusReduce is in SNV_ITH.R

#Reduce MRS data into one.
snvStatusReduce  = snvITH %>%
  filter(Tumor_ID != "Squca") %>%
  group_by(Tumor_ID, Simple) %>%
  summarise(Clonal = mean(Clonal),
            Subclonal = mean(Subclonal),
            Clonal1 = mean(Clonal1),
            Subclonal1 = mean(Subclonal1),
            Clonal2 = mean(Clonal2),
            Subclonal2 = mean(Subclonal2)
  )%>%
  mutate(
    ITH = Subclonal/(Subclonal + Clonal),
    ITH1 = Subclonal1/(Subclonal1 + Clonal1),
    ITH2 = Subclonal2/(Subclonal2 + Clonal2),
    #TMB1 = Clonal1 + Subclonal1,
    TMB = Clonal + Subclonal
  )


summary(snvStatusReduce$Clonal)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.0    47.0    92.5   110.7   141.8   615.0




head(snvStatusReduce)

snvStatusLong = snvStatusReduce %>%
  mutate(ITH = 1 - ITH) %>%
  reshape2::dcast(
    Simple ~ Tumor_ID,
    value.var = "ITH"
  )

p2 = Jerry::ggpoints(snvStatusLong, 
                     x = "Aca", y = "Nec", 
                     annotations = c("cor"),
                     linetypes = "xy.line"
                     ) +
  labs(x = "T index of ACA", y = "T index of NEC")


#Summary of ITHs

snvStatusReduce %>%
  filter(Tumor_ID == "Aca") %>%
  pull(ITH1) %>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.04609 0.20889 0.33245 0.36784 0.53435 0.76389 

snvStatusReduce %>%
  filter(Tumor_ID == "Nec") %>%
  pull(ITH1) %>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03427 0.20394 0.33704 0.36886 0.49485 0.82222 


##########################################################################################
library(ggpmisc)

estimate_type = "nLOH"

#estimate_type = "total2"

segStatusReduce  = segStatus_wider %>%
  filter(Tumor_ID != "Squca") %>%
  group_by(Tumor_ID, Simple) %>%
  summarise(Clonal = mean( get( paste("Clonal", estimate_type, sep = "_") ) ),
            Subclonal = mean( get(  paste("Subclonal", estimate_type, sep = "_") ) ),
            ITH = mean( get(  paste("ITH", estimate_type, sep = "_") ) ),
            GII = mean( get(  paste("GII", estimate_type, sep = "_") ) )
  )

segStatusLong = segStatusReduce %>%
  filter(!Simple %in% c("M16","M23","M29", "M22") ) %>%
  reshape2::dcast(
    Simple ~ Tumor_ID,
    value.var = "ITH"
  )

p1 = Jerry::ggpoints(segStatusLong, x = "Aca", y = "Nec", annotations = c("xy.line","cor")) +
  labs(x = "Prop of subclonal neu-LOH of Aca", y = "Prop of subclonal neu-LOH of Nec")


pdf( sprintf("plot/ITHs_neutal-LOH_CNVs_SNVs.pdf", formate) , width = 6.8, height = 3.9)
p2+ p1
dev.off()


segStatusReduce %>% filter(!Simple %in% c("M23", "M16", "M29")) %>%
  pull(Clonal ) %>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1319  0.2519  0.3575  0.3795  0.5052  0.8170 

#Summary of ITHs

segStatusReduce %>%
  filter(Tumor_ID == "Aca") %>%
  pull(ITH) %>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.00000 0.01033 0.01920 0.06752 0.06691 0.38303       1  

segStatusReduce %>%
  filter(Tumor_ID == "Nec") %>%
  pull(ITH) %>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.00000 0.00595 0.07576 0.15478 0.16681 0.99540       1 


########################################################################################

#Another way to show the data.

segStatusReduce %>%
  select(Tumor_ID, Simple, Subclonal) %>%
  melt(
    id.vars = c("Tumor_ID","Simple")
  ) %>%
  ggplot(aes(x = Tumor_ID, y = value, col = Tumor_ID) ) + theme_pubr() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(position = position_jitterdodge() ) +
  ylim(0, 0.10) +
  labs(x = NULL, y = "")


  #
  #geom_line(aes(group = Simple), alpha = 0.2, linetype = 2, size = 0.5 )


snvStatusReduce %>%
  select(Tumor_ID, Simple, Subclonal) %>%
  melt(
    id.vars = c("Tumor_ID","Simple")
  ) %>%
  ggplot(aes(x = Tumor_ID, y = value, col = variable) ) + theme_pubr() +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(position = position_jitterdodge() ) +
  ylim(0, 300) +
  labs(x = NULL, y = "")





####################################################################################################

#See the correlation between neutral CNV-ITH and SSNV-ITH.

library(ggpubr)
library(rstatix)

#snv ITH information.
head(snvStatusLong)

ITH.data = inner_join( p.surv.Aca.nLOH$survData %>%
                         select(Patient.ID, ITH) %>%
                         setNames(nm = c("Patient.ID","ITH_Aca") ),
                       
                       p.surv.Nec.nLOH$survData %>%
                         select(Patient.ID, ITH) %>%
                         setNames(nm = c("Patient.ID","ITH_Nec") )
) %>%
  inner_join(snvStatusLong,
             by = c("Patient.ID" = "Simple")
  )


ITH.data = data.frame(
  ITH_nLOH = c(ITH.data$ITH_Aca, ITH.data$ITH_Nec),
  ITH_SNV = c(ITH.data$Aca, ITH.data$Nec),
  Cancer = rep( c("Aca","Nec"), each = length(ITH.data$Aca)  )
) %>%
  mutate(
    ITH_nLOH = ifelse(ITH_nLOH == "L", "Low", "High")
  )


wilcox.test = ITH.data %>%
  group_by(Cancer) %>%
  wilcox_test(ITH_SNV  ~ ITH_nLOH ) %>%
  add_xy_position(x = "ITH_nLOH", fun = "mean_sd")

pdf("plot/ITH.SNVs.CNVs.boxplot.pdf", width = 5.69, height = 3.38)

ITH.data %>%
  ggboxplot(
    x = "ITH_nLOH", y = "ITH_SNV", col = "ITH_nLOH",
    xlab = "ITH of Neutral-LOH", ylab = "ITH of SSNVs",
    facet.by = "Cancer", scales = "free",
    add = "jitter", add.params = list(width = 0.1, size = 1, shape = 1) ,
    ggtheme = theme_classic()) +
  stat_pvalue_manual(wilcox.test, label = "p") +
  guides(color = "none"  )

dev.off()

############################################################################################
#Define the ITH using different ways.

#ITH of nLOH.
ITH.data = inner_join( p.surv.Aca.nLOH$survData %>%
                         select(Patient.ID, ITH) %>%
                         setNames(nm = c("Patient.ID","ITH_Aca") ),
                       
                       p.surv.Nec.nLOH$survData %>%
                         select(Patient.ID, ITH) %>%
                         setNames(nm = c("Patient.ID","ITH_Nec") )
) %>%
  inner_join(snvStatusLong,
             by = c("Patient.ID" = "Simple")
  )

#load mutation ITH. 

load("../isma/maftools.data/heterogeneity.index.rda")


mut.het.index = het.index %>% 
  filter(cancer == "MANEC") %>%
  mutate(Patient.ID = mapply(function(x) x[1], str_split(Tumor_Sample_Barcode, "_") ) ) %>%
  left_join(ITH.data) %>%
  filter( !is.na(ITH_Aca) )


mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = MATH, fill = ITH_Aca )) + theme_classic() +
  geom_boxplot()

mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = MATH, fill = ITH_Nec)) + theme_classic() +
  geom_boxplot()


mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = shannon  , fill = ITH_Aca )) + theme_classic() +
  geom_boxplot()

mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = shannon  , fill = ITH_Nec)) + theme_classic() +
  geom_boxplot()



mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = simpson  , fill = ITH_Aca )) + theme_classic() +
  geom_boxplot()

mut.het.index %>%
  ggplot(aes(x = Molecular.Subtype, y = simpson  , fill = ITH_Nec)) + theme_classic() +
  geom_boxplot()


###############################################################################################

neutral.test = readRDS("datasets/neutral.test.rds")

head(ITH.data)

neutral.test.data = neutral.test %>%
  mutate(Patient.ID = mapply(function(x) x[1], str_split(sample, "_") ) ) %>%
  left_join(ITH.data) %>%
  filter( !is.na(ITH_Aca) )

neutral.test.data %>%
  ggplot(aes(x = Tumor_ID , y = Rsq , fill = ITH_Aca )) + theme_classic() +
  geom_boxplot()

neutral.test.data %>%
  ggplot(aes(x = Tumor_ID , y = Rsq , fill = ITH_Nec)) + theme_classic() +
  geom_boxplot()

neutral.test.data %>%
  ggplot(aes(x = Tumor_ID , y = rAUC , fill = ITH_Aca )) + theme_classic() +
  geom_boxplot()

neutral.test.data %>%
  ggplot(aes(x = Tumor_ID , y = rAUC , fill = ITH_Nec)) + theme_classic() +
  geom_boxplot()



neutral.test.data %>%
  filter(Tumor_ID == "Aca") %>%
  ggplot(aes(x = ITH_Aca , fill =  Rsq.test, group =  Rsq.test)) + theme_classic() +
  geom_bar( stat = "count" , position = position_fill(), width = 0.6)


neutral.test.data %>%
  filter(Tumor_ID == "Nec") %>%
  ggplot(aes(x = ITH_Nec , fill =  Rsq.test )) + theme_classic() +
  geom_bar( stat = "count",position = position_fill(), width = 0.6)


neutral.test.data1 = rbind(
  neutral.test.data %>%
    filter(Tumor_ID == "Aca") %>%
    mutate(ITH_group = ITH_Aca) %>%
    select(Rsq, rAUC, Simple, Tumor_ID, Rsq.test, ITH_group),
  
  neutral.test.data %>%
    filter(Tumor_ID == "Nec") %>%
    mutate(ITH_group = ITH_Nec) %>%
    select(Rsq, rAUC, Simple, Tumor_ID, Rsq.test, ITH_group)
)

neutral.test.data1 %>%
  mutate(ITH_group = factor(ITH_group, levels = c("L","H"), labels = c("Low","High") )) %>%
  ggplot(aes(x = ITH_group , fill =  Rsq.test )) + theme_classic() +
  geom_bar( stat = "count",position = position_fill(), width = 0.6) +
  facet_wrap( ~ Tumor_ID )



pdf("plot/Rsp_nLOH-groups.pdf", width = 4.8, height = 3.2)

neutral.test.data1 %>%
  group_by(Simple, Tumor_ID, ITH_group) %>%
  summarise(Rsq = mean(Rsq)) %>%
  mutate(ITH_group = factor(ITH_group, levels = c("L","H"), labels = c("Low","High") )) %>%
  ggplot(aes(x = ITH_group , y = Rsq , col =  ITH_group)) + theme_classic() +
  geom_boxplot() +
  geom_jitter(width = 0.1, position = ) +
  facet_wrap( ~ Tumor_ID ) +
  ylim(0.9, 1)

dev.off()

save.image("Prop.work.image.rda")


##########################################################################################################
#see the stage in the tumors.

stage = inner_join( p.surv.Aca.nLOH$survData %>%
              dplyr::select(Patient.ID, ITH) %>%
              setNames(nm = c("Patient.ID","ITH_Aca") ),
            
            p.surv.Nec.nLOH$survData %>%
              dplyr::select(Patient.ID, ITH) %>%
              setNames(nm = c("Patient.ID","ITH_Nec") )
) %>%
  left_join(
    clinical %>% dplyr::select(Simple , TNM_levels) %>% unique.data.frame(),
    by = c("Patient.ID" = "Simple")
  ) %>%
  mutate(TNM = str_remove(TNM_levels, pattern = "[A|B|C]"))


table(stage$TNM)

stage %>%
  filter( !is.na(ITH_Aca)) %>%
  melt(id.vars = c("Patient.ID","TNM_levels","TNM")) %>%
  mutate(value = factor(value, levels = c("L","H"), labels = c("Low","High") ) ) %>%
  group_by(
    variable, value, TNM 
  ) %>%
  summarise(num = n())

# 1 ITH_Aca  Low   I         3
# 2 ITH_Aca  Low   II        1
# 3 ITH_Aca  Low   III       9
# 4 ITH_Aca  High  II        5
# 5 ITH_Aca  High  III      11
# 6 ITH_Aca  High  IV        1
# 7 ITH_Nec  Low   I         3
# 8 ITH_Nec  Low   II        3
# 9 ITH_Nec  Low   III      11
# 10 ITH_Nec  Low   IV        1
# 11 ITH_Nec  High  II        3
# 12 ITH_Nec  High  III       9

pdf("plot/ITHs-nLOH_Stages-bestcutff-samecutoff.pdf", width = 5, height = 3.5 )

stage %>%
  filter( !is.na(ITH_Aca)) %>%
  melt(id.vars = c("Patient.ID","TNM_levels","TNM")) %>%
  mutate(value = factor(value, levels = c("L","H"), labels = c("Low","High") ) )  %>%
  ggplot(aes(x = value, fill = TNM) ) + theme_classic2() +
  labs(x = "ITH of Neutral-LOH", y = "Percentage of cases") +
  geom_bar(position = position_fill(), width = 0.5 ) +
  facet_wrap(~ variable) +
  scale_fill_brewer() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

dev.off()



p1 = stage %>%
  filter( !is.na(ITH_Aca)) %>%
  mutate(ITH_Aca = factor(ITH_Aca, levels = c("L","H"), labels = c("Low","High") ) ) %>%
  ggplot(aes(x = ITH_Aca, fill = TNM) ) + theme_pubr() +
  labs(x = "ITH of Neutral-LOH", y = "Percentage of cases", title = "Aca") +
  geom_bar(position = position_fill(), width = 0.5 ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p2 = stage %>%
  filter( !is.na(ITH_Nec)) %>%
  mutate(ITH_Nec = factor(ITH_Nec, levels = c("L","H"), labels = c("Low","High") ) ) %>%
  ggplot(aes(x = ITH_Nec, fill = TNM) ) + theme_pubr() +
  labs(x = "ITH of Neutral-LOH", y = "Percentage of cases", title = "Nec") +
  geom_bar(position = position_fill(), width = 0.5 ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p1 + p2






#########################################################################################################

#Plot Heamaps of ITH of SSNVs and CNVs.

head(snvStatusReduce)

#Tumor_ID Simple Clonal Subclonal Clonal1 Subclonal1 Clonal2 Subclonal2   ITH   ITH1  ITH2   TMB
#<chr>    <chr>   <dbl>     <dbl>   <dbl>      <dbl>   <dbl>      <dbl> <dbl>  <dbl> <dbl> <dbl>
#1 Aca      M01      21        123       34       110     21         110  0.854 0.764  0.840   144
#2 Aca      M02     615       3129     3473       258    615         258  0.836 0.0692 0.296  3744

#Heatmaps of subclonal nec, subclonal aca and clonal mutations.  

snvStatus = snvStatusReduce %>%
  pivot_wider(
    id_cols = Simple,
    names_from = Tumor_ID,
    values_from = c(Clonal, Subclonal)
  ) %>%
  mutate(
    Clonal = Clonal_Aca,
    Clonal_Aca = NULL,
         Clonal_Nec = NULL,
         ) %>%
  rowwise() %>%
  mutate(Total = sum(Subclonal_Aca + Subclonal_Nec + Clonal)) %>%
  mutate(
    Subclonal_Aca_Prop = Subclonal_Aca/Total,
    Subclonal_Nec_Prop = Subclonal_Nec/Total,
    Clonal_Prop = Clonal/Total
  )

# Order the samples.

PatientsOrders = snvStatus %>%
  select(Simple, Subclonal_Aca_Prop, Subclonal_Nec_Prop, Clonal_Prop) %>%
  arrange(desc(Clonal_Prop)) %>%
  pull(Simple)



  
p1 = snvStatus %>%
  select(Simple, Subclonal_Aca_Prop, Subclonal_Nec_Prop, Clonal_Prop) %>%
  setNames(nm = c("Simple","Subclonal_Aca","Subclonal_Nec","Clonal")) %>%
  reshape2::melt(
    id.vars = "Simple",
    variable.name = "MutType"
  ) %>%
  mutate(Simple = factor(Simple, levels = PatientsOrders)) %>%
  ggplot(aes(x = Simple, y = value  ,fill = MutType)) + theme_pubr() +
  geom_bar(stat = "identity", position = position_stack() ) +
  colorspace::scale_fill_discrete_diverging() +
  labs(x= NULL, y = "Proportion of Mutations") +
  theme(
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6 ),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

## The ITH of SNVs.

p2= snvStatusReduce %>%
  mutate(Simple = factor(Simple, levels = PatientsOrders)) %>%
  mutate(ITH = 1 - ITH) %>%
  ggplot(aes(x = Simple, y = Tumor_ID, fill = ITH)) + theme_pubr() +
  geom_tile(col = "white") +
  labs(x = NULL, y= "ITH-M") +
  scale_fill_distiller(palette = 'RdBu') +
  #scale_fill_gradient2( ) +
  theme(
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6 ),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_blank()
  ) 



#See the CNVs.

segStatus = segStatusReduce %>%
  pivot_wider(
    id_cols = Simple,
    names_from = Tumor_ID,
    values_from = c(Clonal, Subclonal)
  ) %>%
  mutate(
    Clonal = Clonal_Aca,
    Clonal_Aca = NULL,
    Clonal_Nec = NULL,
  ) %>%
  rowwise() %>%
  mutate(Total = sum(Subclonal_Aca + Subclonal_Nec + Clonal)) %>%
  mutate(
    Subclonal_Aca_Prop = Subclonal_Aca/Total,
    Subclonal_Nec_Prop = Subclonal_Nec/Total,
    Clonal_Prop = Clonal/Total
  )

#M16, M23 and M29
#segStatus$Subclonal_Aca_Prop[ segStatus$Simple %in% c("M23", "M16", "M29") ] = NA
#segStatus$Subclonal_Nec_Prop[ segStatus$Simple %in% c("M23", "M16", "M29") ] = NA
#segStatus$Clonal_Prop[ segStatus$Simple %in% c("M23", "M16", "M29") ] = NA


p3 = segStatus %>%
  select(Simple, Subclonal_Aca_Prop, Subclonal_Nec_Prop, Clonal_Prop) %>%
  mutate(Simple = factor(Simple, levels = PatientsOrders)) %>%
  setNames(nm = c("Simple","Subclonal_Aca","Subclonal_Nec","Clonal")) %>%
  reshape2::melt(
    id.vars = "Simple",
    variable.name = "MutType"
  ) %>%
  ggplot(aes(x = Simple, y = value  ,fill = MutType)) + theme_pubr() +
  geom_bar(stat = "identity", position = position_stack() ) +
  colorspace::scale_fill_discrete_diverging(na.value="grey") +
  labs(x= NULL, y = "Proportion of Neutral LOH") +
  theme(
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6 ),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

## The ITH of SNVs.
#segStatusReduce[segStatusReduce$Simple %in% c("M16", "M23", "M29", "M22"), "ITH" ] = NA

p4 = segStatusReduce %>%
  mutate(Simple = factor(Simple, levels = PatientsOrders)) %>%
  mutate(ITH = ifelse(ITH >= 0.34, 0.34, ITH)) %>%
  mutate(ITH = 1 - ITH ) %>%
  ggplot(aes(x = Simple, y = Tumor_ID, fill = ITH)) + theme_pubr() +
  geom_tile(col = "white") +
  labs(x = NULL, y= "ITH-L") +
  scale_fill_distiller(palette = 'RdBu', na.value="grey") +
  #scale_fill_gradient2( ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6 ),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) 

pdf("plot/ITHs_neutral-LOH_CNVs_SNVs_Heatmaps-3.pdf", width = 6.8, height = 6.8)
p1/p3/p2/p4 + plot_layout(heights = c(1, 1, 0.15, 0.15))
dev.off()


# See correlation between Divergence time from mutations and cnLOH.

ITH = left_join( snvStatusReduce %>%
  mutate(ITH = 1 - ITH) %>%
  pivot_wider(
    id_cols = Simple,
    names_from = Tumor_ID,
    names_prefix = "SNV_",
    values_from = ITH
  ),
segStatusReduce %>%
  mutate(ITH = 1 - ITH) %>%
  pivot_wider(
    id_cols = Simple,
    names_from = Tumor_ID,
    names_prefix = "CNV_",
    values_from = ITH
  )
) %>% as.data.frame() 
  
p1 = ITH %>%
  filter(!Simple %in% c("M29") ) %>%
  Jerry::ggpoints(x = "SNV_Aca", y = "CNV_Aca",
                  annotations = c("cor"),
                  linetypes = "lm.line",
                  xlab = "T index of ACA from mutations", 
                  ylab = "T index of ACA from cnLOH",
                  label.x.cor = 0.3,
                  label.y.cor = 0.3
                  ) +
  theme(
    axis.title = element_text(size = 14)
  )

p2 = ITH %>%
  filter(!Simple %in% c("M29") ) %>%
  Jerry::ggpoints(x = "SNV_Nec", y = "CNV_Nec",
                  annotations = c("cor"),
                  linetypes = "lm.line",
                  xlab = "T index of NEC from mutations", 
                  ylab = "T index of NEC from cnLOH",
                  label.x.cor = 0.3,
                  label.y.cor = 0.3
  ) +
  theme(
    axis.title = element_text(size = 14)
  )




pdf("plot/Tindex-SNV-CNV-components.pdf", width = 7.3, height = 3.8)
p1 + p2
dev.off()


p3 = ITH %>%
  Jerry::ggpoints(x = "SNV_Aca", y = "SNV_Nec",
                  annotations = c("cor"),
                  linetypes = "xy.line",
                  xlab = "T index of ACA", 
                  ylab = "T index of NEC",
                  label.x.cor = 0.1,
                  label.y.cor = 0.9
  ) +
  theme(
    axis.title = element_text(size = 14)
  ) +
  geom_linerange(x = 0.5, ymin = 0, ymax = 0.5, size = 1.2, linetype = 5, col = "#9b9bc5") +
  geom_linerange(y = 0.5, xmin = 0, xmax = 0.5, size = 1.2, linetype = 5, col = "#9b9bc5") 
  
  

# p value = 4.072e-14
  
pdf("plot/Tindex-SNV.pdf", width = 7.3/2, height = 3.8)
p3
dev.off()


#see cnLOH

ITH %>%
  filter(!Simple %in% c("M29") ) %>%
  pairwise_t_test(
    CNV_Aca ~ CNV_Nec
  )


ITH = ITH %>%
  filter(!Simple %in% c("M29"))

t.test(ITH$CNV_Aca, ITH$CNV_Nec, paired = T)

p4 = ITH %>%
  filter(!Simple %in% c("M16","M23","M29", "M22") ) %>%
  mutate( CNV_Aca = 1 - CNV_Aca,
          CNV_Nec = 1 - CNV_Nec
          ) %>%
  Jerry::ggpoints(x = "CNV_Aca", y = "CNV_Nec",
                  annotations = c("cor"),
                  linetypes = "xy.line",
                  xlab = "Prop of private cnLOH in ACA", 
                  ylab = "Prop of private cnLOH in NEC",
                  label.x.cor = 0.1,
                  label.y.cor = 0.9
  ) +
  theme(
    axis.title = element_text(size = 14)
  )


p5 = segStatusReduce %>%
  filter(!Simple %in% c("M29", "M16","M23")) %>%
  ggboxplot(x = "Tumor_ID", y = "Subclonal",
            color = "Tumor_ID",
            palette = "jco",
            add = "jitter") +
  stat_compare_means(method = "wilcox.test", paired = T) +
  labs(x = NULL, y = "Private cnLOH") +
  theme(
    legend.position = "none"
  )


pdf("plot/Private-cnLOH.pdf", width = 7.3/2, height = 3)
p5
dev.off()

ITH %>% filter(Simple %in% c("M03","M09","M12","M14","M15","M33")) %>%
  mutate(SNV = (SNV_Aca+SNV_Nec)/2 )


######################################################################################

#ITH scores.
ITH = left_join(
  snvStatusLong %>%
         setNames( nm = c("Simple","Aca_cnv","Nec_cnv") ) ,
   segStatusLong %>%
      setNames( nm = c("Simple","Aca_snv","Nec_snv") )
  )%>%
  left_join(
    
    inner_join( p.surv.Aca.nLOH$survData %>%
                  dplyr::select(Patient.ID, ITH) %>%
                  setNames(nm = c("Simple","ITH_Aca_cnv") ),
                
                p.surv.Nec.nLOH$survData %>%
                  dplyr::select(Patient.ID, ITH) %>%
                  setNames(nm = c("Simple","ITH_Nec_cnv") )
    ) 
    
  )


write.table(ITH, file = "datasets/ITH.text", row.names = F, quote = F)

#scores by evolution.

scores = read.table("../ITH/MRS_analysis_sequenza/MRS.scores.sequenza.txt", header = T)



ITH_score = left_join(
  ITH, scores
)
  
ITH_score %>%
  na.omit() %>%
  ggplot(aes(x = ITH_Nec_cnv , y = rAUC ) ) + theme_classic() +
  geom_boxplot() +
  geom_jitter(width = 0.3)


ITH_score %>%
  na.omit() %>%
  ggplot(aes(x = ITH_Aca_cnv , y = rAUC  ) ) + theme_classic() +
  geom_boxplot()+
  geom_jitter(width = 0.3)





