
library(reshape2)


load("data/Proportional.rda")


#########################################################################################################

# Figure 4 A-C:

#Estimate intratumor heterogeneity(ITH) to infer the tumor origins (shared common ancestor) and tumor divergence(early vs.late T index).


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

pdf("figures//ITH_cnLOH_SNVs_Heatmaps.pdf", width = 6.8, height = 6.8)
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




pdf("figures/ITH.Tindex.SNV.CNV.pdf", width = 7.3, height = 3.8)
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
  
pdf("figures/ITH.Tindex.SNV.pdf", width = 7.3/2, height = 3.8)
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


pdf("figures/ITH.Private-cnLOH.boxplot.pdf", width = 7.3/2, height = 3)
p5
dev.off()

ITH %>% filter(Simple %in% c("M03","M09","M12","M14","M15","M33")) %>%
  mutate(SNV = (SNV_Aca+SNV_Nec)/2 )


