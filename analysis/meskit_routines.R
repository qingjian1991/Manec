# We want to know how to estimated the divergence routines from the pair-wise comparison 
# We use the pyclone results as the input.

library(MesKit)
library(tidyverse)
library(patchwork)
library(ggpubr)


###########################################################################################
#For MRS data, we first calculate the JSI and T index. We use the two components from the same sampling site to get their index values.


# Fig 5B-G
# Fig 5B: applied JSI index to estimate the divergence levels between ACA and NEC
# Fig 5C-D: The paired JSI in patients with MSR (multi-region sequencing)
# Fig 5E: The relationships between JSI and divergence time (T index)
# Fig 5F: Distributions of JSI in 27 SSR patients.
# Fig 5G: the neutral evolutionary index.


###########################################################################################

maf.split = readMaf(mafFile = sprintf("data/meskit/%s_pyclone.snv.addInfo.maf", "Manec"),
              ccfFile = sprintf("data/meskit/%s_pyclone.ccfinfo.maf", "Manec"),
              clinicalFile  = sprintf("data/meskit/%s_pyclone.clinic.split.maf", "Manec"),
              remove.empty.VAF = FALSE,
              refBuild = "hg19")


#choose the references samples to calculate the index.
patient.id.MRS = setNames( c("Aca_1C","Aca_10A","Aca_1B","Aca_1C","Aca_1A","Aca_9E"),
                           nm = c("M03","M09","M12","M14","M15","M33")
                           )

pdf("meskit/CCF.paired.pyclone.MRS.pdf" , width = 8.7, height = 5.6)

for(i in 1:6){
  Hindex = calRoutines(maf = maf.split,
                     patient.id = names(patient.id.MRS[i]),
                     PrimaryId = patient.id.MRS[i],
                     CCF_cutoff = 0.1,
                     pairByTumor = TRUE
  )
  
  Hindex[[names(patient.id.MRS[i])]]$plist %>%
    wrap_plots(nrow =2) %>% plot()
  
}

dev.off()

# Get the scores.

HIscores = list()

for(i in 1:6){

  patient.id = names(patient.id.MRS[i])
  Tumor_ID = maf.split[[patient.id]]@sample.info$Tumor_ID

  for(j in Tumor_ID){
    Hindex = calRoutines(maf = maf.split,
                       patient.id = patient.id,
                       PrimaryId = j,
                       CCF_cutoff = 0.1,
                       pairByTumor = TRUE
    )
    
    names = names(Hindex[[patient.id]]$stats)
      
    HIscores[[patient.id]][[j]] = purrr::reduce(Hindex[[patient.id]]$stats, rbind) %>%
      as.data.frame() %>%
      mutate(
        names = names,
        patient.id = patient.id
      ) 
  }
  
}


scores = lapply( samples.info$multiple , function(x) 
  purrr::reduce(HIscores[[x]], rbind) %>% 
    as.data.frame() %>%
    mutate(
      Tindex = (Lp +Lm)/(Lp + Lm + La ),
      Tp = La/(Lp + La),
      Tm = La/(Lm + La),
      PrimaryID = mapply(function(x) x[1], str_split(names, "-")),
      TumorID = mapply(function(x) x[2], str_split(names, "-"))
    )  ) 

names(scores) = samples.info$multiple 




#' getJSI.component
#' Get the mean of JSI values according to the components.

getJSI.component = function(JSI.pair, components = c("Aca","Nec") ){
  
  diag(JSI.pair) <- NA
  
  comp.mean = setNames( rep(NA, length(components)) , nm = components)
  
  # get Means of same components.
  for(i in components){
    sample.names = grepl(i , rownames(JSI.pair))
    
    comp.mean[i] = JSI.pair[sample.names,sample.names] %>%
      mean(na.rm = TRUE)
  }
  
  
  comp.mean["Between"] = JSI.pair[grepl(components[1] , rownames(JSI.pair)), !grepl(components[1] , rownames(JSI.pair))] %>%
    mean(na.rm = TRUE)
  
  
  return(comp.mean)
}


#plot Pairs.
plot.JSI = function(patient.id, JSI.pair, col.max = "#BB3058"){
  
  message(patient.id)
  
  colums.orders =  sort(rownames(JSI.pair))
  
  JSI.mean = getJSI.component(JSI.pair)
  
  MesKit:::plotCorr(
    JSI.pair[colums.orders, colums.orders],
    number.col = "black",
    use.circle = TRUE,
    title  = sprintf( "%s\t Aca: %.2f\t Nec: %.2f\t Between: %.2f", patient.id, JSI.mean["Aca"], JSI.mean["Nec"], JSI.mean["Between"] )
  )
  
}


#Plot Paired Tindex

scores.Tindex.pair = lapply( scores , function(x)
  x %>%
    reshape2::dcast(
      PrimaryID ~ TumorID,
      value.var = "Tindex",
      fill = 1
    ) %>%
    column_to_rownames("PrimaryID") %>%
    as.matrix()
  
)

pdf("figures/CCF.paired.pyclone.Tindex.MRS.pdf" , width = 7, height = 5.6)

  for(i in names(scores.Tindex.pair)){
    plot.JSI(i, scores.Tindex.pair[[i]], col.max = "#BB3058") %>% print()
  }

dev.off()

#Plot Paired JSI index 

scores.JSI.pair = lapply( scores , function(x){
  mat  = x %>%
    reshape2::dcast(
      PrimaryID ~ TumorID,
      value.var = "JSI",
      fill = 1
    ) %>%
    column_to_rownames("PrimaryID") %>%
    as.matrix()
  
  patient.id = unique(x$patient.id)
  
  sampleIDs = mapply(function (x) str_c(patient.id, x[2], x[1], sep = "_"), str_split(rownames(mat), "_") )
  
  rownames(mat) = sampleIDs
  colnames(mat) = sampleIDs
  
  mat

  }
)

pdf("figures/CCF.paired.pyclone.JSI.MRS.pdf" , width = 7, height = 5.6)

for(i in names(scores.JSI.pair)){
  plot.JSI(i, scores.JSI.pair[[i]], col.max = "#00A77E") %>% print()
}

dev.off()

### Get the summary of scores

# see scores between Tindex and JSI index
scores.MRS.split = purrr::reduce(scores, rbind) %>%
  mutate(comp1 = mapply(function(x) x[1], str_split(PrimaryID, "_") ),
         comp2 = mapply(function(x) x[1], str_split(TumorID, "_") ),
         site1 = mapply(function(x) x[2], str_split(PrimaryID, "_") ),
         site2 = mapply(function(x) x[2], str_split(TumorID, "_") )
         ) 

p1 = scores.MRS.split %>%
  filter(comp1 != comp2 & site1 == site2) %>%
  mutate(Tindex = 1 -Tindex) %>%
  Jerry::ggpoints(
    x = "Tindex", y = "JSI",
    annotations = c("cor"),
    linetypes = c("lm.line"),
    xlab = "Tindex", ylab = "JSI",
    label.x.cor =  0.2, 
    label.y.cor = 1.1
  ) +
  geom_smooth(method = "lm", linetype =2)

p1

pdf("figures/CCF.JSI_Time.points.pdf", width = 3.7, height = 3.7)
p1
dev.off()



#Get the mean values of 

# see scores between Tindex and JSI index
scores.MRS = scores.MRS.split %>%
  filter(comp1 != comp2 & site1 == site2) %>%
  filter(comp1 == 'Aca') %>%
  group_by(patient.id) %>%
  summarise(
    Lp = round(mean(Lp)),
    Lm = round(mean(Lm)),
    La = round(mean(La)),
    JSI = mean(JSI),
    Tindex = mean(Tindex)
  ) 

################################################################################################

# Similar with the metastatic routines inference where paired CCFs between primary and metastatic tumors were calculated, we calculate the JSI and H index to indicate the early(later) and single(multiple) clonal changes between ACA and NEC.

maf = readMaf(mafFile = sprintf("data/meskit/%s_pyclone.snv.addInfo.maf", "Manec"),
                    ccfFile = sprintf("data/meskit/%s_pyclone.ccfinfo.maf", "Manec"),
                    clinicalFile  = sprintf("data/meskit/%s_clinc.txt", "Manec"),
                    remove.empty.VAF = FALSE,
                    refBuild = "hg19")

Hindex = calRoutines(maf = maf,
                   patient.id = NULL,
                   PrimaryId = "Aca",
                   CCF_cutoff = 0.1,
                   pairByTumor = TRUE,
                   subtitle = "JSI"
)


#Plot every patient in single page.
pdf("figures/CCF.Paired.JSI.SRS.singlePage.pdf", width = 4, height = 3.5)
  lapply(Hindex, function(x) x$plist)
dev.off()

# combined different patients into one page.

plist.CCF = lapply(Hindex[samples.info$single[!samples.info$single %in% c("M16","M23","M29") ]], function(x) x$plist[[1]] + labs(x = NULL, y = NULL)
                     )

pdf( "figures/CCF.Paired.JSI.SRS.pdf", width = 12, height = 7.5)
  wrap_plots( plist.CCF, ncol = 6)
dev.off()




#calculate different scores.

scores.all = lapply(Hindex, function(x) x$stats[[1]]) %>% 
  purrr::reduce(rbind) %>%
  as.data.frame() %>%
  select(Lp, Lm, La, JSI) %>%
  mutate(patient.id = names(Hindex)) %>%
  filter(patient.id %in% samples.info$single) %>%
  rbind(
    scores.MRS %>% select(Lp, Lm, La, JSI, patient.id)
  ) %>%
  mutate(
    Tindex = (Lp +Lm)/(Lp + Lm + La ),
    Tp = Lp/(Lp + La ),
    Tm = Lm/(Lm + La )
  ) %>%
  left_join(
    clinical %>% select(Family, Met_Site, Tumor_Size_cm, Simple, Sex, Age, Adeno, Death, Survival_days,Survival_months, TNM_levels, Adeno, Neuro, Neuro_ki67 ) %>%
      unique.data.frame(),
    by = c("patient.id" = "Simple")
  )  %>%
  mutate(
    TNM_levels = str_remove(TNM_levels, "A|B|C"),
    Family = factor(Family),
    Neuro = as.numeric(str_remove(Neuro, "%")),
    Neuro_ki67 = as.numeric(str_remove(Neuro_ki67, "%")),
    
  ) %>%
  mutate(  
    JSI.group = ifelse(JSI >= 0.3, "High", "Low"),
    JSI.group = ifelse(patient.id %in% c("M03", "M14", "M15"), "Low", JSI.group),
    JSI.group1 = ifelse(patient.id %in% c("M02","M04","M08","M11","M13","M24", "M25", "M09", "M12","M33"), "High","low")
  )


source("analysis/plot_survival.R")

plot.list = plot_survival(
  data = scores.all %>%
    filter(patient.id %in% samples.info$single)  %>%
    mutate(clonality = ifelse(JSI >0.2, "Mono", "Poly") ) ,
  time = "Survival_months",
  status = "Death",
  variables = c("clonality"),
  cutoff = "best",
  minprop = 0.2
)

plot.list$clonality

library(ezcox)

zz = ezcox(
  scores.all %>%
    filter(patient.id %in% samples.info$single)  %>%
    mutate(clonality = ifelse(JSI >0.2, "Mono", "Poly") ) ,
  covariates = c("clonality"),
  controls = NULL,
  time = "Survival_months",
  status = "Death",
  return_models = TRUE
)

mds = get_models(zz)

show_models(mds, drop_controls = TRUE, merge_models = T)

pdf("figures/CCF.Clonality.Survival.pdf", width = 6, height = 6.5)

plot.list$clonality

dev.off()


pdf("figures/CCF.Clonality.Cox.pdf", width = 8, height = 3)

show_models(mds, drop_controls = TRUE, merge_models = T)

dev.off()


#using boxplot to show the relationship between Tp and Tm
scores %>%
  mutate(Simple = names(Hindex) ) %>%
  filter(!Simple %in% samples.info$multiple ) %>%
  reshape2::melt(measure.vars = c("Tp","Tm")) %>%
  ggplot( aes(x = variable, y = value) ) + theme_classic() +
  geom_boxplot()


#################################################################################################

#The following code is copied from Proportional_SCNA_SNV.R

load("data/Proportional.rda")

segStatusReduce = readRDS("data/segStatusReduce.rds")

# see the relationship between mutations and LOH.


library(ggpmisc)

estimate_type = "nLOH"

#Only nLOH is significantly.
#total: contains different kinds of chromosomal changes: Amp, Del and LOH
#estimate_type = "Del"

segStatusReduce  = segStatus_wider %>%
  filter(Tumor_ID != "Squca") %>%
  group_by(Tumor_ID, Simple) %>%
  summarise(Clonal = mean( get( paste("Clonal", estimate_type, sep = "_") ) ),
            Subclonal = mean( get(  paste("Subclonal", estimate_type, sep = "_") ) ),
            ITH = mean( get(  paste("ITH", estimate_type, sep = "_") ) ),
            GII = mean( get(  paste("GII", estimate_type, sep = "_") ) )
  )

segStatusLong = segStatusReduce %>%
  #filter(!Simple %in% c("M16","M23","M29", "M22") ) %>%
  reshape2::dcast(
    Simple ~ Tumor_ID,
    value.var = "ITH"
  )

scores.combined = segStatusReduce %>%
  pivot_wider(
    id_cols = Simple,
    names_from = Tumor_ID,
    values_from = c(Clonal, Subclonal, ITH)
  ) %>%
  right_join(
    scores.all,
    by = c("Simple" = "patient.id")
  ) %>%
  left_join(
    snvStatusReduce %>%
      pivot_wider(
        id_cols = Simple,
        names_prefix = "snv",
        names_from = Tumor_ID,
        values_from = c(Clonal, Subclonal, ITH, TMB)
      )
  ) %>%
  mutate(
    Clonal_Nec = NULL, Clonal_snvNec = NULL
  ) %>%
  dplyr::rename(
    Shared = Clonal_Aca,
    Private_Aca = Subclonal_Aca,
    Private_Nec = Subclonal_Nec
  ) %>%
  as.data.frame()
 


# See survivals 

plot.list = plot_survival(
  data = scores.combined %>%
    filter(Simple != "M29") ,
  time = "Survival_months",
  status = "Death",
  variables = c("Private_Aca","Private_Nec", "Shared"),
  legend.labs = c("High","Low"),
  cutoff = "best",
  minprop = 0.4
)

pdf("figures/CCF.cnLOH.Survival.pdf", width = 4.8, height = 5.3)

plot.list

dev.off()



zz = ezcox(
  scores.combined %>%
    filter(Simple != "M29"),
  covariates =c("Private_Aca","Private_Nec", "Shared"),
  controls = c("Sex","Age","Family"),
  time = "Survival_months",
  status = "Death",
  return_models = TRUE
)

mds = get_models(zz)


show_models(mds, drop_controls = TRUE, merge_models = T)



# The JSI distributions.

p1 = scores.combined %>%
  filter(Simple %in% samples.info$single) %>%
  ggplot( aes(x = JSI) ) + theme_classic() +
  geom_histogram(bins = 10) +
  geom_vline(xintercept = 0.25, linetype=2, size = 1.1 ) +
  labs(y = "# of samples") +
  theme(
    axis.title = element_text(size = 14)
  )

pdf("figures/CCF.JSI.distribution.pdf", width = 4, height = 3)
p1
dev.off()




JSI.value = scores.combined %>%
  filter(Simple %in% samples.info$single) %>%
  pull(JSI) 

table(JSI.value >= 0.2)
#FALSE  TRUE 
#20     7

7/27

write.table(scores.combined, file = "meskit/scores.combined.txt", sep = "\t", row.names = F)


################################################################################
# see neutral index in the two groups.

neutral.test = readRDS("data/neutral.test.rds")

neutral.test %>%
  filter(Simple %in% samples.info$single) %>%
  ggplot( aes(x = Tumor_ID, y = rAUC) ) + theme_classic() +
  geom_boxplot()

neutral.test %>%
  filter(Simple %in% samples.info$single) %>%
  ggboxplot(
    x = "Tumor_ID", y = "Rsq",
    add = "jitter"
  ) +
  stat_compare_means(method = "wilcox")
  

p1 = neutral.test %>%
  filter(Simple %in% samples.info$single) %>%
  filter( !sample %in% c("M34_1G_Aca","M34_1E_Nec") ) %>%
  left_join(
    scores.combined %>%
      select(Simple, JSI.group, JSI.group1, JSI)
  )%>%
  mutate(JSI.group2 = ifelse(JSI >= 0.2, "High", "Low") ) %>%
  ggboxplot(
    x = "JSI.group2", y = "Rsq",
    add = "jitter",
    facet.by = "Tumor_ID",
    col = "JSI.group2"
  ) +
  stat_compare_means(method = "wilcox", label.y = 1) +
  ylim(0.88, NA) +
  geom_hline(yintercept = 0.98, linetype = 2, col = "blue4", size = 1.1 ) +
  #stat_pvalue_manual()
  ggsci::scale_color_d3() +
  labs(y = latex2exp::TeX("R^2 Model Fit"), x = NULL) +
  theme(
    legend.position = "none"
  ) 

pdf("figures//CCF.JSI.NeutralTest.R2.pdf", width = 6.5, height = 4)
p1
dev.off()

