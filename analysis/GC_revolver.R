library(revolver)
library(ctree)

################################################################################
# Read Data
################################################################################


GC = readRDS("data/GC.revolve.PubPriv.all.rds")

# > head(GC)
#                                 Misc patientID variantID cluster is.driver is.clonal
# 1      MIB2:chr1:1562773:1562773:T:G       M01      MIB2       1     FALSE     FALSE
# 2     MMEL1:chr1:2525940:2525940:T:C       M01     MMEL1       3     FALSE     FALSE
# 3    CLSTN1:chr1:9804420:9804420:G:A       M01    CLSTN1       1     FALSE     FALSE
# 4 SPATA21:chr1:16751610:16751610:T:C       M01   SPATA21       4     FALSE      TRUE
# 5   PADI4:chr1:17660989:17660989:A:G       M01     PADI4       4     FALSE      TRUE
# 6    ECE1:chr1:21585151:21585151:G:A       M01      ECE1       3     FALSE     FALSE
#                                 CCF
# 1 M01_1D_Aca:0.660;M01_1D_Nec:0.000;
# 2 M01_1D_Aca:0.000;M01_1D_Nec:0.390;
# 3 M01_1D_Aca:0.530;M01_1D_Nec:0.000;
# 4 M01_1D_Aca:1.000;M01_1D_Nec:0.950;
# 5 M01_1D_Aca:1.000;M01_1D_Nec:1.000;
# 6 M01_1D_Aca:0.000;M01_1D_Nec:0.350;


subtype = "mtree"


#******* Quality controls:

#adjust the clonal ids.
#(1) Merge C1 and C2 for M14
#(2) Merge C1 and C2 for M33
#(3) Removing C2 for M02, too many mutations.


GC = GC %>%
  filter( !( cluster == 2 &  patientID == "M02") ) %>%
  mutate( cluster = ifelse(  ( cluster == 2 &  patientID == "M14") , 1,cluster),
          is.clonal = ifelse(  ( cluster == 1 &  patientID == "M14") , TRUE , is.clonal),
          cluster = ifelse(  ( cluster == 2 &  patientID == "M33") , 1,cluster),
          is.clonal = ifelse(  ( cluster == 1 &  patientID == "M33") , TRUE , is.clonal),
          )



#************ Group some genes family into one label to increase the statistical power.

#drivers = c(
#            "FAT4","FAT3","FAT2","FAT1", 
#            "APC", "APC2", 
#            "LRP1B","LRP1","LRP2",
#            "KMT2D","KMT2A","KMT2B","KMT2C",,
#            "REBB2","ERBB3","ERBB4",
#)


GC %>% as_tibble()
#f-FAT -> FAT1, FAT2, FAT3, FAT4
#f-APC -> APC, APC2
#f-LRP -> LRP1B, "LRP1","LRP2"
#f-KMT2 ->KMT2D, KMT2A, KMT2B, KMT2C,
#f-REBB ->REBB2,REBB3,REBB4

GC = GC %>%
  mutate( variantID = ifelse(variantID %in% c("FAT2","FAT4"), "f-FAT", variantID),
          variantID = ifelse(variantID %in% c("KMT2D","KMT2C"), "f-KMT2", variantID),
          variantID = ifelse(variantID %in% c("ERBB2","ERBB3","ERBB4"), "f-ERBB", variantID),
          variantID = ifelse(variantID %in% c("NOTCH2","NOTCH1"), "f-NOTCH", variantID)
  )

#************ Filtering for duplicated genes

#split data

#get the overlapping drivers.
GC_adj_nest = GC %>% 
  group_by(patientID, variantID) %>%
  mutate(mutNum = n()) %>%
  filter(
    mutNum >1 
  ) %>%
  mutate(mutNum = NULL) %>%
  nest()

#work within list

#unique drivers
mod_fun = function(df){
  
  num.driver = sum(df$is.driver)
  
  if( num.driver >= 2   ){
    #more than one drivers, removing all non-drivers.
    message("more than 2 drivers")
    df = df %>% filter(is.driver) 
  }else if(num.driver == 1){
    #select the drivers.
    message("only 1 drivers")
    df = df %>% filter(is.driver)
    return(df)
  }
  
  num.clonal = sum(df$is.clonal)
  
  df %>%
    mutate(meanCCF = mapply(function(x) mean( na.omit(as.numeric(x) ) ) + rnorm(1, 0.0001)  , str_split(CCF, "[:|;]") )) %>%
    suppressWarnings() %>%
    filter(meanCCF == max(meanCCF)) %>%
    mutate(mutNum = NULL,
           meanCCF = NULL
           )
}


GC_adj_m = GC_adj_nest %>%
  mutate(model = map(data, mod_fun) ) %>%
  select(patientID, variantID, model ) %>%
  unnest()


#combined overlap and non-overlap
GC_new = GC %>% 
  group_by(patientID, variantID) %>%
  mutate(mutNum = n()) %>%
  filter(
    mutNum ==1 
  ) %>%
  rbind(GC_adj_m) %>%
  mutate(mutNum = NULL)

#filtering patients without any drivers.
GC_new = GC_new %>%
  filter( !patientID %in% c("M29", "M30"))

# Constructor
GC_revolver= revolver_cohort(
  GC_new, 
  MIN.CLUSTER.SIZE = 0, 
  annotation = "GC")

#GC_revolver =remove_patients(GC_revolver, patientID = "M02")

#check the cohort
revolver_check_cohort(GC_revolver)


# Driver events that occur in 1 patient
non_recurrent = Stats_drivers(GC_revolver) %>% 
  filter(N_tot <=2) %>% 
  pull(variantID)

# Remove drivers
if( !is.null(length(non_recurrent)) ){
  GC_revolver = remove_drivers(GC_revolver, non_recurrent)
}

revolver:::get_duplicates(GC_revolver)


################################################################################
# Constructing mutation trees
################################################################################

GC_revolver = compute_mutation_trees(GC_revolver, overwrite = TRUE, parallel = FALSE )

#or Constructing clone trees
#GC_revolver = compute_clone_trees(GC_revolver, overwrite = TRUE )


#Fitting models with REVOLVER
#Function revolver_fit implements the 2-steps REVOLVER algorithm to fit the data.
#
#We use the following parameters:
#  
#  initial.solution = NA, to sample random initial solutions for every run of EM;
#n = 3, to repeat the fit 3 times, and retain the one with lower median goodness-of-fit penalty.
#parallel = FALSE, to run serially the fits;


GC_revolver = revolver_fit(
  GC_revolver, 
  parallel = T, 
  n = 10, 
  initial.solution = NA)


#Computing REVOLVER hierarchical clusters
GC_revolver = revolver_cluster(
  GC_revolver, 
  split.method = 'cutreeHybrid',
  min.group.size = 5
  )

saveRDS(GC_revolver, file = sprintf("data/GC_revolver_%s.PubPriv.rm.jack.rds", subtype) )


#evaluate different measures confidence for the clusters using jackknife statistics (revolver_jackknife);;
#Undo for most of time.
#GC_revolver = revolver_jackknife(GC_revolver) 
#saveRDS(GC_revolver, file = sprintf("plot_%s/GC_revolver.PubPriv.rds", subtype) )



#plot

# Plot the heatmaps of REVOLVER"s clusters, as tiles.
pdf( sprintf("figures/GC_%s_driver_clusters.pdf", subtype), width = 8, height = 8)
  plot_clusters(GC_revolver, cutoff_trajectories = 2, cutoff_drivers = 3)
dev.off()

pdf( sprintf("figures/GC_%s_drivers_graph.pdf", subtype), width = 6, height = 6)
  plot_drivers_graph(GC_revolver, min.occurrences = 3)
  plot_drivers_graph(GC_revolver, min.occurrences = 2)
dev.off()


pdf( sprintf("figures/GC_%s_dendrogram.pdf", subtype), width = 10, height = 3.2)
plot_dendrogram(GC_revolver)
dev.off()



plot_DET_index(GC_revolver)

plot(GC_revolver)


p1 = plot_clusters(GC_revolver, cutoff_trajectories = 2, cutoff_drivers = 3)
p2 = plot_dendrogram(GC_revolver)

p1 = p1 + labs(title = NULL)

p2/p1


pdf( sprintf("figures/GC_%s_drivers_clonality.pdf", subtype) , width = 6, height = 6)
plot_drivers_clonality(GC_revolver)
dev.off()


# Plot the index of Divergent Evolutionary Trajectories
pdf( sprintf("figures/GC_%s_drivers_penalty.pdf", subtype) , width = 7, height = 7)
plot_penalty( GC_revolver, alpha_level = 0.5, min.occurrences = 1 )
dev.off()

plot_drivers_occurrence(GC_revolver)

# Plot REVOLVER trees for a patient.
pdf(  sprintf("figures/GC_%s_patient_trees.pdf", subtype), width = 8, height = 6)

for(i in GC_revolver$patients){
  message(i)
  plot_patient_trees(GC_revolver, i) %>% print()
}

dev.off()



# Plot the data histogram for a patient.
plot_patient_CCF_histogram(GC_revolver, "M33")

# Plot the oncoprint for a patient.
plot_patient_oncoprint(GC_revolver, 'M20')

# Plot a number of different measurements for the patient
pdf( sprintf("figures/GC_%s_patient_data.pdf", subtype), width = 8, height = 6)

for(i in GC_revolver$patients[!GC_revolver$patients %in% c("M03","M09","M12","M14","M15","M33") ] ){
  message(i)
  plot_patient_data(GC_revolver, i) %>% print()
}

dev.off()

pdf( sprintf("figures/GC_%s_patient_data_MRS.pdf", subtype), width = 8*2, height = 6*2)

for(i in c("M03","M09","M12","M14","M15","M33") ){
  message(i)
  plot_patient_data(GC_revolver, i) %>% print()
}

dev.off()


# Plot the trajectories detected in at least 30% patients per cluster
pdf( sprintf("figures/GC_%s_trajectories_per_cluster.pdf", subtype), width = 12, height = 8)
plot_trajectories_per_cluster(GC_revolver, min_counts = 2) 
dev.off()



# Plot the patients' jackknife co-clustering probability

#pdf( sprintf("figures/GC_%s_jackknife_coclustering.pdf", subtype), width = 8, height = 8)
#plot_jackknife_coclustering(GC_revolver)
#dev.off()

# Plot the patients' jackknife cluster staability
#pdf( sprintf("figures/GC_%s_jackknife_cluster_stability.pdf", subtype), width = 4, height = 4)
#plot_jackknife_cluster_stability(GC_revolver)
#dev.off()


GC_cluster = Cluster(GC_revolver)

# This returns patient-level statistics like the number of biopsies, overall mutations, drivers,
# clones with drivers, truncal and subclonal mutations.
StatsCohort = Stats_cohort(GC_revolver)


# This returns driver-level statistics like the number of times the driver is clonal,
# subclonal, or found in general, and for quantity normalized by cohort size (i.e., the percentage)
driversStats =  Stats_drivers(GC_revolver)


# This returns patient-level statistics for the trees available in a patient. The tibble reports
# whether the patient has trees annotated, the total number of trees, their minimum and maximum
# scores mutations and the total number of differnet combinations of Information Transfer for 
# the available trees.
treesStats =  Stats_trees(GC_revolver)


# This returns the same table of above, but with some extended information on the fits (like the fit rank, etc)
fitsStats =  Stats_fits(GC_revolver)

