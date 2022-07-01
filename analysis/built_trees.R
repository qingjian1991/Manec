#########################################################################################################

# Built phylogenetic trees from somatic mutations and CNAs.

#Fig 5A: mutation-based trees.
#Fig S9: CNV-based trees.

#The tree is built by MPTevol (https://rpubs.com/cqj_00/MPTevol)

#########################################################################################################

library(MesKit)
library(tidyverse)
library(patchwork)
library(ggsci)

#############################################################################################################
#Read mafs. The mutations is from CHAT-corrected mutations.

load("../isma/maftools/driverGene.rda")

maf = readMaf(mafFile = sprintf("data/meskit/%s_snv.addInfo.maf", "Manec"),
              ccfFile = sprintf("data/meskit/%s_ccfinfo.maf", "Manec"),
              clinicalFile  = sprintf("data/meskit/%s_clinc.txt", "Manec"),
              refBuild = "hg19")

#############################################################

#load samples info
samples.info = readRDS("data/samples.info.rds")

###############################################################################
# We use the MPTevol developed by our group to construct the phylogenetic trees.
#
#
######


##################################################################################
# Mutation-based trees.
#
###
tree.list = list()

min.vaf = setNames(c(0.10, 0.2, 0.10, 0.12, 0.15, 0.15 ), nm = samples.info$multiple )

for(i in samples.info$multiple){
  message(i)
  
  # Set the Group Information
  group <- list(
    Aca = maf[[i]]@sample.info %>% filter(Tumor_ID == "Aca") %>% pull(Tumor_Sample_Barcode),
    Nec = maf[[i]]@sample.info %>% filter(Tumor_ID == "Nec") %>% pull(Tumor_Sample_Barcode)
  )
  
  # set group colors
  group.colors <- setNames(c("#C86166","#476BAE"), nm = c("Aca","Nec") )
  
  
  mutTrees <- plotMutTree(maf,
                          patient.id = i, 
                          group = group,
                          group.colors = group.colors,
                          min.vaf = min.vaf[i],
                          title = i
  )
  
  tree.list[[i]] <-  mutTrees$plot + 
    theme(legend.position = "none")
  
}


pdf("figures/Phylogenetic.mutation-based.trees.pdf", width = 12, height = 7)

wrap_plots(tree.list)

dev.off()



##################################################################################
# CNV-based trees.
#
###

#This distance file is built by MEDICC 
#see: Schwarz, R. F., Trinh, A., Sipos, B., Brenton, J. D., Goldman, N., & Markowetz, F. (2014). Phylogenetic quantification of intra-tumour heterogeneity. PLoS Comput Biol, 10(4), e1003535. doi:10.1371/journal.pcbi.1003535

tree.list = list()

for(i in samples.info$multiple){
  message(i)
  
  # Set the Group Information
  group <- list(
    Aca = maf[[i]]@sample.info %>% filter(Tumor_ID == "Aca") %>% pull(Tumor_Sample_Barcode),
    Nec = maf[[i]]@sample.info %>% filter(Tumor_ID == "Nec") %>% pull(Tumor_Sample_Barcode)
  )
  
  # set group colors
  group.colors <- setNames(c("#C86166","#476BAE"), nm = c("Aca","Nec") )
  
  
  # built trees
  cnaTree <- plotCNAtree(
    dist =  sprintf("data/%s_tree_final.dist", i) ,
    group = group,
    group.colors = group.colors,
    title = i
  )
  
  tree.list[[i]] <-  cnaTree$plot + 
    theme(legend.position = "none")
  
}


pdf("figures/Phylogenetic.CNA-based.trees.pdf", width = 12, height = 7)

wrap_plots(tree.list)

dev.off()


