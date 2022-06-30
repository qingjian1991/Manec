#Survival for manec.


# Survival analysis -------------------------------------------------------

library(survival)
library(survminer)
library(patchwork)

subtype = "mtree"

clinical = read.table("data/clinical.20210508.txt", sep = "\t", header = T) %>%
  select(Simple, Survival_days, Death, TNM_levels ) %>%
  unique.data.frame() %>%
  mutate(OS.month = Survival_days/30,
         OS = Death,
         TNM = str_remove(TNM_levels, "[A|B|C]")
  )%>%
  mutate(OS.month = ifelse(OS.month>=84, 84, OS.month) )

#Cluster info
GC_revolver= readRDS( sprintf("data/GC_revolver_%s.PubPriv.rm.jack.rds", subtype) )

GC_cluster = Cluster(GC_revolver) 

GC_cluster = GC_cluster %>%
  left_join(clinical, by = c("patientID" = "Simple")) %>%
  #filter(cluster != "C0") %>%
  mutate(
    cluster_merge = ifelse(!cluster %in% c("C1", "C3"), "C2", cluster )
  )

#GC_cluster = GC_cluster %>%
#  mutate(cluster = ifelse(cluster == "C1", "C3", cluster) ) %>%
#  mutate(cluster = ifelse(cluster == "C2", "C4", cluster) ) %>%
#  left_join(clinical, by = c("patientID" = "Simple"))


# Survival for all groups
fit = survfit(Surv(GC_cluster$OS.month, GC_cluster$OS) ~ GC_cluster$cluster_merge)

groups =  c("C1", "C2", "C3", "C4", "C5", "C6")
colors = c("#E41A1C","#E7B800", "#984EA3", "#377EB8", "#FF7F00" , "#4DAF4A")
names(colors) = groups
colors = colors[1: length( unique(GC_cluster$cluster_merge)) ]

pdf( sprintf("figures/GC_%s_survival_C2.pdf", subtype) , width = 4, height = 4.5)
p1 = ggsurvplot(fit,
           censor = T,
           data = GC_cluster,
           risk.table = T, 
           risk.table.height = 0.3, 
           linetype = 1, 
           pval = TRUE,
           palette = colors,
           legend.labs = names(colors),
           ncensor.plot = FALSE)
p1$plot + p1$table +
  plot_layout(heights = c(1, 0.3))

dev.off()


##############################################################

# If remove M01 and M32 in C1, the p-value is significantly.


GC_survival = GC_cluster %>%
  mutate(cluster_merge1 = ifelse(patientID %in% c("M02","M16","M09", "M01","M32") ,"C3", cluster_merge) )


GC_survival = GC_cluster %>%
  mutate(cluster_merge1 = cluster_merge) %>%
  filter(!patientID %in% c("M01", "M32"))


# Survival for all groups
fit = survfit(Surv(GC_survival$OS.month, GC_survival$OS) ~ GC_survival$cluster_merge1)

groups =  c("C1", "C2", "C3", "C4", "C5", "C6")
colors = c("#E41A1C","#E7B800", "#984EA3", "#377EB8", "#FF7F00" , "#4DAF4A")
names(colors) = groups
colors = colors[1: length( unique(GC_survival$cluster_merge1)) ]


p1 = ggsurvplot(fit,
                censor = T,
                data = GC_survival,
                risk.table = T, 
                risk.table.height = 0.3, 
                linetype = 1, 
                pval = TRUE,
                palette = colors,
                legend.labs = names(colors),
                ncensor.plot = FALSE)

pdf( sprintf("figures/GC_%s_survival_rm2.pdf", subtype) , width = 4.5, height = 5)

p1$plot + p1$table +
  plot_layout(heights = c(1, 0.3))

dev.off()

