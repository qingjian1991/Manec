## The Main analysis

### 1. intratumor heterogeneity (ITH)
**Scripts**:  

`Proportional_SCNA_SNV.R`

Figure 4 A-C: Estimate intratumor heterogeneity(ITH) to infer the tumor origins (shared common ancestor) and tumor divergence(early vs. late divergence using T index).


### 2. JSI index

The level of divergence between ACA and NEC is estimated by the Jaccard similarity index of their mutations, similar to the level of divergence between primary and metastatic tumors. The JSI level is calculated using the function `calRoutines` in `MPTevol`(https://github.com/qingjian1991/MPTevol)


**Scripts**: 
`meskit_routines.R` Fig 5B-G

Fig 5B: applied JSI index to estimate the divergence levels between ACA and NEC

Fig 5C-D: The paired JSI in patients with MSR (multi-region sequencing)

Fig 5E: The relationships between JSI and divergence time (T index)

Fig 5F: Distributions of JSI in 27 SSR patients.

Fig 5G: the neutral evolutionary index.


### 3. Whole Genome Doubling (WGD)

The WGD status of each sample is determined by `EstimateClonality` (https://github.com/qingjian1991/EstimateClonality).

**Scripts**: 

`MCN_analysis.R`

Fig 3B: The scatter plot of tumor ploidy and the fraction of autosomal tumor genomes with MCN >=2

Fig 3D: The timing of WGD with genome loss.


`Segs_LOH.R`

Fig 3C and Supplementary Figures:

Make comparisons between MANEC-ACA, MANEC-NEC, STAD-CIN:GD and STAD-CIN:noGD. including:

(1) proportional of copy neutral-LOH, haploid LOH (Fig 3C)

(2) clonal and subclonal LOH (Fig S4C)

(3) MATH scores (tumor heterogeneity scores) (Fig S4D)


`Segs_essential.R`

Fig 3F: calculate the relationship between the proportions of copy number change before (MCN = 1) and after WGD (MCN>=2) with the gene density (including essential genes, tumor suppersor genes and oncogenes ) in each arm.

### 4. Estimate the driver order.

The temporal orders of driver genes are determined by `revolver` (https://github.com/caravagnalab/revolver).

**Scripts**: 
`GC_revolver.R `

Fig 6B-D


`GC_revolver_surivival.R`

Supplementary Fig S13


### 5. Built phylogenetic trees

Constructing of a phylogenetic tree based on somatic mutations or copy number alterations using MPTevol (https://rpubs.com/cqj_00/MPTevol)

**Scripts**: 
`built_trees.R`


1. Fig 5A: mutation-based trees.

2. Fig S9: CNV-based trees.
