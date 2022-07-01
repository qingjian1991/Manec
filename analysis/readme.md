## The Main analysis
In this fold, the main analysis are including: 

(1) Estimating the CCF values based on CHAT

(2) The comparison between different regions.

(3) The phylogenetic tress of samples and the mutations clsuters across samples.

(4) The clonal orders and the clonal evolution history.


### 1. Estimating the CCF values based on CHAT
**Scripts**: CCF_adjust.R

### 2. JSI index
**Scripts**: 
`meskit_routines.R` Fig 5B-G

Fig 5B: applied JSI index to estimate the divergence levels between ACA and NEC

Fig 5C-D: The paired JSI in patients with MSR (multi-region sequencing)

Fig 5E: The relationships between JSI and divergence time (T index)

Fig 5F: Distributions of JSI in 27 SSR patients.

Fig 5G: the neutral evolutionary index.


### 3. Whole Genome Doubling (WGD)
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
**Scripts**: 
`GC_revolver.R `

Fig 6B-D


`GC_revolver_surivival.R`

Supplementary Fig S13

