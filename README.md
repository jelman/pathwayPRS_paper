### pathwayPRS_paper
Code used in the manuscript investigating associations between pathway-specific PRS, Alzheimer's disease diagnosis and ATN biomarkers

### Authors
* **Nicholas J. Schork** - *The Translational Genomics Research Institute*
* **Jeremy A. Elman** - *University of California San Diego*

### Citation
Schork, NJ, Elman,. Pathway-specific polygenic risk scores correlate with clinical status and Alzheimer’s-related biomarkers, 01 March 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2583037/v1]

### Prerequisites:
1. [MAGMA](https://ctg.cncr.nl/software/magma)
2. [PRSice](https://choishingwan.github.io/PRSice/)
3. [Cytoscape](https://cytoscape.org/)
4. [EnrichmentMap](https://apps.cytoscape.org/apps/enrichmentmap)
5. [AutoAnnotate](https://apps.cytoscape.org/apps/autoannotate)
6. [ADNIMERGE](https://adni.bitbucket.io/reference/ADNIMERGE-package.html)
7. [R](https://www.r-project.org/)
9. [PLINK](https://www.cog-genomics.org/plink/)


### Data requirements
1. Summary statistics from the Kunkle et al. 2019 GWAS of AD (PMID: 30820047) are available from [NIAGADS](https://www.niagads.org/igap-rv-summary-stats-kunkle-p-value-data)
2. GO (no iea), KEGG, and Reactome gene sets from [Bader Lab](http://download.baderlab.org/EM_Genesets/February_01_2022/Human/) (February 2, 2022)
3. Download the following files from the [MAGMA webite](https://ctg.cncr.nl/software/magma):
   - Gene locations, build 37
   - SNP synonyms, dbSNP 151
   - 1000 Genomes EUR reference data 
4. Gene boundaries from [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz)
5. ADNIMERGE R Package from the [ADNI website](http://adni.loni.usc.edu/).  


### Code
1. `create_genesets.sh` - Clean, format, and concatenate gene sets from Bader Lab
2. `run_magma.sh` - Run MAGMA gene analysis on AD GWAS summary statistics
3. `magma2enrichmentMap.R` - Convert MAGMA output to EnrichmentMap input
4. `createPathwayClusters.R` - Create pathway clusters using Cytoscape apps EnrichmentMap and AutoAnnotate
5. `cytoscape2gmt.R` - Convert Cytoscape cluster output to GMT file for use in PRSice
6. `run_PRSice.sh` - Calculate pathway-specific PRS for ADNI participants using PRSice and AD GWAS summary statistics
7. `prsAssociation.R` - Run linear regression models to test associations between pathway-specific PRS and ATN biomarkers

