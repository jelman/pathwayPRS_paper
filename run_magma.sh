#!/bin/sh

# Prerequisites: 
# Download file from MAGMA webite https://ctg.cncr.nl/software/magma
#   1. Gene locations, build 37
#   2. SNP synonyms, dbSNP 151
#   2. 1000 Genomes EUR reference data 
# Download Kunkle summary statistics from https://www.niagads.org/igap-rv-summary-stats-kunkle-p-value-data

# Set data directory
DATADIR=${DATADIR}

# Create version of gene location file with gene symbols 
awk '{print $6, $2, $3, $4, $5, $1}' OFS='\t' ${DATADIR}/MAGMA/NCBI37.3/NCBI37.3.gene.loc > ${DATADIR}/MAGMA/NCBI37.3/NCBI37.3.gene.symbol.loc

# Get list of SNPs in Kunkle GWAS
zcat ${DATADIR}/Kunkle_2019_IGAP2/Kunkle_etal_Stage1_results.txt.gz | awk '{ print $3"\t"$$1"\t"$2 }' > ${DATADIR}/Kunkle_2019_IGAP2/kunkle_snplocs.txt

# Annote snps
magma --annotate window=35,10 \
--snp-loc ${DATADIR}/Kunkle_2019_IGAP2/kunkle_snplocs.txt \
--gene-loc ${DATADIR}/MAGMA/NCBI37.3/NCBI37.3.gene.symbol.loc \
--out ${DATADIR}/MAGMA/Kunkle/kunkle_annot

# Gene analysis (SNPwise-mean)
magma --bfile ${DATADIR}/MAGMA/g1000_eur \
--pval ${DATADIR}/Kunkle_2019_IGAP2/Kunkle_etal_Stage1_results.txt use=3,8 N=63926 \
--gene-annot ${DATADIR}/MAGMA/Kunkle/kunkle_annot.genes.annot \
--out ${DATADIR}/MAGMA/Kunkle/kunkle_magma

# Gene set analysis
magma --gene-results ${DATADIR}/MAGMA/Kunkle/kunkle_magma.genes.raw \
--set-annot ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt \
--out ${DATADIR}/MAGMA/kunkle_magma_Trimmed_10-1000
