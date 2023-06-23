#!/bin/sh

# Prerequisites: Download GTF with with gene boundaries for use with PRSet from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz

# Set data directory
DATADIR=/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data

# Run PRSet to calculte pathway-specific PRS
/usr/local/PRSice/PRSice_linux \
    --a1 Effect_allele \
    --a2 Non_Effect_allele \
    --allow-inter  \
    --bar-levels 1 \
    --base ${DATADIR}/Kunkle_2019_IGAP2/Kunkle_etal_Stage1_results.txt \
    --beta  \
    --binary-target T \
    --clump-kb 500kb \
    --clump-p 1.000000 \
    --clump-r2 0.200000 \
    --fastscore  \
    --feature exon,gene,protein_coding,CDS \
    --gtf ${DATADIR}/gencode.v26lift37.annotation.gtf.gz \
    --ignore-fid  \
    --info 0.5 \
    --ld ${DATADIR}/MAGMA/g1000_eur \
    --ld-geno 0.1 \
    --ld-maf 0.01 \
    --maf 0.01 \
    --model add \
    --msigdb ${DAADIR}/cytoscape/Kunkle_Magma_BaderLab_Trimmed_10-1000_q25_o5_ClusteredPathwayGeneSets_Symbols.gmt \
    --num-auto 22 \
    --out ${DATADIR}/PRS/ADNI/ADNI_Kunkle_ClusteredPathways_q25_o5_p1 \
    --print-snp  \
    --pvalue Pvalue \
    --score avg \
    --seed 1234 \
    --set-perm 10000 \
    --snp MarkerName \
    --stat Beta \
    --target ${DATADIR}/ADNI_genetic_data/ADNI_allchr_EUR \
    --use-ref-maf  \
    --wind-3 10kb \
    --wind-5 35kb


## Create pathway PRS excluding APOE region (Chr19:45,116,911–46,318,605) ##

# Get list of SNPs in Kunkle GWAS excluding APOE region
awk '!(($2 == 19) && ($3 >= 45116911 && $3 <= 46318605))' ${DATADIR}/Kunkle_2019_IGAP2/kunkle_snplocs.txt > ${DATADIR}/Kunkle_2019_IGAP2/kunkle_snplocs_NoAPOE.txt

# Run PRSet to calculte pathway-specific PRS excluding APOE region
/usr/local/PRSice/PRSice_linux \
    --a1 Effect_allele \
    --a2 Non_Effect_allele \
    --allow-inter  \
    --bar-levels 1 \
    --base ${DATADIR}/Kunkle_2019_IGAP2/Kunkle_etal_Stage1_results.txt \
    --beta  \
    --binary-target T \
    --clump-kb 500kb \
    --clump-p 1.000000 \
    --clump-r2 0.200000 \
    --fastscore  \
    --feature exon,gene,protein_coding,CDS \
    --gtf ${DATADIR}/gencode.v26lift37.annotation.gtf.gz \
    --ignore-fid  \
    --info 0.5 \
    --ld ${DATADIR}/MAGMA/g1000_eur \
    --ld-geno 0.1 \
    --ld-maf 0.01 \
    --maf 0.01 \
    --model add \
    --msigdb ${DAADIR}/cytoscape/Kunkle_Magma_BaderLab_Trimmed_10-1000_q25_o5_ClusteredPathwayGeneSets_Symbols.gmt \
    --num-auto 22 \
    --out ${DATADIR}/PRS/ADNI/ADNI_Kunkle_ClusteredPathways_q25_o5_p1 \
    --print-snp  \
    --pvalue Pvalue \
    --score avg \
    --seed 1234 \
    --set-perm 10000 \
    --snp MarkerName \
    --stat Beta \
    --target ${DATADIR}/ADNI_genetic_data/ADNI_allchr_EUR \
    --use-ref-maf  \
    --wind-3 10kb \
    --wind-5 35kb \
    --x-range 19:45116911-46318605