#!/bin/sh

# Set data directory
DATADIR=/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data

# Prerequisite: Download GO (no iea), KEGG, and Reactome gene sets from Bader Lab http://download.baderlab.org/EM_Genesets/February_01_2022/Human/

# Find and replace "\t\t" from KEGG file
sed -i 's/\t\t//g' ${DATADIR}/GeneSets/Human_KEGG_February_01_2022_symbol.gmt
# Sort GO file so that IDs in in multiple categories show up with BP first
sort -t'%' -s -k 1,1q ${DATADIR}/Genesets/Human_GOALL_no_GO_iea_February_01_2022_symbol.gmt > ${DATADIR}/Genesets/Human_GOALL_no_GO_iea_February_01_2022_symbol_sorted.gmt
# merge gene set files into one
cat ${DATADIR}/GeneSets/Human_GOALL_no_GO_iea_February_01_2022_symbol_sorted.gmt ${DATADIR}/GeneSets/Human_KEGG_February_01_2022_symbol.gmt ${DATADIR}/GeneSets//Human_Reactome_February_01_2022_symbol.gmt | sed '/^$/d' > ${DATADIR}/GeneSets//Human_GO_REACTOME_KEGG_February_01_2022_symbol.gmt
# Remove first part of ID field
sed -i 's/.*%//g' ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol.gmt 

# Filter for gene sets with 10-1000 genes and remove duplicated IDs
awk -F "\t" '(NF>12 && NF<1004) {print}' ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol.gmt > ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt
sort -buk1,1 ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt > ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000_unique.gmt
mv ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000_unique.gmt ${DATADIR}/GeneSets/Human_GO_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt