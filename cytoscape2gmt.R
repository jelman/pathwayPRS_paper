library(dplyr)
library(tidyr)

# Set data directory
data.dir = "/home/jelman/netshare/BETELSHARE/PROJ/AD_genetic_subtypes/data"

# Read in gene set file
geneset.file = file.path(data.dir, "Genesets", "Human_GOBP_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt")
df = read.csv(geneset.file, stringsAsFactors = FALSE)

# Format select columns of interest and rename
df = df %>% dplyr::select(ClusterID=X__mclCluster, GenesSymbols=EnrichmentMap..Genes, ClusterName=name, nSets=cluster.node.count, nGenes=EnrichmentMap..gs_size) %>%
  mutate(ClusterName = gsub(" ", "_", ClusterName))
# Remove duplicate gene symbols
allSymbols = unique(unlist(strsplit(df$GenesSymbols, "|", fixed = T)))

#Create gmt format file
gmt = paste(df$ClusterName, df$GenesSymbols, sep="|")
gmt = gsub("|", "\t", gmt, fixed=TRUE)

# Write out gmt file
gmt.file = file.path(data.dir, "cytoscape", "Kunkle_Magma_BaderLab_Trimmed_10-1000_q25_o5_ClusteredPathwayGeneSets_Symbols.gmt")
write.table(gmt, file=gmt.file, quote=FALSE, row.names = FALSE, col.names = FALSE)