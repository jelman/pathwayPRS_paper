library(dplyr)
library(tidylog)
library(readr)

# Data directory
data.dir = "/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data"
# Set magma file name
magma.file = file.path(data.dir, "MAGMA", "kunkle_magma_Trimmed_10-1000.gsa.out")

# Load magma output
print(paste0("Reading MAGMA output from ", magma.file))
df = read.table(magma.file, sep="", comment.char = "#", header=TRUE)

# Get FDR corrected p-values
print("Calculating FDR corrected p-values")
df$P_fdr = p.adjust(df$P, method="fdr")

# Load gene sets to get descriptions
geneset.file = file.path(data.dir, "Genesets", "Human_GOBP_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt")
print(paste0("Reading gene set file from ", geneset.file))
gs = read_delim(geneset.file, delim='\t', col_names=c("ID", "Description"), col_types = cols_only(ID = 'c', Description = 'c' ))

# Join gene set descriptions with MAGMA results
print("Joining gene set descriptions with MAGMA results")
newdf = df %>% left_join(gs, by=c("VARIABLE"="ID"))
newdf = newdf %>% dplyr::select(ID=VARIABLE, Description, P, P_fdr)

# Write out results
outfile.name = file.path(data.dir, "MAGMA", "kunkle_magma_enrichment_results_Trimmed_1000.txt")
print(paste0("Writing MAGMA results to ", outfile.name))
write_delim(newdf, outfile.name, delim = "\t")
