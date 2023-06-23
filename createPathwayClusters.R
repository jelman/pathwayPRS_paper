library(RCy3)
# installApp("EnrichmentMap Pipeline Collection") # installs apps if not installed

# Set data directory
data.dir = "/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data"
# Set gmt filename
gmt.file = file.path(data.dir, "Genesets", "Human_GOBP_REACTOME_KEGG_February_01_2022_symbol_Trimmed_10-1000.gmt")

# Set enrichment filename
enrichments.file = file.path(data.dir, "MAGMA", "kunkle_magma_enrichment_results_Trimmed_1000.txt")

# Set output filename
nodetable.file = file.path(data.dir, "cytoscape", "Kunkle_Magma_BaderLab_Trimmed_10-1000_q25_o5_SummaryNetwork-ClustAndUnclust.csv")

# Run EnrichmentMap build command
em_command = paste('enrichmentmap build analysisType="generic"',
                    "gmtFile=", gmt.file, 
                    "enrichmentsDataset1=", enrichments.file,
                    "pvalue=", 1,
                    "qvalue=", .25,
                    "similaritycutoff=",0.5,
                    "coefficients=","OVERLAP")
print(em_command)
commandsGET(em_command)

# Run the AutoAnnotate command to get pathway clusters
aa_command = paste("autoannotate annotate-clusterBoosted",
                    "clusterAlgorithm=MCL",
                    "labelColumn=EnrichmentMap::GS_DESCR",
                    "createSingletonClusters=true")
print(aa_command)
commandsGET(aa_command)

# Create summary network
aa_command = paste("autoannotate summary", 
                    "includeUnclustered=true")
print(aa_command)
commandsGET(aa_command)

# Export node table
aa_command = paste("table export",
                    "options=CSV",
                    "outputFile=", nodetable.file,
                    "table=", "AutoAnnotate - Summary Network default node")
print(aa_command)
commandsGET(aa_command)