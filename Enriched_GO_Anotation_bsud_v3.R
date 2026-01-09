setwd("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ALL")

#BiocManager::install("edgeR")
BiocManager::install("topGO")

library(edgeR)
library(topGO)

#### Generate GO files for up and down regulated genes ####

# Set the input and output directory paths
input_dir <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/UP_n_DOWN_regulated_pw"

if (!dir.exists("GO_analysis")) {
  dir.create("GO_analysis")
  message("Directory 'GO_analysis' created.")
} else {
  message("Directory 'GO_analysis' already exists.")
}
output_dir <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/GO_analysis"

# List all .tab files in the input directory
file_list <- list.files(path = input_dir, pattern = "\\.tab$", full.names = TRUE)

# Loop through each file
for (file_path in file_list) {
  # Read the data
  data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Extract geneID for files ending in "_p0.05_GO.tab"
  gene_ids_0.05 <- data$geneID
  # Define the output file path for "_p0.05_GO.tab"
  output_file_name_0.05 <- paste0(sub("\\.tab$", "_p0.05_GO.tab", basename(file_path)))
  output_file_path_0.05 <- file.path(output_dir, output_file_name_0.05)
  
  # Write to a new file without the column header
  write.table(gene_ids_0.05, file = output_file_path_0.05, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Extract geneID for files ending in "_p0.001_GO.tab" where p0.001 column is 1
  if ("p0.001" %in% colnames(data)) {
    gene_ids_0.001 <- data$geneID[data$p0.001 == 1]
    
    # Define the output file path for "_p0.001_GO.tab"
    output_file_name_0.001 <- paste0(sub("\\.tab$", "_p0.001_GO.tab", basename(file_path)))
    output_file_path_0.001 <- file.path(output_dir, output_file_name_0.001)
    
    # Write to a new file without the column header, only if there are values
    if (length(gene_ids_0.001) > 0) {
      write.table(gene_ids_0.001, file = output_file_path_0.001, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    }
  }
}

# if any _GO.tab only has one gene in it (or none) it will cause an error that stops the loop.
### Remove these files from directory and run only with those with a minimum of two genes and it should work (if they have GO terms)

if (!dir.exists("GO_analysis_too_few_genes")) {
  dir.create("GO_analysis_too_few_genes")
  message("Directory 'GO_analysis_too_few_genes' created.")
} else {
  message("Directory 'GO_analysis_too_few_genes' already exists.")
}

# too few dir
too_few_dir <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/GO_analysis_too_few_genes"

# List all *_GO.tab files in the output dir
go_files <- list.files(output_dir, pattern = "_GO\\.tab$", full.names = TRUE)

# Loop through files and move those with <= 1 gene
for (f in go_files) {
  # Count number of lines (genes)
  n_lines <- length(readLines(f))
  
  if (n_lines <= 1) {
    file.rename(
      from = f,
      to   = file.path(too_few_dir, basename(f))
    )
  }
}

# clean up
rm(list = ls())

#### RUN GO analysis ####


# Define directory paths
go_analysis_dir <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/GO_analysis"
full_genome_path <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/Full_Genome.tab"
full_genome_id_path <- "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/Full_Genome.id"

# Get list of files in GO_analysis directory, remove "_GO.tab" suffix
Files <- gsub("_GO.tab", "", list.files(path = go_analysis_dir, pattern = "_GO.tab$"))


# Read annotation and gene IDs
annotationGO <- readMappings(full_genome_path)
totalGenes <- readLines(full_genome_id_path)


# determine what p value threshold for significance used for Fisher test
pvalue_good <- 0.05

# Define GO ontologies to test
go_classes <- c("BP", "CC", "MF") 

# Loop over each file in Files
for (spe in Files) {
  # Construct the full file path for the "_GO.tab" file
  Interes_genes <- file.path(go_analysis_dir, paste(spe, "_GO.tab", sep = ""))
  
  # Check if the file exists before reading
  if (!file.exists(Interes_genes)) {
    message("File does not exist: ", Interes_genes)
    next  # Skip to the next iteration if file is missing
  }
  
  # Read genes of interest from the file
  genes.interes <- readLines(Interes_genes)
  
  # Generate gene list factor
  geneList <- factor(as.integer(totalGenes %in% genes.interes))
  names(geneList) <- totalGenes
  
  # Loop over each GO class
  for (go_c in go_classes) {
    sampleGOdata <- new("topGOdata", ontology = go_c, allGenes = geneList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = annotationGO)
    resultFisher <- runTest(sampleGOdata, statistic = "fisher", algorithm = "Weight01")
    
    # Check for significant GO terms (CHECK WHERE FISHER )
    sig_go <- sum(score(resultFisher) <= pvalue_good)
    if (sig_go > 0) {
      GOtermsdeInteres.data.frame <- GenTable(sampleGOdata, P_Value = resultFisher, orderBy = "Weight01Fisher", ranksOf = "Weight01Fisher", topNodes = sig_go)
    } else {
      GOtermsdeInteres.data.frame <- data.frame()
    }
    
    # Define output filename and write table
    print_tab <- file.path(go_analysis_dir, paste(spe, "_", go_c, ".tab", sep = ""))
    write.table(GOtermsdeInteres.data.frame, print_tab, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# take tab produced from this and assess go terms - these are interpreted as:
# the enriched GO terms among the query list (up and down regulated) vs the background list
# What "background" means depends on what you put on that file. Could have been the entire genome or the genes with no read mappings.
# if list of GO terms is too big, one could process the results with ReviGO http://revigo.irb.hr/
# Revigo can take long lists of Gene Ontology terms and summarize them by removing redundant GO terms



