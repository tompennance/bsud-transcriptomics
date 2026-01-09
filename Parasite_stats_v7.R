# parasite stats but looking at which genes specifically aligned for samples #

# HISAT2_Parasite alignment stats
library(ggplot2)
#install.packages("dunn.test")
library(dunn.test)
#install.packages("ARTool")
library(ARTool)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)
library(readr)

setwd("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/Parasite_transcripts_analysis")

#### IMPORT AND FORMAT DATA ####

# get table that provides overall summary of reads aligned to S. mansoni genome
# read table generated from script: "/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis/Summary_alignment_stats_COMPLETE_script.R"
Sm_alignment_stats <- read.table(file ="/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/data/alignment_stats/Summary_HISAT2_Alignment_Stats_with_raw_read_data.tab", sep = "\t", header = TRUE)
# convert to factors
Sm_alignment_stats$Time = factor(Sm_alignment_stats$Time)
Sm_alignment_stats$Treatment = factor(Sm_alignment_stats$Treatment)

Sm_featureCounts <- read.csv(file = "/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/data/alignment_stats/Sm_featureCounts_formatted.csv", header = TRUE, sep = ",", row.names = "Geneid")
# transpose datadrame to make rows = samples and columns = genes

Sm_featureCounts_transposed <- t(Sm_featureCounts)
Sm_featureCounts_transposed <- as.data.frame(Sm_featureCounts_transposed)

# combine Sm_alignment_stats and Sm_featureCounts
# Convert row names of Sm_featureCounts into a new column (e.g., "SAMPLE")
Sm_featureCounts_transposed$SAMPLE <- rownames(Sm_featureCounts_transposed)
# Merge Sm_alignment_stats and Sm_featureCounts_transposed by matching the 'SAMPLE' column in both
merged_df <- merge(Sm_featureCounts_transposed, Sm_alignment_stats, by = "SAMPLE")

# tidy enviroment
rm(Sm_featureCounts_transposed)

# Import S. mansoni genome annotation ID's

Sm_gene_annotations <- read_delim("/Users/tpennance/Documents/Papers/Schistosoma_mansoni/Genome/gene_annotations.tab", 
                                  delim = "\t", 
                                  col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
                                  comment = "#")

# rework attributes column to define a column per variable
Sm_gene_annotations <- Sm_gene_annotations %>%
  mutate(
    GeneID = str_extract(attributes, "(?<=ID=gene:)Smp_\\d+(?=;)"),
    Name = str_extract(attributes, "(?<=Name=)Smp_\\d+(?=;)"),
    biotype = str_extract(attributes, "(?<=biotype=)[^;]+"),
    description = str_extract(attributes, "(?<=description=)[^;]+"),
    description_source = str_extract(attributes, "(?<=description_source=)[^;]+"),
    description_source_acc = str_extract(attributes, "(?<=description_source_acc=)[^;]+"),
    locus = str_extract(attributes, "(?<=locus=)[^;]+"),
    previous_stable_id = str_extract(attributes, "(?<=previous_stable_id=)[^;]+")
  )

# add S. mansoni UniProt data
Sm_gene_UniProtKB <- read_xlsx("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/Parasite_transcripts_analysis/UniProt_Data/idmapping_2024_11_10.xlsx")
colnames(Sm_gene_UniProtKB)[1] <- "description_source_acc"
# Clean column names: replace spaces with underscores and remove any brackets
colnames(Sm_gene_UniProtKB) <- colnames(Sm_gene_UniProtKB) %>%
  gsub(" ", "_", .) %>%          # Replace spaces with underscores
  gsub("\\[|\\]|\\(|\\)", "", .) # Remove brackets (both [] and ())

Sm_gene_annotations_UniProt <- Sm_gene_annotations %>%
  left_join(Sm_gene_UniProtKB, by = "description_source_acc")

# tidy enviroment
rm(Sm_gene_annotations)
rm(Sm_gene_UniProtKB)
rm(Sm_alignment_stats)

# add Smp annotation to feature count file 
Sm_featureCounts$GeneID <- rownames(Sm_featureCounts)
Sm_featureCounts_annotated <- Sm_featureCounts %>%
  left_join(Sm_gene_annotations_UniProt, by = "GeneID")

# tidy
rm(Sm_featureCounts)

#### REMOVAL OF CONTROL TRANSCRIPTS and generating all sample gene summary ####

# Identify gene columns (those starting with 'Smp')
gene_cols <- grep("^Smp", names(merged_df), value = TRUE)

# Create functions to count presence (≥1 read) and sum reads for each gene
count_presence <- function(x) sum(x >= 1)
sum_reads <- function(x) sum(x)

# Group by Treatment and Time, then summarize gene presence and total reads
# Could modify this to resistance / susceptible as determined from other dataset
gene_summary <- merged_df %>%
  group_by(Treatment, Time) %>%
  summarise(across(all_of(gene_cols), 
                   list(SamplesWithReads = count_presence,
                        TotalReads = sum_reads)),
            .groups = "drop") %>%
  pivot_longer(cols = -c(Treatment, Time),
               names_to = c("GeneID", ".value"),
               names_pattern = "(.+)_(SamplesWithReads|TotalReads)")
write_xlsx(gene_summary, "")

# determine the transcripts aligning in control snails
Control_transcripts <- gene_summary %>%
  filter(Treatment == "Control" & TotalReads >= 1) %>% 
  group_by(GeneID) %>%
  summarise(
    TotalReads = sum(TotalReads, na.rm = TRUE),
    SamplesWithReads = sum(SamplesWithReads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalReads))

# merge control transcripts with S. mansoni annotations and make a table for all
Control_transcripts <- merge(Control_transcripts, Sm_gene_annotations_UniProt, by = "GeneID", all.x = TRUE)
write_xlsx(Control_transcripts, "/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/Parasite_transcripts_analysis/Transcripts_aligning_to_control_snails_v2.xlsx")

control_gene_ids <- Control_transcripts$GeneID

# Remove columns from merged_df that match the gene IDs in control_gene_ids
merged_df_no_control_transcripts <- merged_df[, !colnames(merged_df) %in% control_gene_ids]
Sm_featureCounts_annotated_no_control_transcripts <- Sm_featureCounts_annotated %>%
  filter(!GeneID %in% control_gene_ids)
write_xlsx(Sm_featureCounts_annotated_no_control_transcripts, path = "Sm_featureCounts_annotated_no_control_transcripts.xlsx")

# create new gene columns (those starting with 'Smp') WITHOUT control transcripts
gene_cols <- grep("^Smp", names(merged_df_no_control_transcripts), value = TRUE)

# make a new gene summary with no control transcripts
gene_summary <- merged_df %>%
  group_by(Treatment, Time) %>%
  summarise(across(all_of(gene_cols), 
                   list(SamplesWithReads = count_presence,
                        TotalReads = sum_reads)),
            .groups = "drop") %>%
  pivot_longer(cols = -c(Treatment, Time),
               names_to = c("GeneID", ".value"),
               names_pattern = "(.+)_(SamplesWithReads|TotalReads)")

# tidy environment
rm(Control_transcripts)
rm(control_gene_ids)
rm(gene_cols)

##### PLOT READ COUNTS WITH ALL TRANSCRIPTS ####

# Step 1: Sum read counts for each sample and extract Treatment and Time information
sample_totals <- Sm_featureCounts_annotated %>%
  select(starts_with("KC"), starts_with("KK"), starts_with("KN")) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "SAMPLE", values_to = "Total_Reads") %>%
  mutate(
    Treatment = case_when(
      str_starts(SAMPLE, "KK") ~ "UNMKENYA",
      str_starts(SAMPLE, "KN") ~ "NMRI",
      str_starts(SAMPLE, "KC") ~ "Control"
    ),
    Time = case_when(
      str_detect(SAMPLE, "KK8") ~ "8",
      str_detect(SAMPLE, "KC8") ~ "8",
      str_detect(SAMPLE, "KN8") ~ "8",
      str_detect(SAMPLE, "KK24") ~ "24",
      str_detect(SAMPLE, "KC24") ~ "24",
      str_detect(SAMPLE, "KN24") ~ "24",
      str_detect(SAMPLE, "KC72") ~ "72",
      str_detect(SAMPLE, "KK72") ~ "72",
      str_detect(SAMPLE, "KN72") ~ "72"
    )
  ) %>%
  arrange(Time, Treatment)  # Arrange by Treatment, Time, and SAMPLE

sample_totals <- sample_totals %>%
  mutate(
    Treatment_Time = factor(
      interaction(Treatment, Time),
      levels = c("Control.8", "UNMKENYA.8", "NMRI.8", 
                 "Control.24", "UNMKENYA.24", "NMRI.24", 
                 "Control.72", "UNMKENYA.72", "NMRI.72")
    ),
    SAMPLE = factor(SAMPLE, levels = SAMPLE[order(Treatment_Time)])
  )

# Plot with the reordered SAMPLE factor
pdf('Sm_Read_count_per_sample.pdf', width=20, height=10, useDingbats=FALSE)
p <- ggplot(sample_totals, aes(x = SAMPLE, y = Total_Reads, fill = Treatment)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample ID", y = "Total Read Counts", fill = "Treatment-Time") +
  theme_minimal() +
  ylim(0,5050) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("Control" = "black", "NMRI" = "red", "UNMKENYA" = "lightgreen"))  # Set custom colors
print(p)
dev.off()

# tidy enviroment 
rm(p)
rm(sample_totals)

##### PLOT READ COUNTS WITHOUT CONTROL TRANSCRIPTS ####
sample_totals_no_control <- Sm_featureCounts_annotated_no_control_transcripts %>%
  dplyr::select(starts_with("KC"), starts_with("KK"), starts_with("KN")) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "SAMPLE", values_to = "Total_Reads") %>%
  mutate(
    Treatment = case_when(
      str_starts(SAMPLE, "KK") ~ "UNMKENYA",
      str_starts(SAMPLE, "KN") ~ "NMRI",
      str_starts(SAMPLE, "KC") ~ "Control"
    ),
    Time = case_when(
      str_detect(SAMPLE, "KK8") ~ "8",
      str_detect(SAMPLE, "KC8") ~ "8",
      str_detect(SAMPLE, "KN8") ~ "8",
      str_detect(SAMPLE, "KK24") ~ "24",
      str_detect(SAMPLE, "KC24") ~ "24",
      str_detect(SAMPLE, "KN24") ~ "24",
      str_detect(SAMPLE, "KC72") ~ "72",
      str_detect(SAMPLE, "KK72") ~ "72",
      str_detect(SAMPLE, "KN72") ~ "72"
    )
  ) %>%
  arrange(Time, Treatment)  # Arrange by Treatment, Time, and SAMPLE

sample_totals_no_control <- sample_totals_no_control %>%
  mutate(
    Treatment_Time = factor(
      interaction(Treatment, Time),
      levels = c("Control.8", "UNMKENYA.8", "NMRI.8", 
                 "Control.24", "UNMKENYA.24", "NMRI.24", 
                 "Control.72", "UNMKENYA.72", "NMRI.72")
    ),
    SAMPLE = factor(SAMPLE, levels = SAMPLE[order(Treatment_Time)])
  )

# Plot with the reordered SAMPLE factor
pdf('Sm_Read_count_per_sample_no_control_transcripts.pdf', width=20, height=10, useDingbats=FALSE)
p <- ggplot(sample_totals_no_control, aes(x = SAMPLE, y = Total_Reads, fill = Treatment)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample ID", y = "Total Read Counts", fill = "Treatment-Time") +
  theme_minimal() +
  ylim(0,5050) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("Control" = "black", "NMRI" = "red", "UNMKENYA" = "lightgreen"))  # Set custom colors
print(p)
dev.off()


# tidy enviroment 
rm(p)
rm(sample_totals_no_control)


#### Create summary for each gene and experimental group ####

# NMRI 8 hours
Sm_transcripts_NMRI_8 <- gene_summary %>%
  filter(Treatment == "NMRI", Time == "8") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_NMRI_8) <- c("GeneID", "Samples_NMRI_8", "Reads_NMRI_8")
# UNMKENYA 8 hours
Sm_transcripts_UNMKENYA_8 <- gene_summary %>%
  filter(Treatment == "UNMKENYA", Time == "8") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_UNMKENYA_8) <- c("GeneID", "Samples_UNMKENYA_8", "Reads_UNMKENYA_8")
# NMRI 24 hours
Sm_transcripts_NMRI_24 <- gene_summary %>%
  filter(Treatment == "NMRI", Time == "24") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_NMRI_24) <- c("GeneID", "Samples_NMRI_24", "Reads_NMRI_24")
# UNMKENYA 24 hours
Sm_transcripts_UNMKENYA_24 <- gene_summary %>%
  filter(Treatment == "UNMKENYA", Time == "24") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_UNMKENYA_24) <- c("GeneID", "Samples_UNMKENYA_24", "Reads_UNMKENYA_24")
# NMRI 72 hours
Sm_transcripts_NMRI_72 <- gene_summary %>%
  filter(Treatment == "NMRI", Time == "72") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_NMRI_72) <- c("GeneID", "Samples_NMRI_72", "Reads_NMRI_72")
# UNMKENYA 72 hours
Sm_transcripts_UNMKENYA_72 <- gene_summary %>%
  filter(Treatment == "UNMKENYA", Time == "72") %>%
  group_by(GeneID) %>%
  summarise(TotalSamplesWithReads = sum(SamplesWithReads),
            GrandTotalReads = sum(TotalReads)) %>%
  arrange(desc(TotalSamplesWithReads), desc(GrandTotalReads))
colnames(Sm_transcripts_UNMKENYA_72) <- c("GeneID", "Samples_UNMKENYA_72", "Reads_UNMKENYA_72")

# merge datasets
merged_data <- Sm_transcripts_NMRI_8 %>%
  inner_join(Sm_transcripts_NMRI_24, by = "GeneID")
merged_data <- merged_data %>%
  inner_join(Sm_transcripts_NMRI_72, by = "GeneID")
merged_data <- merged_data %>%
  inner_join(Sm_transcripts_UNMKENYA_8, by = "GeneID")
merged_data <- merged_data %>%
  inner_join(Sm_transcripts_UNMKENYA_24, by = "GeneID")
merged_data <- merged_data %>%
  inner_join(Sm_transcripts_UNMKENYA_72, by = "GeneID")

# remove genes with only 0's?
merged_data <- merged_data %>%
  filter(if_any(-1, ~. != 0))

# append gene annotation information
Sm_transcript_ALL_no_control_transcripts <- merged_data %>%
  left_join(Sm_gene_annotations_UniProt, by = "GeneID")

write_xlsx(Sm_transcript_ALL_no_control_transcripts, path = "/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/Parasite_transcripts_analysis/Sm_transcript_ALL_no_control_transcripts.xlsx")

# tidy environment
rm(merged_data)
rm(Sm_transcripts_NMRI_24)
rm(Sm_transcripts_NMRI_72)
rm(Sm_transcripts_NMRI_8)
rm(Sm_transcripts_UNMKENYA_8)
rm(Sm_transcripts_UNMKENYA_24)
rm(Sm_transcripts_UNMKENYA_72)

#### plot reads per gene across time ####

# Step 1: Gather read count data for each Time and Treatment combination
read_counts_long <- Sm_transcript_ALL_no_control_transcripts %>%
  select(GeneID, starts_with("Reads")) %>% 
  pivot_longer(
    cols = starts_with("Reads"),
    names_to = c("Treatment", "Time"),
    names_pattern = "Reads_(\\w+)_(\\d+)",
    values_to = "ReadCount"
  )

# Step 3: Plot ReadCount by GeneID, Time, and Treatment
# Aggregate to get mean ReadCount per GeneID, Treatment, and Time
read_counts <- read_counts_long %>%
  group_by(GeneID, Treatment, Time) %>%
  summarise(ReadCount, .groups = "drop")

# Plot each gene’s read count over time for each treatment
pdf('Sm_Read_count_per_transcript_line_graph.pdf', width=10, height=10, useDingbats=FALSE)
p <- ggplot(read_counts, aes(x = as.numeric(Time), y = ReadCount, color = GeneID, group = GeneID)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Treatment) +
  labs(
    title = "Schistosoma mansoni transcripts across time and treatment",
    x = "Time (hours)",
    y = "Read Count"
  ) +
  scale_x_continuous(breaks = c(8, 24, 72)) +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
dev.off()

# tidy enviroment 
rm(read_counts)
rm(read_counts_long)
rm(p)

#### Summarize GO terms per sample containing Sm transcripts function ####

# Pick which samples to plot - either do all KN and KK snails using below:
#sample_ids <- colnames(Sm_featureCounts_annotated_no_control_transcripts) %>%
#  .[starts_with(c("KK", "KN"), vars = .)] 

# or select those with a minimum transcript number
# Step 1: Calculate total reads for each column that starts with "KK" or "KN"
total_reads_per_sample <- Sm_featureCounts_annotated_no_control_transcripts %>%
  select(starts_with("KK"), starts_with("KN")) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))
write_xlsx(total_reads_per_sample, "total_Sm_reads_per_sample_no_control.xlsx")

# Step 2: Filter for columns with total reads >= 100
sample_ids_of_Sm_infected <- names(total_reads_per_sample)[total_reads_per_sample >= 100]

# Now `sample_ids` contains only columns starting with "KK" or "KN" that have total reads >= 100
print(sample_ids_of_Sm_infected)

### Summarise BIOLOGICAL PROCESS GO terms ###

# Define a function to generate and save GO term plots for a list of sample IDs
plot_GO_terms_list <- function(data, sample_ids_of_Sm_infected, top_n = 20, output_dir = "GO_term_plots") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop through each sample ID in the list
  for (sample_id in sample_ids_of_Sm_infected) {
    # Step 1: Replace blank GO terms with "Unknown" and separate multiple GO terms into individual rows
    Sm_featureCounts_long <- data %>%
      select(GeneID, Gene_Ontology_biological_process, !!sym(sample_id)) %>%
      rename(Sample_Reads = !!sym(sample_id)) %>%  # Rename the column for generality
      mutate(
        Gene_Ontology_biological_process = ifelse(
          Gene_Ontology_biological_process == "" | is.na(Gene_Ontology_biological_process),
          "Unknown",
          Gene_Ontology_biological_process
        )
      ) %>%  # Replace blanks or NAs with "Unknown"
      mutate(GO_Term = strsplit(Gene_Ontology_biological_process, ";")) %>%  # Split GO terms by ';'
      unnest(GO_Term) %>%  # Convert list to long format
      mutate(GO_Term = trimws(GO_Term))  # Remove any extra whitespace
    
    # Step 2: Summarize the total read counts for each GO term and concatenate GeneIDs
    GO_counts <- Sm_featureCounts_long %>%
      group_by(GO_Term) %>%
      summarise(
        Total_Reads = sum(Sample_Reads, na.rm = TRUE),
        GeneIDs = paste(unique(GeneID), collapse = ";")  # Concatenate unique GeneIDs with ";"
      ) %>%
      arrange(desc(Total_Reads))
    
    # Step 3: Filter top GO_Terms
    GO_counts_top <- GO_counts %>%
      filter(GO_Term != "Unknown") %>%  # Exclude "Unknown" terms
      slice_max(Total_Reads, n = top_n)  # Select the top GO terms by Total_Reads
    
    # Step 4: Generate the plot
    p <- ggplot(GO_counts_top, aes(x = reorder(GO_Term, Total_Reads), y = Total_Reads)) +
      geom_bar(stat = "identity") +
      coord_flip() +  # Flip coordinates for better readability
      labs(
        x = "GO Term (Biological Process)",
        y = "Total Read Counts",
        title = paste("Top", top_n, "Read Counts by GO Term for Sample", sample_id)
      ) +
      theme_minimal()
    
    # Step 5: Save the plot as a PDF
    pdf_filename <- paste0(output_dir, "/", sample_id, "_Sm_top", top_n, "_GO_BP_terms.pdf")
    pdf(pdf_filename, width = 10, height = 10, useDingbats = FALSE)
    print(p)
    dev.off()
    
    message("Plot saved for sample ", sample_id, " at ", pdf_filename)
  }
}

# generate and save plots
plot_GO_terms_list(Sm_featureCounts_annotated_no_control_transcripts, sample_ids_of_Sm_infected)



### generate a function to plot MOLECULAR FUNCTIONS Gene_Ontology_molecular_function

plot_GO_MF_terms_list <- function(data, sample_ids_of_Sm_infected, top_n = 20, output_dir = "GO_term_plots") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop through each sample ID in the list
  for (sample_id in sample_ids_of_Sm_infected) {
    # Step 1: Replace blank GO terms with "Unknown" and separate multiple GO terms into individual rows
    Sm_featureCounts_long <- data %>%
      select(GeneID, Gene_Ontology_molecular_function, !!sym(sample_id)) %>%
      rename(Sample_Reads = !!sym(sample_id)) %>%  # Rename the column for generality
      mutate(
        Gene_Ontology_molecular_function = ifelse(
          Gene_Ontology_molecular_function == "" | is.na(Gene_Ontology_molecular_function),
          "Unknown",
          Gene_Ontology_molecular_function
        )
      ) %>%  # Replace blanks or NAs with "Unknown"
      mutate(GO_Term = strsplit(Gene_Ontology_molecular_function, ";")) %>%  # Split GO terms by ';'
      unnest(GO_Term) %>%  # Convert list to long format
      mutate(GO_Term = trimws(GO_Term))  # Remove any extra whitespace
    
    # Step 2: Summarize the total read counts for each GO term and concatenate GeneIDs
    GO_counts <- Sm_featureCounts_long %>%
      group_by(GO_Term) %>%
      summarise(
        Total_Reads = sum(Sample_Reads, na.rm = TRUE),
        GeneIDs = paste(unique(GeneID), collapse = ";")  # Concatenate unique GeneIDs with ";"
      ) %>%
      arrange(desc(Total_Reads))
    
    # Step 3: Filter top GO_Terms
    GO_counts_top <- GO_counts %>%
      filter(GO_Term != "Unknown") %>%  # Exclude "Unknown" terms
      slice_max(Total_Reads, n = top_n)  # Select the top GO terms by Total_Reads
    
    # Step 4: Generate the plot
    p <- ggplot(GO_counts_top, aes(x = reorder(GO_Term, Total_Reads), y = Total_Reads)) +
      geom_bar(stat = "identity") +
      coord_flip() +  # Flip coordinates for better readability
      labs(
        x = "GO Term (Molecular Function)",
        y = "Total Read Counts",
        title = paste("Top", top_n, "Read Counts by GO Term for Sample", sample_id)
      ) +
      theme_minimal()
    
    # Step 5: Save the plot as a PDF
    pdf_filename <- paste0(output_dir, "/", sample_id, "_Sm_top", top_n, "_GO_MF_terms.pdf")
    pdf(pdf_filename, width = 10, height = 10, useDingbats = FALSE)
    print(p)
    dev.off()
    
    message("Plot saved for sample ", sample_id, " at ", pdf_filename)
  }
}

# generate and save plots
plot_GO_MF_terms_list(Sm_featureCounts_annotated_no_control_transcripts, sample_ids_of_Sm_infected)

