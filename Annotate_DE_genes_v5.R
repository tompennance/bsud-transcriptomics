## Annotate DE genes version 3 ## 4-26-24 ##
## This script developed to create a more comprehensive dataset for the significantly DE genes, including:
## pw comparison data, logFC values, FDR etc for the significant genes bu across every time point
## ignores the extra comparisons not used (i.e. C8vsC24 etc)

setwd("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/")

# packages required

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")
#library(biomaRt)
library(tidyverse)

library(stringr)
library(dplyr)

## Load data generated in edgeR_bsud.R script for significantly expressed genes:
ALL_DE_pwcomp_p0.05 <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/ALL_DE_pwcomp_p0.05.tab", header = TRUE, row.names = 1)
ALL_DE_pwcomp_p0.001 <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/ALL_DE_pwcomp_p0.001.tab", header = TRUE, row.names = 1)

## Load edgeR data for the pw comparisons of all relevant comparisons

C8vsK8_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsK8_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C8vsK8_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsK8")
colnames(C8vsK8_all_genes) <- new_column_names
colnames(C8vsK8_all_genes)[1] <- "geneID"

C8vsN8_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsN8_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C8vsN8_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsN8")
colnames(C8vsN8_all_genes) <- new_column_names
colnames(C8vsN8_all_genes)[1] <- "geneID"

C8vsK8Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsK8Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C8vsK8Sm_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsK8Sm")
colnames(C8vsK8Sm_all_genes) <- new_column_names
colnames(C8vsK8Sm_all_genes)[1] <- "geneID"

C8vsN8Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsN8Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C8vsN8Sm_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsN8Sm")
colnames(C8vsN8Sm_all_genes) <- new_column_names
colnames(C8vsN8Sm_all_genes)[1] <- "geneID"

C24vsK24_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsK24_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C24vsK24_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24")
colnames(C24vsK24_all_genes) <- new_column_names
colnames(C24vsK24_all_genes)[1] <- "geneID"

C24vsN24_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsN24_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C24vsN24_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsN24")
colnames(C24vsN24_all_genes) <- new_column_names
colnames(C24vsN24_all_genes)[1] <- "geneID"

C24vsK24Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsK24Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C24vsK24Sm_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24Sm")
colnames(C24vsK24Sm_all_genes) <- new_column_names
colnames(C24vsK24Sm_all_genes)[1] <- "geneID"

C72vsK72_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsK72_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C72vsK72_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72")
colnames(C72vsK72_all_genes) <- new_column_names
colnames(C72vsK72_all_genes)[1] <- "geneID"

C72vsN72_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsN72_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C72vsN72_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsN72")
colnames(C72vsN72_all_genes) <- new_column_names
colnames(C72vsN72_all_genes)[1] <- "geneID"

C72vsK72Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsK72Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
new_column_names <- colnames(C72vsK72Sm_all_genes)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72Sm")
colnames(C72vsK72Sm_all_genes) <- new_column_names
colnames(C72vsK72Sm_all_genes)[1] <- "geneID"

## Create dataframe of just geneID for the p0.05 genes
DE_gene_list_p0.05 <- ALL_DE_pwcomp_p0.05[, 1, drop = FALSE]
## Create column just for letting know if gene in p0.001 list
DE_gene_list_p0.05$p0.001 <- ifelse(ALL_DE_pwcomp_p0.05$geneID %in% ALL_DE_pwcomp_p0.001$geneID, 1, 0)

## Merge DE gene data
temp_merge_result_df <- merge(DE_gene_list_p0.05, C8vsK8_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C8vsK8Sm_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C24vsK24_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C24vsK24Sm_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C72vsK72_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C72vsK72Sm_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C8vsN8_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C8vsN8Sm_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C24vsN24_all_genes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, C72vsN72_all_genes, by = "geneID", all.x = TRUE)

#### Annotation data and description data to merge with ####

#Javiers annotation data
Annotation_ALL_DE_pwcomp_p0.05 <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/expresion_annotation_Javier/Annotation_ALL_DE_pwcomp_p0.05.tab", header = TRUE)
names(Annotation_ALL_DE_pwcomp_p0.05)[names(Annotation_ALL_DE_pwcomp_p0.05) == "GenID"] <- "geneID"

# New TMR table created by Javier 4-8-2024, listing TMRs and regions for largest isofrom
New_TMR_Table <- read.csv("/Users/tpennance/Documents/tpennanceWU/PostDoc/Biom_sud_annotation/TMHMM/New_TMR_Table.csv")
names(New_TMR_Table)[names(New_TMR_Table) == "Gen"] <- "geneID"
New_TMR_Table <- New_TMR_Table[, 1:8]

# Pennance et al. 2024 Supp Table 2 Key Gene orthologs (immune related genes)
KeyGenes <- read.csv("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_old/Supp_Table_2_Key_genes.csv", header = TRUE)
names(KeyGenes)[names(KeyGenes) == "Biomphalaria.sudanica.gene"] <- "geneID"
KeyGenes <- KeyGenes[, c(1, 7)]
KeyGenes[, 1] <- gsub("\\s*\\([^\\)]+\\)", "", KeyGenes[, 1])
names(KeyGenes)[names(KeyGenes) == "Candidate.gene...genomic.region..and.study.ies..where.first.identified..DOI."] <- "CandidateGene"

# BB02 orthologs
Bsud_BB02_orthologs <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/orthologs/Bsud_BB02_orthologs_v2.txt", header = TRUE)
# remove isoform information from BGLB protein ID's
Bsud_BB02_orthologs$GeneA_part <- str_extract(Bsud_BB02_orthologs$GeneA, "^[^-]+")
names(Bsud_BB02_orthologs)[names(Bsud_BB02_orthologs) == "GeneB"] <- "geneID"
names(Bsud_BB02_orthologs)[names(Bsud_BB02_orthologs) == "GeneA_part"] <- "BglabBB02_geneID"
names(Bsud_BB02_orthologs)[names(Bsud_BB02_orthologs) == "GeneA"] <- "BglabBB02_isoformID"

# Bp505 orthologs
Bsud_Bp505_orthologs <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/orthologs/Bsud_Bp505_orthologs_v2.txt", header = TRUE)
# remove isoform information from BGLB protein ID's
Bsud_Bp505_orthologs$GeneA_part <- str_extract(Bsud_Bp505_orthologs$GeneA, "^[^-]+")
names(Bsud_Bp505_orthologs)[names(Bsud_Bp505_orthologs) == "GeneB"] <- "geneID"
names(Bsud_Bp505_orthologs)[names(Bsud_Bp505_orthologs) == "GeneA_part"] <- "Bp505_geneID"
names(Bsud_Bp505_orthologs)[names(Bsud_Bp505_orthologs) == "GeneA"] <- "Bp505_isoformID"

# Pennance et al. 2024 Supp Table 8 Location signals, Transmembrane domains, InterPro signals
LocationSignals <- read.csv("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_old/Supp_Table_8_Location_Signals.csv", header =TRUE)
# remove isoform infomration and duplicated entries
LocationSignals$ProtID <- sub("\\.(\\d+)\\..*", ".\\1", LocationSignals$ProtID)
LocationSignals <- LocationSignals %>%
  mutate(ProtID = sub("\\.(\\d+)\\..*", ".\\1", ProtID)) %>%  # Extract numeric part
  distinct(ProtID, .keep_all = TRUE)  # Remove duplicated rows based on ProtID
names(LocationSignals)[names(LocationSignals) == "ProtID"] <- "geneID"
LocationSignals <- LocationSignals[, c(1, 2)]
names(LocationSignals)[names(LocationSignals) == "Type"] <- "LocationSignals_Type"


# Pennance et al. 2024 Supp Table 14 CREPs/FREPs
FREPsCREPs <- read.csv("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_old/Supp_Table_14_Bsud_CREP_FREP_Summary.csv", header = TRUE)
FREPsCREPs$ProtID <- sub("\\.(\\d+)\\..*", ".\\1", FREPsCREPs$ProtID)
names(FREPsCREPs)[names(FREPsCREPs) == "ProtID"] <- "geneID"
FREPsCREPs <- FREPsCREPs[, c(1, 3)]

# Pennance et al. 2024 Supp Table 17 Hyperdiverse genes
DiverseGenes <- read.csv("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_old/Supp_Table_17_Hyperdiverse_genes.csv", header = TRUE)
DiverseGenes <- DiverseGenes[, c(1, 9, 10, 12, 13, 14, 32)]

# Lu et al. 2022 Table S3. Expression summary of the previous genome wide mapping studies (GWMS) and other 6 feature-specific studies 
LuTabS3 <- read.csv("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_old/Lu-2022-Table_S3.csv", header = TRUE)
names(LuTabS3)[names(LuTabS3) == "VB.ID"] <- "BglabBB02_geneID"

# Linkage groups from B. sudanica
LGsBsud <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/Biom_sud_annotation/Linkage_maps/Bsud111_ccs_assembly_contig_lengths_Bu_ordered_tweaked.txt", header = TRUE, sep = "\t")
LGsBsud <- LGsBsud[, c(1,4)]

#### Interpro data ####
Supp_File_3_Interpro <- read_xlsx("Supp_File_3_Interpro.xlsx")

# Step 1: Rename "Protein_ID" to "geneID"
Supp_File_3_Interpro <- Supp_File_3_Interpro %>%
  rename(geneID = Protein_ID)

# Step 2: Remove everything after the second '.' in "geneID"
Supp_File_3_Interpro <- Supp_File_3_Interpro %>%
  mutate(geneID = str_replace(geneID, "^([^.]+\\.[^.]+)\\..*", "\\1"))

# Step 3: Filter rows with non-missing InterPro_Accession and create Supp_File_3_Interpro_2
Supp_File_3_Interpro_2 <- Supp_File_3_Interpro %>%
  filter(
    !is.na(InterPro_Accession) &
      geneID != "-" &
      InterPro_Accession != "-" &
      InterPro_Accession_description != "-"
  ) %>%
  dplyr::select(geneID, InterPro_Accession, InterPro_Accession_description)

# Step 4: Concatenate unique InterPro_Accession and InterPro_Accession_description per geneID
Supp_File_3_Interpro_2 <- Supp_File_3_Interpro_2 %>%
  group_by(geneID) %>%
  summarise(
    InterPro_Accession = paste(unique(InterPro_Accession), collapse = ","),
    InterPro_Accession_description = paste(unique(InterPro_Accession_description), collapse = ";")
  ) %>%
  ungroup()

rm(Supp_File_3_Interpro)


### MERGE DATA ####

# Merge dataframes
temp_merge_result_df <- merge(temp_merge_result_df, Annotation_ALL_DE_pwcomp_p0.05, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, LocationSignals, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, New_TMR_Table, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, KeyGenes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, Bsud_BB02_orthologs, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, Bsud_Bp505_orthologs, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, FREPsCREPs, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, DiverseGenes, by = "geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, LuTabS3, by = "BglabBB02_geneID", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, LGsBsud, by = "Contig", all.x = TRUE)
temp_merge_result_df <- merge(temp_merge_result_df, Supp_File_3_Interpro_2, by = "geneID", all.x = TRUE)


#### CHECK: Remove columns that are not significant DE genes for CvsN or CvsK (shouldn't be any) ####
# Subset the dataframe to keep rows where at least one of the specified columns has a value less than 0.05
# check colnames
colnames(temp_merge_result_df)[startsWith(colnames(temp_merge_result_df), "FDR")]
temp_merge_result_df <- temp_merge_result_df[rowSums(temp_merge_result_df[, c("FDR_C8vsK8", "FDR_C8vsK8Sm", "FDR_C24vsK24", "FDR_C24vsK24Sm", "FDR_C72vsK72", "FDR_C72vsK72Sm", "FDR_C8vsN8", "FDR_C8vsN8Sm", "FDR_C24vsN24", "FDR_C72vsN72")] < 0.05, na.rm = TRUE) > 0, ]
# should stay the same size in our data.


#### Add interpretation columns ####

## Add the new 'DE_result_CvNMRI' column

# Add if significantly up or down regulated at 8 hours 
temp_merge_result_df <- mutate(temp_merge_result_df,
                            sig_DE_N8 = case_when(
                              logFC_C8vsN8 < 0 & FDR_C8vsN8 < 0.05 ~ "Up",
                              logFC_C8vsN8 > 0 & FDR_C8vsN8 < 0.05 ~ "Down",
                              TRUE ~ "Not_Sig"
                            ))

temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_N8Sm = case_when(
                                 logFC_C8vsN8Sm < 0 & FDR_C8vsN8Sm < 0.05 ~ "Up",
                                 logFC_C8vsN8Sm > 0 & FDR_C8vsN8Sm < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

# Add if significantly up or down regulated at 24 hours 
temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_N24 = case_when(
                                 logFC_C24vsN24 < 0 & FDR_C24vsN24 < 0.05 ~ "Up",
                                 logFC_C24vsN24 > 0 & FDR_C24vsN24 < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

# Add if significantly up or down regulated at 72 hours
temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_N72 = case_when(
                                 logFC_C72vsN72 < 0 & FDR_C72vsN72 < 0.05 ~ "Up",
                                 logFC_C72vsN72 > 0 & FDR_C72vsN72 < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

# Add the new 'DE_result_CvUNMK' column

# Add if significantly up or down regulated at 8 hours 
temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K8 = case_when(
                                 logFC_C8vsK8 < 0 & FDR_C8vsK8 < 0.05 ~ "Up",
                                 logFC_C8vsK8 > 0 & FDR_C8vsK8 < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K8Sm = case_when(
                                 logFC_C8vsK8Sm < 0 & FDR_C8vsK8Sm < 0.05 ~ "Up",
                                 logFC_C8vsK8Sm > 0 & FDR_C8vsK8Sm < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

# Add if significantly up or down regulated at 24 hours 
temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K24 = case_when(
                                 logFC_C24vsK24 < 0 & FDR_C24vsK24 < 0.05 ~ "Up",
                                 logFC_C24vsK24 > 0 & FDR_C24vsK24 < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K24Sm = case_when(
                                 logFC_C24vsK24Sm < 0 & FDR_C24vsK24Sm < 0.05 ~ "Up",
                                 logFC_C24vsK24Sm > 0 & FDR_C24vsK24Sm < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

# Add if significantly up or down regulated at 72 hours
temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K72 = case_when(
                                 logFC_C72vsK72 < 0 & FDR_C72vsK72 < 0.05 ~ "Up",
                                 logFC_C72vsK72 > 0 & FDR_C72vsK72 < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))

temp_merge_result_df <- mutate(temp_merge_result_df,
                               sig_DE_K72Sm = case_when(
                                 logFC_C72vsK72Sm < 0 & FDR_C72vsK72Sm < 0.05 ~ "Up",
                                 logFC_C72vsK72Sm > 0 & FDR_C72vsK72Sm < 0.05 ~ "Down",
                                 TRUE ~ "Not_Sig"
                               ))


# Move the new columns to the first positions
temp_merge_result_df <- dplyr::select(temp_merge_result_df, geneID, Contig, LG, Start, End, Length_aa, N_TMRs, p0.001, sig_DE_N8, sig_DE_N8Sm, sig_DE_N24, sig_DE_N72, sig_DE_K8, sig_DE_K8Sm, sig_DE_K24, sig_DE_K24Sm, sig_DE_K72, sig_DE_K72Sm, everything())

write.table(temp_merge_result_df, file = "BSUD_transcriptomics_sig_DE_genes_edgeR_v4.8.2.tab", sep = "\t", col.names = TRUE, row.names = FALSE)
write_xlsx(temp_merge_result_df, path = "BSUD_transcriptomics_sig_DE_genes_edgeR_v4.8.2.xlsx")
