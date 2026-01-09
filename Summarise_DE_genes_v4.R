# Script to summarise data of DE genes generated in script Annotate_DE_genes_vX.R
library(ggplot2)
library(tidyverse)
library(dplyr)
#install.packages("ggrepel")  # Install ggrepel if you haven't already
library(ggrepel)
setwd("/Users/tpennance/Documents/tpennanceWU/R_folder/mRNA/mRNA_analysis_v2024/")

bsudDEdata <- read.table("BSUD_transcriptomics_sig_DE_genes_edgeR_v4.8.2.tab", header = TRUE, sep = "\t")


#### Summary Tables ####



####  GENES LogFC PW COMPARISONS ####
if (!dir.exists("UP_n_DOWN_regulated_pw")) {
  dir.create("UP_n_DOWN_regulated_pw")
  message("Directory 'UP_n_DOWN_regulated_pw' created.")
} else {
  message("Directory 'UP_n_DOWN_regulated_pw' already exists.")
}

### C8vN8 ###

# Select  genes with the lowest values in logFC_C8vsN8
C8vsN8_upregulated <- bsudDEdata %>%
  filter(logFC_C8vsN8 < 0) %>%
  filter(FDR_C8vsN8 < 0.05) %>%
  arrange(logFC_C8vsN8)
write.table(C8vsN8_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsN8_upegulated.tab", sep = "\t", col.names = NA)

# Select  genes with the highest values in logFC_C8vsN8
C8vsN8_downregulated <- bsudDEdata %>%
  filter(logFC_C8vsN8 > 0) %>%
  filter(FDR_C8vsN8 < 0.05) %>%
  arrange(desc(logFC_C8vsN8))
write.table(C8vsN8_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsN8_downregulated.tab", sep = "\t", col.names = NA)


### C8vK8 ###

C8vsK8_upregulated <- bsudDEdata %>%
  filter(logFC_C8vsK8 < 0) %>%
  filter(FDR_C8vsK8 < 0.05) %>%
  arrange(logFC_C8vsK8)
write.table(C8vsK8_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsK8_upegulated.tab", sep = "\t", col.names = NA)

# Select  genes with the highest values in logFC_C8vsN8
C8vsK8_downregulated <- bsudDEdata %>%
  filter(logFC_C8vsK8 > 0) %>%
  filter(FDR_C8vsK8 < 0.05) %>%
  arrange(desc(logFC_C8vsK8))
write.table(C8vsN8_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsK8_downregulated.tab", sep = "\t", col.names = NA)

### C8vK8Sm ###

C8vsK8Sm_upregulated <- bsudDEdata %>%
  filter(logFC_C8vsK8Sm < 0) %>%
  filter(FDR_C8vsK8Sm < 0.05) %>%
  arrange(logFC_C8vsK8Sm)
write.table(C8vsK8Sm_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsK8Sm_upegulated.tab", sep = "\t", col.names = NA)

C8vsK8Sm_downregulated <- bsudDEdata %>%
  filter(logFC_C8vsK8Sm > 0) %>%
  filter(FDR_C8vsK8Sm < 0.05) %>%
  arrange(desc(logFC_C8vsK8Sm))
write.table(C8vsN8_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C8vsK8Sm_downregulated.tab", sep = "\t", col.names = NA)




### C24vN24 ###

# Select  genes with the lowest values (upregulated) in logFC_C24vsN24
C24vsN24_upregulated <- bsudDEdata %>%
  filter(logFC_C24vsN24 < 0) %>%
  filter(FDR_C24vsN24 < 0.05) %>%
  arrange(logFC_C24vsN24)
write.table(C24vsN24_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsN24_upregulated.tab", sep = "\t", col.names = NA)

# Select  genes with the highest values (downregulated) in logFC_C24vsN24
C24vsN24_downregulated <- bsudDEdata %>%
  filter(logFC_C24vsN24 > 0) %>%
  filter(FDR_C24vsN24 < 0.05) %>%
  arrange(desc(logFC_C24vsN24))
write.table(C24vsN24_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsN24_downregulated.tab", sep = "\t", col.names = NA)


### C24vK24 ###

# Select  genes with the lowest values (upregulated) in logFC_C24vsK24
C24vsK24_upregulated <- bsudDEdata %>%
  filter(logFC_C24vsK24 < 0) %>%
  filter(FDR_C24vsK24 < 0.05) %>%
  arrange(logFC_C24vsK24)
write.table(C24vsK24_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsK24_upregulated.tab", sep = "\t", col.names = NA)

# Select  genes with the highest values (downregulated) in logFC_C24vsK24
C24vsK24_downregulated <- bsudDEdata %>%
  filter(logFC_C24vsK24 > 0) %>%
  filter(FDR_C24vsK24 < 0.05) %>%
  arrange(desc(logFC_C24vsK24))
write.table(C24vsK24_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsK24_downregulated.tab", sep = "\t", col.names = NA)

### C24vK24Sm ###

# Select  genes with the lowest values (upregulated) in logFC_C24vsK24Sm
C24vsK24Sm_upregulated <- bsudDEdata %>%
  filter(logFC_C24vsK24Sm < 0) %>%
  filter(FDR_C24vsK24Sm < 0.05) %>%
  arrange(logFC_C24vsK24Sm)
write.table(C24vsK24Sm_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsK24Sm_upregulated.tab", sep = "\t", col.names = NA)

# Select  genes with the highest values (downregulated) in logFC_C24vsK24Sm
C24vsK24Sm_downregulated <- bsudDEdata %>%
  filter(logFC_C24vsK24Sm > 0) %>%
  filter(FDR_C24vsK24Sm < 0.05) %>%
  arrange(desc(logFC_C24vsK24Sm))
write.table(C24vsK24Sm_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C24vsK24Sm_downregulated.tab", sep = "\t", col.names = NA)


### C72vN72 ###

# Select  genes with the lowest values (upregulated) in logFC_C72vsN72
C72vsN72_upregulated <- bsudDEdata %>%
  filter(logFC_C72vsN72 < 0) %>%
  filter(FDR_C72vsN72 < 0.05) %>%
  arrange(logFC_C72vsN72)
write.table(C72vsN72_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsN72_upregulated.tab", sep = "\t", col.names = NA)


# Select  genes with the highest values (downregulated) in logFC_C72vsN72
C72vsN72_downregulated <- bsudDEdata %>%
  filter(logFC_C72vsN72 > 0) %>%
  filter(FDR_C72vsN72 < 0.05) %>%
  arrange(desc(logFC_C72vsN72))
write.table(C72vsN72_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsN72_downregulated.tab", sep = "\t", col.names = NA)


### C72vK72 ###

# Select  genes with the lowest values (upregulated) in logFC_C72vsK72
C72vsK72_upregulated <- bsudDEdata %>%
  filter(logFC_C72vsK72 < 0) %>%
  filter(FDR_C72vsK72 < 0.05) %>%
  arrange(logFC_C72vsK72)
write.table(C72vsK72_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsK72_upregulated.tab", sep = "\t", col.names = NA)

# Select genes with the highest values (downregulated) in logFC_C72vsK72
C72vsK72_downregulated <- bsudDEdata %>%
  filter(logFC_C72vsK72 > 0) %>%
  filter(FDR_C72vsK72 < 0.05) %>%
  arrange(desc(logFC_C72vsK72))
write.table(C72vsK72_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsK72_downregulated.tab", sep = "\t", col.names = NA)

### C72vK72Sm ###

# Select  genes with the lowest values (upregulated) in logFC_C72vsK72Sm
C72vsK72Sm_upregulated <- bsudDEdata %>%
  filter(logFC_C72vsK72Sm < 0) %>%
  filter(FDR_C72vsK72Sm < 0.05) %>%
  arrange(logFC_C72vsK72Sm)
write.table(C72vsK72Sm_upregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsK72Sm_upregulated.tab", sep = "\t", col.names = NA)

# Select genes with the highest values (downregulated) in logFC_C72vsK72Sm
C72vsK72Sm_downregulated <- bsudDEdata %>%
  filter(logFC_C72vsK72Sm > 0) %>%
  filter(FDR_C72vsK72Sm < 0.05) %>%
  arrange(desc(logFC_C72vsK72Sm))
write.table(C72vsK72Sm_downregulated, file= "UP_n_DOWN_regulated_pw/logFC_C72vsK72Sm_downregulated.tab", sep = "\t", col.names = NA)





#### VOLCANO PLOTS OF DE GENES FDR<0.05 ####

if (!dir.exists("volcano_plots")) {
  dir.create("volcano_plots")
  message("Directory 'volcano_plots' created.")
} else {
  message("Directory 'volcano_plots' already exists.")
}

C8vsK8_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsK8_2024_11_11_test.tab", header = TRUE, sep = "\t")
C8vsK8Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsK8Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
C8vsN8_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsN8_2024_11_11_test.tab", header = TRUE, sep = "\t")
C8vsN8Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C8vsN8Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
C24vsK24_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsK24_2024_11_11_test.tab", header = TRUE, sep = "\t")
C24vsK24Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsK24Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
C24vsN24_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C24vsN24_2024_11_11_test.tab", header = TRUE, sep = "\t")
C72vsK72_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsK72_2024_11_11_test.tab", header = TRUE, sep = "\t")
C72vsK72Sm_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsK72Sm_2024_11_11_test.tab", header = TRUE, sep = "\t")
C72vsN72_all_genes <- read.table("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/PW_DE_tables/Results_C72vsN72_2024_11_11_test.tab", header = TRUE, sep = "\t")

# volcano plot for C8vK8 all genes
pdf("volcano_plots/C8vsK8_volcano_plot.pdf")
ggplot(C8vsK8_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsK8_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs K8)", y = "-log10(FDR)", title = "C8vsK8 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C8vsK8_all_genes[order(C8vsK8_all_genes$FDR), ], 5),
                  aes(label = X), size = 3) +
  #geom_text_repel(data = C8vsK8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C8vsK8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C8vK8Sm all genes
pdf("volcano_plots/C8vsK8Sm_volcano_plot.pdf")
ggplot(C8vsK8Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsK8Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs K8Sm)", y = "-log10(FDR)", title = "C8vsK8Sm DE") +
  theme_minimal() +
  geom_text_repel(data = head(C8vsK8Sm_all_genes[order(C8vsK8Sm_all_genes$FDR), ], 5),
                  aes(label = X), size = 3) +
  #geom_text_repel(data = C8vsK8Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C8vsK8Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C8vN8 all genes
pdf("volcano_plots/C8vsN8_volcano_plot.pdf")
ggplot(C8vsN8_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsN8_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs N8)", y = "-log10(FDR)", title = "C8vsN8 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C8vsN8_all_genes[order(C8vsN8_all_genes$FDR), ], 4),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C8vsN8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  geom_text_repel(data = C8vsN8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C8vN8Sm all genes
#NA because no significant differences or enough comparisons (only one N8Sm sample)

# volcano plot for C24vK24 all genes
pdf("volcano_plots/C24vsK24_volcano_plot.pdf")
ggplot(C24vsK24_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsK24_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs K24)", y = "-log10(FDR)", title = "C24vsK24 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C24vsK24_all_genes[order(C24vsK24_all_genes$FDR), ], 10),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C24vsK24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  geom_text_repel(data = C24vsK24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C24vK24Sm all genes
pdf("volcano_plots/C24vsK24Sm_volcano_plot.pdf")
ggplot(C24vsK24Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsK24Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs K24Sm)", y = "-log10(FDR)", title = "C24vsK24Sm DE") +
  theme_minimal() +
  geom_text_repel(data = head(C24vsK24Sm_all_genes[order(C24vsK24Sm_all_genes$FDR), ], 10),
                  aes(label = X), size = 3) +
  #geom_text_repel(data = C24vsK24Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(8, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C24vsK24Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-2, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C24vN24 all genes
pdf("volcano_plots/C24vsN24_volcano_plot.pdf")
ggplot(C24vsN24_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsN24_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs N24)", y = "-log10(FDR)", title = "C24vsN24 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C24vsN24_all_genes[order(C24vsN24_all_genes$FDR), ], 15),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C24vsN24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C24vsN24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vK72 all genes
pdf("volcano_plots/C72vsK72_volcano_plot.pdf")
ggplot(C72vsK72_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsK72_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs K72)", y = "-log10(FDR)", title = "C72vsK72 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C72vsK72_all_genes[order(C72vsK72_all_genes$FDR), ], 4),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C72vsK72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  geom_text_repel(data = C72vsK72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vK72Sm all genes
pdf("volcano_plots/C72vsK72Sm_volcano_plot.pdf")
ggplot(C72vsK72Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsK72Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs K72Sm)", y = "-log10(FDR)", title = "C72vsK72Sm DE") +
  theme_minimal() +
  geom_text_repel(data = head(C72vsK72Sm_all_genes[order(C72vsK72Sm_all_genes$FDR), ], 4),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C72vsK72Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  geom_text_repel(data = C72vsK72Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vN72 all genes
pdf("volcano_plots/C72vsN72_volcano_plot.pdf")
ggplot(C72vsN72_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsN72_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs N72)", y = "-log10(FDR)", title = "C72vsN72 DE") +
  theme_minimal() +
  geom_text_repel(data = head(C72vsN72_all_genes[order(C72vsN72_all_genes$FDR), ], 4),
                  aes(label = X), size = 3) +
  geom_text_repel(data = C72vsN72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  geom_text_repel(data = C72vsN72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
                  aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  guides(color = "none")  # Remove legend
dev.off()


#### scaled volcano plots ####

# Define axis limits
x_limits <- c(-12, 12)
y_limits <- c(0, 11)

# Volcano plot for C8 vs K8
pdf("volcano_plots/C8vsK8_volcano_plot_scaled.pdf")
ggplot(C8vsK8_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsK8_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs K8)", y = "-log10(FDR)", title = "C8vsK8 DE") +
  theme_minimal() +
  xlim(x_limits) + ylim(y_limits) +  # Set fixed axis limits
  #geom_text_repel(data = head(C8vsK8_all_genes[order(C8vsK8_all_genes$FDR), ], 5),
  #                aes(label = X), size = 3) +
  guides(color = "none")  # Remove legend
dev.off()

# Volcano plot for C8 vs K8Sm
pdf("volcano_plots/C8vsK8Sm_volcano_plot_scaled.pdf")
ggplot(C8vsK8Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsK8Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs K8Sm)", y = "-log10(FDR)", title = "C8vsK8Sm DE") +
  theme_minimal() +
  xlim(x_limits) + ylim(y_limits) +  # Set fixed axis limits
  #geom_text_repel(data = head(C8vsK8Sm_all_genes[order(C8vsK8Sm_all_genes$FDR), ], 5),
  #                aes(label = X), size = 3) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C8vN8 all genes
pdf("volcano_plots/C8vsN8_volcano_plot_scaled.pdf")
ggplot(C8vsN8_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C8vsN8_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C8 vs N8)", y = "-log10(FDR)", title = "C8vsN8 DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C8vsN8_all_genes[order(C8vsN8_all_genes$FDR), ], 4),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C8vsN8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C8vsN8_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C8vN8Sm all genes
#NA because no significant differences or enough comparisons (only one N8Sm sample)

# volcano plot for C24vK24 all genes
pdf("volcano_plots/C24vsK24_volcano_plot_scaled.pdf")
ggplot(C24vsK24_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsK24_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs K24)", y = "-log10(FDR)", title = "C24vsK24 DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C24vsK24_all_genes[order(C24vsK24_all_genes$FDR), ], 10),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C24vsK24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C24vsK24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C24vK24Sm all genes
pdf("volcano_plots/C24vsK24Sm_volcano_plot_scaled.pdf")
ggplot(C24vsK24Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsK24Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs K24Sm)", y = "-log10(FDR)", title = "C24vsK24Sm DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C24vsK24Sm_all_genes[order(C24vsK24Sm_all_genes$FDR), ], 10),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C24vsK24Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(8, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C24vsK24Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-2, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C24vN24 all genes
pdf("volcano_plots/C24vsN24_volcano_plot_scaled.pdf")
ggplot(C24vsN24_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C24vsN24_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C24 vs N24)", y = "-log10(FDR)", title = "C24vsN24 DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C24vsN24_all_genes[order(C24vsN24_all_genes$FDR), ], 15),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C24vsN24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C24vsN24_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vK72 all genes
pdf("volcano_plots/C72vsK72_volcano_plot_scaled.pdf")
ggplot(C72vsK72_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsK72_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs K72)", y = "-log10(FDR)", title = "C72vsK72 DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C72vsK72_all_genes[order(C72vsK72_all_genes$FDR), ], 4),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C72vsK72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C72vsK72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vK72Sm all genes
pdf("volcano_plots/C72vsK72Sm_volcano_plot_scaled.pdf")
ggplot(C72vsK72Sm_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsK72Sm_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs K72Sm)", y = "-log10(FDR)", title = "C72vsK72Sm DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C72vsK72Sm_all_genes[order(C72vsK72Sm_all_genes$FDR), ], 4),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C72vsK72Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C72vsK72Sm_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

# volcano plot for C72vN72 all genes
pdf("volcano_plots/C72vsN72_volcano_plot_scaled.pdf")
ggplot(C72vsN72_all_genes, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(C72vsN72_all_genes$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Add FDR threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +  # Add logFC threshold lines
  labs(x = "logFC (C72 vs N72)", y = "-log10(FDR)", title = "C72vsN72 DE") +
  theme_minimal() +
  #geom_text_repel(data = head(C72vsN72_all_genes[order(C72vsN72_all_genes$FDR), ], 4),
  #                aes(label = X), size = 3) +
  #geom_text_repel(data = C72vsN72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(4, logFC) %>% mutate(label_type = "high_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  #geom_text_repel(data = C72vsN72_all_genes %>% filter(-log10(FDR) > -log10(0.05)) %>% top_n(-4, logFC) %>% mutate(label_type = "low_logFC"),
  #                aes(label = X, color = label_type), size = 3) +
  scale_color_manual(values = c("high_logFC" = "blue", "low_logFC" = "green")) +
  xlim(x_limits) + ylim(y_limits) +
  guides(color = "none")  # Remove legend
dev.off()

#### Line graphs for logFC change in particular genes over time ####

if (!dir.exists("logFC_time_change_plots")) {
  dir.create("logFC_time_change_plots")
  message("Directory 'logFC_time_change_plots' created.")
} else {
  message("Directory 'logFC_time_change_plots' already exists.")
}

## C8vK8 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C8vsK8_upregulated$geneID[1:4]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsK8_upegulated_top4.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 4 Genes upregulated by logFC at C8vsK8") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C8vsK8_downregulated$geneID[1:1]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsK8_downregulated_BSUD.21951.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "BSUD.21951 downregulated by logFC at C8vsK8") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()

## C24vK24 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C24vsK24_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsK24_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C24vsK24") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C24vsK24_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsK24_downregulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C24vsK24") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()


## C72vK72 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C72vsK72_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsK72_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C72vsK72") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C72vsK72_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsK72_downregulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C72vsK72") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()




## C8vN8 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
# BSUD.3525 not listed in the 
top_10_genes <- C8vsN8_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
 
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsN8_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C8vsN8") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()




# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C8vsN8_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsN8_downregulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C8vsN8") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()




## C24vN24 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
# BSUD.3525 not listed in the 
top_10_genes <- C24vsN24_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsN24_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C24vsN24") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
# BSUD.3525 not listed in the 
top_10_genes <- C24vsN24_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsN24_downegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C24vsN24") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()




## C72vN72 ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C72vsN72_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsN72_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C72vsN72") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()




# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C72vsN72_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsN8_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsN24_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsN72_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsN8_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsN24_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsN72_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsN72_downegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C72vsN72") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



## C8vK8Sm ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C8vsK8Sm_upregulated$geneID[1:1]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsK8Sm_upegulated_BSUD.5984.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "BSUD.5984 upregulated by logFC at C8vsK8Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C8vsK8Sm_downregulated$geneID[1:3]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C8vsK8Sm_downregulated_top_3.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 3 genes downregulated by logFC at C8vsK8Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()

## C24vK24Sm ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C24vsK24Sm_upregulated$geneID[1:8]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsK24Sm_upegulated_top8.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C24vsK24Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C24vsK24Sm_downregulated$geneID[1:2]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C24vsK24Sm_downregulated_top2.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 2 Genes downregulated by logFC at C24vsK24Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()


## C72vK72Sm ##

# Select top 10 genes that are SIGNIFICANTLY upregulated FDR p<0.05
top_10_genes <- C72vsK72Sm_upregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsK72Sm_upegulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes upregulated by logFC at C72vsK72Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()



# Select top 10 genes that are SIGNIFICANTLY DOWNregulated FDR p<0.05
top_10_genes <- C72vsK72Sm_downregulated$geneID[1:10]
# Create a color palette for the genes
gene_colors <- rainbow(length(top_10_genes))
# Initialize an empty dataframe to store the data for plotting
plot_data <- data.frame()

# Iterate over each gene to extract logFC and FDR values from respective dataframes and create plot data
for (i in 1:length(top_10_genes)) {
  gene_id <- top_10_genes[i]
  
  # Extract logFC and FDR values for the gene at each time point
  logFC_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$logFC
  logFC_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$logFC
  logFC_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$logFC
  FDR_8 <- subset(C8vsK8Sm_all_genes, X == gene_id)$FDR
  FDR_24 <- subset(C24vsK24Sm_all_genes, X == gene_id)$FDR
  FDR_72 <- subset(C72vsK72Sm_all_genes, X == gene_id)$FDR
  
  # Combine logFC and FDR values into a single dataframe
  gene_data <- data.frame(Time = factor(c("8 hours", "24 hours", "72 hours"), levels = c("8 hours", "24 hours", "72 hours")),
                          logFC = c(logFC_8, logFC_24, logFC_72),
                          FDR_0.05 = c(FDR_8, FDR_24, FDR_72) < 0.05,  # Convert to boolean
                          Gene_ID = gene_id,
                          Color = gene_colors[i])
  
  # Append gene data to plot_data
  plot_data <- rbind(plot_data, gene_data)
}

# Plot the line plot
pdf("logFC_time_change_plots/C72vsK72Sm_downregulated_top10.pdf")
ggplot(plot_data, aes(x = Time, y = logFC, color = Gene_ID, group = Gene_ID)) +
  geom_line() +
  geom_point(aes(shape = FDR_0.05), size = 3) +  # Differentiate points based on significance
  scale_shape_manual(values = c(4, 16)) +  # 16: Circle (significant), 4: Cross (insignificant)
  scale_color_manual(values = gene_colors) +
  labs(x = "Time", y = "logFC", title = "Top 10 Genes downregulated by logFC at C72vsK72Sm") +
  theme_bw() +  # Remove grey background
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "dotted"))  # Insert y-axis grid lines
dev.off()


#### VENN diagrams gene similarities ####

if (!dir.exists("venn_diagrams")) {
  dir.create("venn_diagrams")
  message("Directory 'venn_diagrams' created.")
} else {
  message("Directory 'venn_diagrams' already exists.")
}

#setwd("venn_diagrams")
install.packages("VennDiagram")
library(VennDiagram)


# Define lists of upregulated and downregulated genes for each comparison

# For 8-hour comparison (sig_DE_N8 vs sig_DE_K8)
upregulated_N8 <- bsudDEdata %>% filter(sig_DE_N8 == "Up") %>% pull(geneID)
upregulated_K8 <- bsudDEdata %>% filter(sig_DE_K8 == "Up") %>% pull(geneID)
downregulated_N8 <- bsudDEdata %>% filter(sig_DE_N8 == "Down") %>% pull(geneID)
downregulated_K8 <- bsudDEdata %>% filter(sig_DE_K8 == "Down") %>% pull(geneID)

# For 24-hour comparison (sig_DE_N24 vs sig_DE_K24)
upregulated_N24 <- bsudDEdata %>% filter(sig_DE_N24 == "Up") %>% pull(geneID)
upregulated_K24 <- bsudDEdata %>% filter(sig_DE_K24 == "Up") %>% pull(geneID)
downregulated_N24 <- bsudDEdata %>% filter(sig_DE_N24 == "Down") %>% pull(geneID)
downregulated_K24 <- bsudDEdata %>% filter(sig_DE_K24 == "Down") %>% pull(geneID)

# For 72-hour comparison (sig_DE_N72 vs sig_DE_K72)
upregulated_N72 <- bsudDEdata %>% filter(sig_DE_N72 == "Up") %>% pull(geneID)
upregulated_K72 <- bsudDEdata %>% filter(sig_DE_K72 == "Up") %>% pull(geneID)
downregulated_N72 <- bsudDEdata %>% filter(sig_DE_N72 == "Down") %>% pull(geneID)
downregulated_K72 <- bsudDEdata %>% filter(sig_DE_K72 == "Down") %>% pull(geneID)

# Generate Venn Diagrams

# Upregulated genes at 8 hours
pdf("venn_diagrams/venn.plot.up_8.pdf")
venn.plot.up_8 <- venn.diagram(
  x = list(N8 = upregulated_N8, K8 = upregulated_K8),
  category.names = c("sig_DE_N8", "sig_DE_K8"),
  filename = NULL,
  main = "Upregulated Genes (8 Hours)"
)
grid.draw(venn.plot.up_8)
dev.off()

# Downregulated genes at 8 hours
#pdf("venn_diagrams/venn.plot.down_8.pdf")
#venn.plot.down_8 <- venn.diagram(
##  x = list(N8 = downregulated_N8, K8 = downregulated_K8),
#  category.names = c("sig_DE_N8", "sig_DE_K8"),
#  filename = NULL,
#  main = "Downregulated Genes (8 Hours)"
#)
#grid.draw(venn.plot.down_8)
#dev.off

# Repeat the above process for 24 and 72 hours

# Upregulated genes at 24 hours
pdf("venn_diagrams/venn.plot.up_24.pdf")
venn.plot.up_24 <- venn.diagram(
  x = list(N24 = upregulated_N24, K24 = upregulated_K24),
  category.names = c("sig_DE_N24", "sig_DE_K24"),
  filename = NULL,
  main = "Upregulated Genes (24 Hours)"
)
grid.draw(venn.plot.up_24)
dev.off()

# Downregulated genes at 24 hours
pdf("venn_diagrams/venn.plot.down_24.pdf")
venn.plot.down_24 <- venn.diagram(
  x = list(N24 = downregulated_N24, K24 = downregulated_K24),
  category.names = c("sig_DE_N24", "sig_DE_K24"),
  filename = NULL,
  main = "Downregulated Genes (24 Hours)"
)
grid.draw(venn.plot.down_24)
dev.off()

# Upregulated genes at 72 hours
pdf("venn_diagrams/venn.plot.up_72.pdf")
venn.plot.up_72 <- venn.diagram(
  x = list(N72 = upregulated_N72, K72 = upregulated_K72),
  category.names = c("sig_DE_N72", "sig_DE_K72"),
  filename = NULL,
  main = "Upregulated Genes (72 Hours)"
)
grid.draw(venn.plot.up_72)
dev.off()

# Downregulated genes at 72 hours
pdf("venn_diagrams/venn.plot.down_72.pdf")
venn.plot.down_72 <- venn.diagram(
  x = list(N72 = downregulated_N72, K72 = downregulated_K72),
  category.names = c("sig_DE_N72", "sig_DE_K72"),
  filename = NULL,
  main = "Downregulated Genes (72 Hours)"
)
grid.draw(venn.plot.down_72)
dev.off()

rm(venn.plot.up_8, venn.plot.down_8, venn.plot.up_24, venn.plot.down_24, venn.plot.up_72, venn.plot.down_72)
# Remove all files starting with "upregulated_" or "downregulated_"
file.remove(list.files(pattern = "upregulated_*"))
file.remove(list.files(pattern = "downregulated_.*"))



#### Histogram of up and down genes over comparison ####

if (!dir.exists("histogram")) {
  dir.create("histogram")
  message("Directory 'histogram' created.")
} else {
  message("Directory 'histogram' already exists.")
}

# Assuming columns: geneID, logFC_C8vsK8, logFC_C8vsN8, logFC_C24vsK24, logFC_C24vsN24, logFC_C72vsK72, logFC_C72vsN72

# Define a function to categorize logFC values into bins
categorize_logFC <- function(logFC) {
  case_when(
    logFC >= 4 ~ ">4",
    logFC >= 2 & logFC < 4 ~ "2-4",
    logFC >= 1 & logFC < 2 ~ "1-2",
    logFC > 0 & logFC < 1 ~ "0-1",
    logFC <= -4 ~ "<-4",
    logFC <= -2 & logFC > -4 ~ "-4 to -2",
    logFC <= -1 & logFC > -2 ~ "-2 to -1",
    logFC < 0 & logFC > -1 ~ "-1 to 0",
    TRUE ~ "Not_Sig"
  )
}

# Apply categorization for each comparison
bsudDEdata_long <- bsudDEdata %>%
  pivot_longer(cols = starts_with("logFC"), names_to = "Comparison", values_to = "logFC") %>%
  mutate(Category = categorize_logFC(logFC)) %>%
  filter(Category != "Not_Sig")

# Count the number of genes in each bin for each comparison
gene_counts <- bsudDEdata_long %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  ungroup()

# Plotting the bar chart

pdf("histogram/histogram.pdf")
ggplot(gene_counts, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "blue", "lightblue", "lightgreen", "pink", "red", "darkred", "black")) +
  labs(
    title = "Distribution of Upregulated and Downregulated Genes by LogFC Range",
    x = "Comparison",
    y = "Number of Genes",
    fill = "LogFC Range"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#### histogram version 2 ####

# Prepare data by filtering out non-significant genes
bsudDEdata_filtered <- bsudDEdata %>%
  filter(sig_DE_N8 != "Not_Sig" | sig_DE_K8 != "Not_Sig" |
           sig_DE_N24 != "Not_Sig" | sig_DE_K24 != "Not_Sig" |
           sig_DE_N72 != "Not_Sig" | sig_DE_K72 != "Not_Sig") %>%
  pivot_longer(
    cols = starts_with("logFC"),
    names_to = "Comparison",
    values_to = "logFC"
  ) %>%
  mutate(Category = categorize_logFC(logFC)) %>%
  filter(Category != "Not_Sig")

# Filter out "Not_Sig" genes based on corresponding sig_DE... column
bsudDEdata_filtered <- bsudDEdata_filtered %>%
  filter((Comparison == "logFC_C8vsK8" & sig_DE_K8 != "Not_Sig") |
           (Comparison == "logFC_C8vsN8" & sig_DE_N8 != "Not_Sig") |
           (Comparison == "logFC_C24vsK24" & sig_DE_K24 != "Not_Sig") |
           (Comparison == "logFC_C24vsN24" & sig_DE_N24 != "Not_Sig") |
           (Comparison == "logFC_C72vsK72" & sig_DE_K72 != "Not_Sig") |
           (Comparison == "logFC_C72vsN72" & sig_DE_N72 != "Not_Sig"))

# Count the number of genes in each bin for each comparison
gene_counts <- bsudDEdata_filtered %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  ungroup()

# Plotting the bar chart
pdf("histogram/histogram2.pdf")
ggplot(gene_counts, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of significantly Upregulated and Downregulated Genes by LogFC Range",
    x = "Comparison",
    y = "Number of Genes",
    fill = "LogFC Range"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
