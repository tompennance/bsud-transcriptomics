setwd("/Users/tpennance/Documents/tpennanceWU/PostDoc/RESEARCH_PROPOSALS/Transcriptomics/analysis/EdgeR_results_2024/")
# Analysis of differential expression between the different exposures of B. sudanica (KEMRI) to S. mansoni (UNMKenya v NMRI v Control)
# Queries
# - Is it convenient to use other approximations for dispersion (trends/tags)?
# - Normalize by GC%?

# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("edgeR")

# load dplyr etc.
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)

# 1) Load edgeR
library(edgeR)

# 2) Load the data in R:
raw_counts <- read.csv("All_reads.csv", row.names="GenID")

# 2.5) get column names and print group names for organising
colnames <- colnames(raw_counts)

colnames

# 3) Experimental groups for basic study design: Control 8, Control 24, Control 72, Kenya 8, Kenya 24, Kenya 72, NMRI 8, NMRI 24, NMRI 72
#Groups <- factor(c("C24","C24","C24","C24","C24","C24","C24","C72","C72","C72","C72","C72","C72","C72","C72","C72","C72","C8","C8","C8","C8","C8","C8","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K24","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K72","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","K8","N24","N24","N24","N24","N24","N72","N72","N72","N72","N72","N72","N72","N72","N8","N8","N8","N8","N8","N8","N8","N8","N8"))

# 3) Experimental groups accounting for Sm transcripts study design:
# Control 8, Control 24, Control 72, Kenya 8, Kenya 8 +ve, Kenya 24, Kenya 24 +ve, Kenya 72, Kenya 72 +ve, NMRI 8, NMRI 8 +ve, NMRI 24, NMRI 72
# From list sample_ids_of_Sm_infected (all with >100 transcripts and having biological processes related to survival)

# "KK24_10" "KK24_12" "KK24_17" "KK24_7"  "KK72_12" "KK72_13" "KK72_14" "KK72_15" "KK72_16" "KK72_1"  "KK72_20"
# "KK72_3"  "KK72_4"  "KK72_5"  "KK72_6"  "KK72_8"  "KK8_11"  "KK8_2"   "KK8_4"   "KK8_6"   "KK8_8"   "KK8_9"  
# "KN8_1" 

Groups <- factor(c("C24","C24","C24","C24","C24","C24","C24","C72","C72","C72","C72","C72","C72","C72","C72","C72","C72","C8","C8","C8","C8","C8","C8","K24","K24Sm","K24","K24Sm","K24","K24","K24","K24","K24Sm","K24","K24","K24","K24","K24","K24","K24","K24Sm","K24","K24","K72Sm","K72","K72","K72Sm","K72Sm","K72Sm","K72Sm","K72Sm","K72","K72","K72Sm","K72Sm","K72Sm","K72Sm","K72Sm","K72","K72Sm","K72","K8","K8Sm","K8","K8","K8","K8","K8","K8","K8","K8Sm","K8","K8","K8Sm","K8","K8Sm","K8","K8Sm","K8Sm","N24","N24","N24","N24","N24","N72","N72","N72","N72","N72","N72","N72","N72","N8Sm","N8","N8","N8","N8","N8","N8","N8","N8"))


# 4) Production of the object to be used by edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 5) Library filtering:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 6) Normalized
edgeR_Input<-calcNormFactors(edgeR_Input)

# 6.5) Experimental matrix design
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples) 
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 7) Differential expression calculations
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 7.5) PLOT PCA
pdf('BSUD_transcriptomics_ALL_PCA.pdf',width=10,height=6.5,useDingbats=FALSE)

library(RColorBrewer)
# Define color palettes
blues <- brewer.pal(3, "Blues")
oranges <- brewer.pal(6, "Oranges")
greens <- brewer.pal(4, "Greens")
# Combine the palettes
group_colors <- c(blues, oranges, greens)
group_levels <- levels(edgeR_Input$samples$group)
# Adjust the right margin to make space for the legend
par(mar = c(5, 4, 4, 8))  # Increase the right margin
plotMDS(edgeR_Input, col = group_colors[edgeR_Input$samples$group], pch = 19)
#legend("right", legend = group_levels, col = group_colors, pch = 19)
legend("right", inset = c(-0.13, 0), legend = group_levels, col = group_colors, pch = 19, xpd = TRUE)
dev.off()


# 8) Binomial Distribution Fitting
fit <- glmQLFit(edgeR_Input,design_matrix)
#residual deviance plot
plotQLDisp(fit)

# plotMDS: Multidimensional scaling plot 
#col.cond <- c(rep("red",4),rep("blue",4),rep("green4",4)) # Colors of the conditions (The columns are taken in order)
#png("MDS_plot_8hrs.png", width = 650, height = 650)
#plotMDS(edgeR_Input, col= col.cond, cex=5, main="Control 8 Vs Kenya 8 Vs NMRI 8", pch=".")
#legend(x=1.8,y=0.4,c("C8","K8","N8"),cex=0.9,col=c("red","blue","green4"),pch=".") # Legend (I have not managed to make it look good)
#dev.off()

# tidy enviroment
rm(group_colors)
rm(group_levels)
rm(blues)
rm(greens)
rm(oranges)

#### Biological Coefficient of Variation (BCV) ####
# diagnostic tool for understanding the variability in your RNA-seq data

if (!dir.exists("BCV_analysis")) {
  dir.create("BCV_analysis")
  message("Directory 'BCV_analysis' created.")
} else {
  message("Directory 'BCV_analysis' already exists.")
}

# generate plot - As gene expression increases (right side of the plot), variability tends to decrease, and the BCV stabilizes
pdf('BCV_analysis/Bsud_BCV_plot.pdf',width=10,height=6.5,useDingbats=FALSE)
plotBCV(edgeR_Input)
dev.off()

# Step 1: Calculate the common dispersion
edgeR_Input <- estimateCommonDisp(edgeR_Input)
common_dispersion <- edgeR_Input$common.dispersion
common_BCV <- sqrt(common_dispersion)  # BCV is the square root of common dispersion
print(paste("Common BCV:", round(common_BCV, 3)))
# "Common BCV: 0.553"

# Step 2: Calculate trended and tagwise dispersions
edgeR_Input <- estimateTrendedDisp(edgeR_Input)
edgeR_Input <- estimateTagwiseDisp(edgeR_Input)

# Extract trended BCV (for comparison across expression levels)
trended_BCV <- sqrt(edgeR_Input$trended.dispersion)
# For a quick look at the range of trended BCV values:
print(paste("Trended BCV range:", round(range(trended_BCV), 3)))
# "Trended BCV range: 0.315" "Trended BCV range: 2.105"

# Extract tagwise BCV (gene-specific BCV for each gene)
tagwise_BCV <- sqrt(edgeR_Input$tagwise.dispersion)

# Optional: Get a summary of tagwise BCVs (10th, 50th, and 90th percentiles)
summary_BCV <- quantile(tagwise_BCV, probs = c(0.1, 0.5, 0.9))
print(paste("Tagwise BCV percentiles - 10%:", round(summary_BCV[1], 3),
            "50% (median):", round(summary_BCV[2], 3),
            "90%:", round(summary_BCV[3], 3)))
# "Tagwise BCV percentiles - 10%: 0.338 50% (median): 0.594 90%: 1.083"

## Investigate genes with particularly high BCV (i.e. top 1%)

# Define a threshold for "high" BCV, e.g., top 1%
high_BCV_threshold <- quantile(tagwise_BCV, 0.99)  # Top 1% of BCV values

# Identify indices of genes with high BCV
high_BCV_indices <- which(tagwise_BCV >= high_BCV_threshold)

# Extract the GeneID and BCV values for high BCV genes
high_BCV_genes <- data.frame(
  GeneID = rownames(edgeR_Input$counts)[high_BCV_indices],  # GeneIDs from the DGEList
  BCV = tagwise_BCV[high_BCV_indices]  # Corresponding BCV values
)

write_xlsx(high_BCV_genes, path = "BCV_analysis/Bsud_genes_top_1perc_BCV.xlsx")

# tidy enviroment
rm(high_BCV_indices)
rm(high_BCV_threshold)
rm(summary_BCV)
rm(tagwise_BCV)
rm(trended_BCV)
rm(common_BCV)
rm(common_dispersion)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


#### DE of Bsud genes 8 hours ####



# 9) differential expression for comparing across 8 hours for all 23598 genes:
# make folder for DE_tables
if (!dir.exists("PW_DE_tables")) {
  dir.create("PW_DE_tables")
  message("Directory 'PW_DE_tables' created.")
} else {
  message("Directory 'PW_DE_tables' already exists.")
}
# Control 8 (C8) vs Kenya 8 (K8)
C8vsK8 <- makeContrasts(C8-K8, levels=design_matrix)
C8vsK8_Test <- glmQLFTest(fit, contrast=C8vsK8) # The contrast argument in this case requests a statistical test of the null hypothesis that C8-K8 is equal to zero.
topTags(C8vsK8_Test, 23598)-> C8vsK8_2024_11_11
write.table(C8vsK8_2024_11_11, file="PW_DE_tables/Results_C8vsK8_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 8 (C8) vs Kenya 8 +ve (K8Sm)
C8vsK8Sm <- makeContrasts(C8-K8Sm, levels=design_matrix)
C8vsK8Sm_Test <- glmQLFTest(fit, contrast=C8vsK8Sm) # The contrast argument in this case requests a statistical test of the null hypothesis that C8-K8 is equal to zero.
topTags(C8vsK8Sm_Test,23598)-> C8vsK8Sm_2024_11_11
write.table(C8vsK8Sm_2024_11_11, file="PW_DE_tables/Results_C8vsK8Sm_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 8 (C8) vs NMRI 8 (N8)
C8vsN8 <- makeContrasts(C8-N8, levels=design_matrix)
C8vsN8_Test <- glmQLFTest(fit, contrast=C8vsN8)
topTags(C8vsN8_Test,23598)-> C8vsN8_2024_11_11
write.table(C8vsN8_2024_11_11, file="PW_DE_tables/Results_C8vsN8_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 8 (C8) vs NMRI 8 +ve (N8Sm)
C8vsN8Sm <- makeContrasts(C8-N8Sm, levels=design_matrix)
C8vsN8Sm_Test <- glmQLFTest(fit, contrast=C8vsN8Sm)
topTags(C8vsN8Sm_Test,23598)-> C8vsN8Sm_2024_11_11
write.table(C8vsN8Sm_2024_11_11, file="PW_DE_tables/Results_C8vsN8Sm_2024_11_11_test.tab", sep = "\t", col.names = NA)

# excluding below comparison for exposed snails
# Kenya 8 (K8) vs NMRI 8 (N8)
#K8vsN8 <- makeContrasts(K8-N8, levels=design_matrix)
#K8vsN8_Test <- glmQLFTest(fit, contrast=K8vsN8)
#topTags(K8vsN8_Test,10088)-> K8vsN8_2024_11_11
#write.table(K8vsN8_2024_11_11, file="PW_DE_tables/Results_K8vsN8_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Summary and results for less stringent FDR <5% (matching Lu et al. 2022)
# This generates slightly larger lists of up and down regulated genes for most comparisons
summary(decideTests(C8vsK8_Test, p.value=0.05))
summary(decideTests(C8vsN8_Test, p.value=0.05))
summary(decideTests(C8vsK8Sm_Test, p.value=0.05))
summary(decideTests(C8vsN8Sm_Test, p.value=0.05))
# Summary and results for stringent FDR <0.1%
summary(decideTests(C8vsK8_Test, p.value=0.001))
summary(decideTests(C8vsN8_Test, p.value=0.001))
#summary(decideTests(K8vsN8_Test, p.value=0.001))
summary(decideTests(C8vsK8Sm_Test, p.value=0.001))
summary(decideTests(C8vsN8Sm_Test, p.value=0.001))

# Plot log-fold change against log-counts per million, with DE genes highlighted:
if (!dir.exists("PW_logFC_logCPM_plots")) {
  dir.create("PW_logFC_logCPM_plots")
  message("Directory 'PW_logFC_logCPM_plots' created.")
} else {
  message("Directory 'PW_logFC_logCPM_plots' already exists.")
}
pdf('PW_logFC_logCPM_plots/C8vsK8_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C8vsK8_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

pdf('PW_logFC_logCPM_plots/C8vsN8_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C8vsN8_Test)
abline(h=c(-1, 1), col="blue") 
dev.off()

pdf('PW_logFC_logCPM_plots/C8vsK8Sm_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C8vsK8Sm_Test)
abline(h=c(-1, 1), col="blue") 
dev.off()

pdf('PW_logFC_logCPM_plots/C8vsN8Sm_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C8vsN8Sm_Test)
abline(h=c(-1, 1), col="blue") 
dev.off()

#plotMD(K8vsN8_Test)
#abline(h=c(-1, 1), col="blue")

# Merge FDR 0.1% (p<0.001) test results for those with significantly up and down regulated genes

tempmatrix1 <- topTags(C8vsK8_Test, p.value = 0.001,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsK8")
colnames(tempmatrix1) <- new_column_names

tempmatrix2 <- topTags(C8vsN8_Test, p.value = 0.001,23598)$table
tempmatrix2 <- cbind(geneID = rownames(tempmatrix2), tempmatrix2)
# Change column names
new_column_names <- colnames(tempmatrix2)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsN8")
colnames(tempmatrix2) <- new_column_names

# no significant up or down genes for C8vsK8Sm or C8vsN8Sm

#join matrices for this comparison
#temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C8N8K8_p0.001 <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')


# Merge FDR 5% (p<0.05) test results for those with significantly up and down regulated genes

tempmatrix1 <- topTags(C8vsK8_Test, p.value = 0.05,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsK8")
colnames(tempmatrix1) <- new_column_names

tempmatrix2 <- topTags(C8vsN8_Test, p.value = 0.05,23598)$table
tempmatrix2 <- cbind(geneID = rownames(tempmatrix2), tempmatrix2)
# Change column names
new_column_names <- colnames(tempmatrix2)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsN8")
colnames(tempmatrix2) <- new_column_names

tempmatrix3 <- topTags(C8vsK8Sm_Test, p.value = 0.05,23598)$table
tempmatrix3 <- cbind(geneID = rownames(tempmatrix3), tempmatrix3)
# Change column names to include which pairwise comparison results represent
new_column_names <- colnames(tempmatrix3)
new_column_names[-1] <- paste0(new_column_names[-1], "_C8vsK8Sm")
colnames(tempmatrix3) <- new_column_names

# no up or down for C8vsN8Sm

temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C8N8K8_p0.05 <- full_join(as.data.frame(temp_result_df), as.data.frame(tempmatrix2), by = 'geneID')

# tidy enviroment
rm(new_column_names)
rm(tempmatrix1)
rm(tempmatrix2)
rm(tempmatrix3)
rm(temp_result_df)


#### DE of Bsud genes 24 hours ####
# make folder for DE_tables
if (!dir.exists("PW_DE_tables")) {
  dir.create("PW_DE_tables")
  message("Directory 'PW_DE_tables' created.")
} else {
  message("Directory 'PW_DE_tables' already exists.")
}

# Control 24 (C24) vs Kenya 24 (K24)
C24vsK24 <- makeContrasts(C24-K24, levels=design_matrix)
C24vsK24_Test <- glmQLFTest(fit, contrast=C24vsK24)
topTags(C24vsK24_Test,23598)-> C24vsK24_2024_11_11
write.table(C24vsK24_2024_11_11, file="PW_DE_tables/Results_C24vsK24_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 24 (C24) vs Kenya 24 +ve (K24Sm)
C24vsK24Sm <- makeContrasts(C24-K24Sm, levels=design_matrix)
C24vsK24Sm_Test <- glmQLFTest(fit, contrast=C24vsK24Sm)
topTags(C24vsK24Sm_Test,23598)-> C24vsK24Sm_2024_11_11
write.table(C24vsK24Sm_2024_11_11, file="PW_DE_tables/Results_C24vsK24Sm_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 24 (C24) vs NMRI 24 (N24)
C24vsN24 <- makeContrasts(C24-N24, levels=design_matrix)
C24vsN24_Test <- glmQLFTest(fit, contrast=C24vsN24)
topTags(C24vsN24_Test,23598)-> C24vsN24_2024_11_11
write.table(C24vsN24_2024_11_11, file="PW_DE_tables/Results_C24vsN24_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Kenya 24 (K24) vs NMRI 24 (N24)
#K24vsN24 <- makeContrasts(K24-N24, levels=design_matrix)
#K24vsN24_Test <- glmQLFTest(fit, contrast=K24vsN24)
#topTags(K24vsN24_Test,10088)-> K24vsN24_2022_01_22
#write.table(K24vsN24_2022_01_22, file="Results_K24vsN24_2022_01_22_test.tab", sep = "\t", col.names = NA)

# Summary and results
summary(decideTests(C24vsK24_Test, p.value=0.001))
summary(decideTests(C24vsN24_Test, p.value=0.001))
summary(decideTests(C24vsK24Sm_Test, p.value=0.001))
summary(decideTests(C24vsK24_Test, p.value=0.05))
summary(decideTests(C24vsN24_Test, p.value=0.05))
summary(decideTests(C24vsK24Sm_Test, p.value=0.05))


# Plot log-fold change against log-counts per million, with DE genes highlighted:
if (!dir.exists("PW_logFC_logCPM_plots")) {
  dir.create("PW_logFC_logCPM_plots")
  message("Directory 'PW_logFC_logCPM_plots' created.")
} else {
  message("Directory 'PW_logFC_logCPM_plots' already exists.")
}


pdf('PW_logFC_logCPM_plots/C24vsK24_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C24vsK24_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

pdf('PW_logFC_logCPM_plots/C24vsK24Sm_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C24vsK24Sm_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

pdf('PW_logFC_logCPM_plots/C24vsN24_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C24vsN24_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

# Merge results FDR <0.1%

tempmatrix1 <- topTags(C24vsK24_Test, p.value = 0.001,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24")
colnames(tempmatrix1) <- new_column_names

# no signif for C24vsK24Sm at 0.1% level
#tempmatrix2 <- topTags(C24vsK24Sm_Test, p.value = 0.001,23598)$table
#tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
#new_column_names <- colnames(tempmatrix1)
#new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24Sm")
#colnames(tempmatrix1) <- new_column_names

tempmatrix3 <- topTags(C24vsN24_Test, p.value = 0.001,23598)$table
tempmatrix3 <- cbind(geneID = rownames(tempmatrix3), tempmatrix3)
new_column_names <- colnames(tempmatrix3)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsN24")
colnames(tempmatrix3) <- new_column_names

#join matrices for this comparison
#temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C24N24K24_p0.001 <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix3), by = 'geneID')

# Merge results FDR <5%

tempmatrix1 <- topTags(C24vsK24_Test, p.value = 0.05,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24")
colnames(tempmatrix1) <- new_column_names

tempmatrix2 <- topTags(C24vsK24Sm_Test, p.value = 0.05,23598)$table
tempmatrix2 <- cbind(geneID = rownames(tempmatrix2), tempmatrix2)
new_column_names <- colnames(tempmatrix2)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsK24Sm")
colnames(tempmatrix2) <- new_column_names

tempmatrix3 <- topTags(C24vsN24_Test, p.value = 0.05,23598)$table
tempmatrix3 <- cbind(geneID = rownames(tempmatrix3), tempmatrix3)
new_column_names <- colnames(tempmatrix3)
new_column_names[-1] <- paste0(new_column_names[-1], "_C24vsN24")
colnames(tempmatrix3) <- new_column_names

#join matrices for this comparison
temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C24N24K24_p0.05 <- full_join(as.data.frame(temp_result_df), as.data.frame(tempmatrix3), by = 'geneID')


# tidy enviroment
rm(new_column_names)
rm(tempmatrix1)
rm(tempmatrix2)
rm(tempmatrix3)
rm(temp_result_df)


#### DE of Bsud genes 72 hours ####
# make folder for DE_tables
if (!dir.exists("PW_DE_tables")) {
  dir.create("PW_DE_tables")
  message("Directory 'PW_DE_tables' created.")
} else {
  message("Directory 'PW_DE_tables' already exists.")
}

# Control 72 (C72) vs Kenya 72 (K72)
C72vsK72 <- makeContrasts(C72-K72, levels=design_matrix)
C72vsK72_Test <- glmQLFTest(fit, contrast=C72vsK72)
topTags(C72vsK72_Test,23598)-> C72vsK72_2024_11_11
write.table(C72vsK72_2024_11_11, file="PW_DE_tables/Results_C72vsK72_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 72 (C72) vs Kenya 72 +ve (K72Sm)
C72vsK72Sm <- makeContrasts(C72-K72Sm, levels=design_matrix)
C72vsK72Sm_Test <- glmQLFTest(fit, contrast=C72vsK72Sm)
topTags(C72vsK72Sm_Test,23598)-> C72vsK72Sm_2024_11_11
write.table(C72vsK72Sm_2024_11_11, file="PW_DE_tables/Results_C72vsK72Sm_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Control 72 (C72) vs NMRI 72 (N72)
C72vsN72 <- makeContrasts(C72-N72, levels=design_matrix)
C72vsN72_Test <- glmQLFTest(fit, contrast=C72vsN72)
topTags(C72vsN72_Test,23598)-> C72vsN72_2024_11_11
write.table(C72vsN72_2024_11_11, file="PW_DE_tables/Results_C72vsN72_2024_11_11_test.tab", sep = "\t", col.names = NA)

# Kenya 72 (K72) vs NMRI 72 (N72)
#K72vsN72 <- makeContrasts(K72-N72, levels=design_matrix)
#K72vsN72_Test <- glmQLFTest(fit, contrast=K72vsN72)
#topTags(K72vsN72_Test,10088)-> K72vsN72_2022_01_22
#write.table(K72vsN72_2022_01_22, file="Results_K72vsN72_2022_01_22_test.tab", sep = "\t", col.names = NA)

# Summary and results
summary(decideTests(C72vsK72_Test, p.value=0.001))
summary(decideTests(C72vsN72_Test, p.value=0.001))
summary(decideTests(C72vsK72Sm_Test, p.value=0.001))
summary(decideTests(C72vsK72_Test, p.value=0.05))
summary(decideTests(C72vsN72_Test, p.value=0.05))
summary(decideTests(C72vsK72Sm_Test, p.value=0.05))

# Plot log-fold change against log-counts per million, with DE genes highlighted:
if (!dir.exists("PW_logFC_logCPM_plots")) {
  dir.create("PW_logFC_logCPM_plots")
  message("Directory 'PW_logFC_logCPM_plots' created.")
} else {
  message("Directory 'PW_logFC_logCPM_plots' already exists.")
}

pdf('PW_logFC_logCPM_plots/C72vsK72_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C72vsK72_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

pdf('PW_logFC_logCPM_plots/C72vsK72Sm_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C72vsK72Sm_Test)
abline(h=c(-1, 1), col="blue")
dev.off()

pdf('PW_logFC_logCPM_plots/C72vsN72_logFC_logCPM_plot.pdf',width=5,height=5,useDingbats=FALSE)
plotMD(C72vsN72_Test)
abline(h=c(-1, 1), col="blue") 
dev.off()


# Merge results FDR <0.1%

tempmatrix1 <- topTags(C72vsK72_Test, p.value = 0.001,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72")
colnames(tempmatrix1) <- new_column_names

tempmatrix2 <- topTags(C72vsN72_Test, p.value = 0.001,23598)$table
tempmatrix2 <- cbind(geneID = rownames(tempmatrix2), tempmatrix2)
new_column_names <- colnames(tempmatrix2)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsN72")
colnames(tempmatrix2) <- new_column_names

tempmatrix3 <- topTags(C72vsK72Sm_Test, p.value = 0.001,23598)$table
tempmatrix3 <- cbind(geneID = rownames(tempmatrix3), tempmatrix3)
new_column_names <- colnames(tempmatrix3)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72Sm")
colnames(tempmatrix3) <- new_column_names

#join matrices for this comparison
temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C72N72K72_p0.001 <- full_join(as.data.frame(temp_result_df), as.data.frame(tempmatrix3), by = 'geneID')

# Merge results FDR <5%

tempmatrix1 <- topTags(C72vsK72_Test, p.value = 0.05,23598)$table
tempmatrix1 <- cbind(geneID = rownames(tempmatrix1), tempmatrix1)
new_column_names <- colnames(tempmatrix1)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72")
colnames(tempmatrix1) <- new_column_names

tempmatrix2 <- topTags(C72vsN72_Test, p.value = 0.05,23598)$table
tempmatrix2 <- cbind(geneID = rownames(tempmatrix2), tempmatrix2)
new_column_names <- colnames(tempmatrix2)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsN72")
colnames(tempmatrix2) <- new_column_names

tempmatrix3 <- topTags(C72vsK72Sm_Test, p.value = 0.05,23598)$table
tempmatrix3 <- cbind(geneID = rownames(tempmatrix3), tempmatrix3)
new_column_names <- colnames(tempmatrix3)
new_column_names[-1] <- paste0(new_column_names[-1], "_C72vsK72Sm")
colnames(tempmatrix3) <- new_column_names

#join matrices for this comparison
temp_result_df <- full_join(as.data.frame(tempmatrix1), as.data.frame(tempmatrix2), by = 'geneID')
result_C72N72K72_p0.05 <- full_join(as.data.frame(temp_result_df), as.data.frame(tempmatrix3), by = 'geneID')

# tidy enviroment
rm(new_column_names)
rm(tempmatrix1)
rm(tempmatrix2)
rm(tempmatrix3)
rm(temp_result_df)

####### MERGE ALL FDR 0.01% RESULTS DATA TOGETHER #########

temp_merge_result_df <- full_join(as.data.frame(result_C8N8K8_p0.001), as.data.frame(result_C24N24K24_p0.001), by = 'geneID')
ALL_DE_pwcomp_p0.001 <- full_join(as.data.frame(temp_merge_result_df), as.data.frame(result_C72N72K72_p0.001), by = 'geneID')

write.table(ALL_DE_pwcomp_p0.001, file= "PW_DE_tables/ALL_DE_pwcomp_p0.001.tab", sep = "\t", col.names = NA)

####### MERGE ALL FDR 5% RESULTS DATA TOGETHER #########

temp_merge_result_df <- full_join(as.data.frame(result_C8N8K8_p0.05), as.data.frame(result_C24N24K24_p0.05), by = 'geneID')
ALL_DE_pwcomp_p0.05 <- full_join(as.data.frame(temp_merge_result_df), as.data.frame(result_C72N72K72_p0.05), by = 'geneID')

write.table(ALL_DE_pwcomp_p0.05, file= "PW_DE_tables/ALL_DE_pwcomp_p0.05.tab", sep = "\t", col.names = NA)



