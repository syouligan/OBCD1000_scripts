#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Annotate 1000 genes with osteocyte signature information
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/OBCD1000/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/OBCD1000/")
  place <- "wolfpack"
}

# Libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(colorblindr)
library(gplots)
library(forcats)

# Load data
genes1000 <- read.csv('data/1000_genes_list_data_annotated.csv', header = TRUE, row.names = 1)

# Load osteocyte gene activity data
geneActivity <- read.csv("~/cloudstor/projects/bone_comparison/project_results/osteocyte_signature/Dysplasia_analysis_2019/Gene_activity_table_collated_datasets_dysplasia_organs_Nosology2019.csv", header = TRUE, row.names = 1)

# Annotate with osteocyte expression data
idx <- match(genes1000$Ensembl, rownames(geneActivity))
genes1000$Description <- geneActivity$Description [idx]
genes1000$Biotype <- geneActivity$Biotype [idx]
genes1000$Skeletal_GO_Term <- geneActivity$Skeletal_GO_Term [idx]
genes1000$Skeletal_MP_Term <- geneActivity$Skeletal_MP_Term [idx]
genes1000$Skeletal_annotation <- geneActivity$Skeletal_annotation [idx]
genes1000$Bone_transcriptome <- geneActivity$Bone_transcriptome [idx]
genes1000$LFC <- geneActivity$LFC [idx]
genes1000$Osteocyte_signature <- geneActivity$Bone_specific_genes [idx]
genes1000$NSGD <- geneActivity$Nosology2019 [idx]

# Annotate with IMPC data
IMPC_skele <- read.table("~/cloudstor/scott_projects/osteochondro_signature/data/IMPC_20210525/skeleton_phenotype.tsv", header = TRUE, sep = '\t')
idx <- match(genes1000$GeneSymbol, IMPC_skele$Gene)
genes1000$IMPC_skeleton_significant <- genes1000$GeneSymbol %in% IMPC_skele$Gene
genes1000$IMPC_phenotype <- IMPC_skele$Phenotype [idx]

# Plot OBCD gene characteritsics
# --------------------------------------------------------------------------
biotypes <- c("protein_coding", "lincRNA", "processed_transcript", "antisense", "miRNA")

# Plot gene biotype proportions
data.frame(table(genes1000$Biotype)) %>%
  mutate(Biotype = factor(Var1, levels = biotypes)) %>%
  ggplot(aes(x = "", y = Freq, fill = Biotype)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_OkabeIto(darken = -0.4) +
  theme_void() +
  coord_polar("y", start=0) +
  ggsave("project_results/Biotype_proportions.pdf")

# Venn diagram overlap between phenotyped genes and NSGD
GOI_lists <- list("OBCD1000" = genes1000$Ensembl, "NSGD" = rownames(geneActivity[geneActivity$Nosology2019, ]))
pdf("project_results/OBCD1000_NSGD_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Venn diagram overlap between phenotyped genes and osteocyte transcriptome
GOI_lists <- list("OBCD1000" = genes1000$Ensembl, "OsteocyteTranscriptome" = rownames(geneActivity[geneActivity$Bone_transcriptome & geneActivity$Biotype %in% biotypes, ]))
pdf("project_results/OBCD1000_OcyTranscriptome_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Venn diagram overlap between phenotyped genes and osteocyte signature
GOI_lists <- list("OBCD1000" = genes1000$Ensembl, "OsteocyteSignature" = rownames(geneActivity[geneActivity$Bone_specific_genes & geneActivity$Biotype %in% biotypes, ]))
pdf("project_results/OBCD1000_OcySig_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Characterise OBCD and IMPC significant genes
# --------------------------------------------------------------------------

# Plot overlap between IMPC significant phenotypes and OBCD significant phenotypes
GOI_lists <- list("OBCDsig" = genes1000[ genes1000$Score > 0,"GeneSymbol"], "IMPCsig" = genes1000[ genes1000$IMPC_skeleton_significant, "GeneSymbol"])
pdf("project_results/OBCD1000_OBCDsig_IMPCsig_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Plot phenotype frequency for each population
lists <- attributes(ItemsList)$intersections

data.frame(table(genes1000[genes1000$GeneSymbol %in% lists$`OBCDsig:IMPCsig`, 'IMPC_phenotype'])) %>%
  mutate(IMPC_OBCD_overlap = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x = IMPC_OBCD_overlap, y = Freq)) +
  geom_bar(stat="identity", width = 0.7, fill = 'steelblue') +
  theme_minimal() +
  coord_flip() +
  ggsave("project_results/IMPC_OBCD_overlap_phenotypes_proportions.pdf")

data.frame(table(genes1000[genes1000$GeneSymbol %in% lists$`IMPCsig`, 'IMPC_phenotype'])) %>%
  mutate(IMPC_only = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x = IMPC_only, y = Freq)) +
  geom_bar(stat="identity", width = 0.7, fill = 'steelblue') +
  theme_minimal() +
  coord_flip() +
  ggsave("project_results/IMPC_only_phenotypes_proportions.pdf")

# Plot OBCD and Osteocyte transcriptome and signature genes overlap
# --------------------------------------------------------------------------

# Subset to genes with significant phenotype
genes1000_sig <- genes1000[genes1000$Score > 0, ]

# Venn diagram overlap between OBCD significant genes and osteocyte transcriptome
GOI_lists <- list("OBCDsig" = genes1000_sig$Ensembl, "OsteocyteTranscriptome" = rownames(geneActivity[geneActivity$Bone_transcriptome & geneActivity$Biotype %in% biotypes, ]))
pdf("project_results/OBCDsig_OcyTranscriptome_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Venn diagram overlap between each OBCD significant phenotype and osteocyte transcriptome
phenotypes <- c("Cortical_Structure", "Cortical_Strength", "Trabecular_Structure", "Trabecular_Strength", "Cortical_Quality", "Trabecular_Quality", "Mahalanobis")
for (i in phenotypes) {
  GOI_lists <- list("OBCDsig" = genes1000_sig[genes1000_sig[ ,i] == 1,'Ensembl'], "OsteocyteTranscriptome" = rownames(geneActivity[geneActivity$Bone_transcriptome & geneActivity$Biotype %in% biotypes, ]))
  pdf(paste0("project_results/OBCDsig_", i,"_OcyTranscriptome_overlap.pdf"))
  ItemsList <- venn(GOI_lists, show.plot = TRUE)
  dev.off()
}

# Venn diagram overlap between OBCD significant genes and osteocyte signature
GOI_lists <- list("OBCDsig" = genes1000_sig$Ensembl, "OsteocyteSignature" = rownames(geneActivity[geneActivity$Bone_specific_genes & geneActivity$Biotype %in% biotypes, ]))
pdf("project_results/OBCDsig_OcySig_overlap.pdf")
ItemsList <- venn(GOI_lists, show.plot = TRUE)
dev.off()

# Venn diagram overlap between each OBCD significant phenotype and osteocyte signature
phenotypes <- c("Cortical_Structure", "Cortical_Strength", "Trabecular_Structure", "Trabecular_Strength", "Cortical_Quality", "Trabecular_Quality", "Mahalanobis")
for (i in phenotypes) {
  GOI_lists <- list("OBCDsig" = genes1000_sig[genes1000_sig[ ,i] == 1,'Ensembl'], "OsteocyteSignature" = rownames(geneActivity[geneActivity$Bone_specific_genes & geneActivity$Biotype %in% biotypes, ]))
  pdf(paste0("project_results/OBCDsig_", i,"_OcySig_overlap.pdf"))
  ItemsList <- venn(GOI_lists, show.plot = TRUE)
  dev.off()
}

# Calculate overall hypergeometric enrichment in osteocyte transcriptome and osteocyte signature for OBCD significant genes
# --------------------------------------------------------------------------

# Define universe as genes with biotypes in the 1000 lines
universe <- geneActivity[geneActivity$Biotype %in% biotypes, ]

# Calculate overall enrichment of OBCD significant genes in the Osteocyte Signature
Overlap <- length(unique(genes1000_sig[genes1000_sig$Osteocyte_signature, 'Ensembl']))
Group1 <- length(unique(genes1000_sig$Ensembl))
Group2 <- length(unique(rownames(geneActivity[geneActivity$Bone_specific_genes & geneActivity$Biotype %in% biotypes, ])))
Universe <- length(unique(rownames(universe)))
Total_enrichment <- phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)
Total_expected <- Group1*Group2/Universe
Total_foldChange <- Overlap/Total_expected
Percent_signature <- Overlap/Group1*100

OcySigHyp <- c('OBCD_phenotype' = 'All_phenoytpes', 'Overlap' = Overlap, 'Group1' = Group1, 'Group2' = Group2, 'Universe' = Universe, 'Total_enrichment' = Total_enrichment, 'Total_expected' = Total_expected, 'Total_foldChange' = Total_foldChange, 'Percent_signature' = Percent_signature, 'Set' = 'Osteocyte_signature')

# Calculate overall enrichment of OBCD significant genes in the Osteocyte Transcriptome
Overlap <- length(unique(genes1000_sig[genes1000_sig$Bone_transcriptome, 'Ensembl']))
Group1 <- length(unique(genes1000_sig$Ensembl))
Group2 <- length(unique(rownames(geneActivity[geneActivity$Bone_transcriptome & geneActivity$Biotype %in% biotypes, ])))
Universe <- length(unique(rownames(universe)))
Total_enrichment <- phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)
Total_expected <- Group1*Group2/Universe
Total_foldChange <- Overlap/Total_expected
Percent_expressed <- Overlap/Group1*100

OcyTranscriptomeHyp <- c('OBCD_phenotype' = 'All_phenoytpes', 'Overlap' = Overlap, 'Group1' = Group1, 'Group2' = Group2, 'Universe' = Universe, 'Total_enrichment' = Total_enrichment, 'Total_expected' = Total_expected, 'Total_foldChange' = Total_foldChange, 'Percent_expressed' = Percent_expressed, 'Set' = 'Osteocyte_transcriptome')

# Calculate overall enrichment of OBCD significant genes in the Osteocyte Transcriptome
Overlap <- length(unique(genes1000_sig[genes1000_sig$NSGD, 'Ensembl']))
Group1 <- length(unique(genes1000_sig$Ensembl))
Group2 <- length(unique(rownames(geneActivity[geneActivity$Nosology2019 & geneActivity$Biotype %in% biotypes, ])))
Universe <- length(unique(rownames(universe)))
Total_enrichment <- phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)
Total_expected <- Group1*Group2/Universe
Total_foldChange <- Overlap/Total_expected
Percent_NSGD <- Overlap/Group1*100

NSGDHyp <- c('OBCD_phenotype' = 'All_phenoytpes', 'Overlap' = Overlap, 'Group1' = Group1, 'Group2' = Group2, 'Universe' = Universe, 'Total_enrichment' = Total_enrichment, 'Total_expected' = Total_expected, 'Total_foldChange' = Total_foldChange, 'Percent_NSGD' = Percent_NSGD, 'Set' = 'Nosology')


# Calculate hypergeometric enrichment of osteocyte transcriptome and osteocyte signature in specific phenotypes for OBCD significant genes
# --------------------------------------------------------------------------

# Convert to long format
rm(genes1000_sig_long)
for(i in phenotypes) {
  if(exists('genes1000_sig_long')){
    genes1000_sig_long_tmp <- genes1000_sig[genes1000_sig[ ,i] == 1, ]
    genes1000_sig_long_tmp$OBCD_phenotype <- i
    genes1000_sig_long <- rbind(genes1000_sig_long, genes1000_sig_long_tmp)
  } else {
    genes1000_sig_long <- genes1000_sig[genes1000_sig[ ,i] == 1, ]
    genes1000_sig_long$OBCD_phenotype <- i
  }
}

# Calculate hypergeometric enrichment of Osteocyte Signature genes in each group
OcySigHyperEnrich <- genes1000_sig_long %>%
  dplyr::select(OBCD_phenotype, Ensembl, Osteocyte_signature) %>%
  distinct() %>%
  group_by(OBCD_phenotype) %>%
  summarise(Overlap = sum(Osteocyte_signature), # Number of signature gene in group
            Group1 = n(), # Total number in group
            Group2 = length(unique(rownames(geneActivity[geneActivity$Bone_specific_genes & geneActivity$Biotype %in% biotypes, ]))), # Total number in signature
            Universe = length(unique(rownames(universe))) ) %>% # Total number of genes in universe
  mutate(Percent = Overlap/Group1*100,
         Pvalue = phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)) %>% # hypergeometric enrichment
  mutate(AdjustedPValue = p.adjust(Pvalue, method = 'BH'),
         Expected = Group1*Group2/Universe) %>%
  mutate(FoldChange = Overlap/Expected) %>%
  mutate(Set = 'Osteocyte_signature') 

# Calculate hypergeometric enrichment of Osteocyte Transcriptome genes in each group
OcyTranscriptomeHyperEnrich <- genes1000_sig_long %>%
  dplyr::select(OBCD_phenotype, Ensembl, Bone_transcriptome) %>%
  distinct() %>%
  group_by(OBCD_phenotype) %>%
  summarise(Overlap = sum(Bone_transcriptome), # Number of expressed gene in group
            Group1 = n(), # Total number in group
            Group2 = length(unique(rownames(geneActivity[geneActivity$Bone_transcriptome & geneActivity$Biotype %in% biotypes, ]))), # Total number expressed in osteocytes
            Universe = length(unique(rownames(universe))) ) %>% # Total number of genes in universe
  mutate(Percent = Overlap/Group1*100,
         Pvalue = phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)) %>% # hypergeometric enrichment
  mutate(AdjustedPValue = p.adjust(Pvalue, method = 'BH'),
         Expected = Group1*Group2/Universe) %>%
  mutate(FoldChange = Overlap/Expected) %>%
  mutate(Set = 'Osteocyte_transcriptome') 

# Calculate hypergeometric enrichment of NSGD genes in each group
NSGDHyperEnrich <- genes1000_sig_long %>%
  dplyr::select(OBCD_phenotype, Ensembl, NSGD) %>%
  distinct() %>%
  group_by(OBCD_phenotype) %>%
  summarise(Overlap = sum(NSGD), # Number of expressed gene in group
            Group1 = n(), # Total number in group
            Group2 = length(unique(rownames(geneActivity[geneActivity$Nosology2019 & geneActivity$Biotype %in% biotypes, ]))), # Total number expressed in osteocytes
            Universe = length(unique(rownames(universe))) ) %>% # Total number of genes in universe
  mutate(Percent = Overlap/Group1*100,
         Pvalue = phyper(Overlap-1, Group2, Universe-Group2, Group1, lower.tail= FALSE)) %>% # hypergeometric enrichment
  mutate(AdjustedPValue = p.adjust(Pvalue, method = 'BH'),
         Expected = Group1*Group2/Universe) %>%
  mutate(FoldChange = Overlap/Expected) %>%
  mutate(Set = 'NSGD') 


# Generate bubble plot of enrichment
ggplot2() +
  


