#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Annotate 1000 genes with ensembl ids and human gene ids (hg19 and hg38)
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/OBCD1000/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/OBCD1000/")
  place <- "wolfpack"
}

library("biomaRt")
library("openxlsx")
library("org.Hs.eg.db")
library('edgeR')

# Load and prepare data
# --------------------------------------------------------------------------

# First pass labelling with Ensembl IDs
# GeneMap <- read.csv("data/GeneSymbols_1000_genes.csv", row.names = NULL, header = TRUE)
# geneActivity <- read.csv('data/Chondrocyte_signature_limma.csv', row.names = 1, header = TRUE)
# 
# idx <- match(GeneMap$GeneSymbol, geneActivity$GeneSymbol)
# GeneMap$Ensembl <- rownames(geneActivity) [idx]

# Save for curation
# write.csv(GeneMap, "data/GeneSymbols_1000_genes_ensembl_not_curated.csv", row.names = FALSE)

# Load gene list with curated mouse ensembl ids
genes1000_annot <- read.csv('data/GeneSymbols_1000_genes_ensembl_curated.csv', header = TRUE)

# Annotate with human symbols
# --------------------------------------------------------------------------

# Find human ensembl ids for protein coding genes
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2014.archive.ensembl.org")
human_orths <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = genes1000_annot$Ensembl,
                     mart = ensembl)

# Annotate mouse ENSEMBL ids with human ENSEMBL Ids
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl", host = "feb2021.archive.ensembl.org")
human_orths_new <- getBM(attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name', 'hsapiens_homolog_orthology_confidence'),
                     filters = 'ensembl_gene_id',
                     values = genes1000_annot$Ensembl,
                     mart = ensembl)

# Add human gene symbols to geneActivity dataframe
idx <- match(genes1000_annot$Ensembl, human_orths$ensembl_gene_id)
genes1000_annot$Human_Ensembl_ID_hg19 <- human_orths$hsapiens_homolog_ensembl_gene [idx]

idx <- match(genes1000_annot$Ensembl, human_orths_new$ensembl_gene_id)
genes1000_annot$Human_Ensembl_ID_hg38 <- human_orths_new$hsapiens_homolog_ensembl_gene [idx]
genes1000_annot$Human_geneSymbol_hg38 <- human_orths_new$hsapiens_homolog_associated_gene_name [idx]

# Write annotation data
write.csv(genes1000_annot, "data/GeneSymbols_1000_genes_ensembl_curated_human.csv", row.names = FALSE)

# Annotate data with gene information
genes1000 <- read.csv('data/1000_genes_list_data.csv', header = TRUE)

all_data <- merge(genes1000_annot, genes1000, by.y = 'X', by.x = 'GeneSymbol', all.x = TRUE)
rownames(all_data) <- all_data$Ensembl

write.csv(all_data, "data/1000_genes_list_data_annotated.csv", row.names = TRUE)


