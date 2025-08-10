############################################################
## KNOWLEDGE-BASED ML: UPREGULATED GENES PIPELINE
## Author: Your Name
## Date: Sys.Date()
## Description:
##   Integrates statistical evidence with biological knowledge
##   to rank upregulated genes and interpret their functional impact.
############################################################

# -------------------------------
# STEP 0: Load required libraries
# -------------------------------
# If not installed, uncomment and run:
# install.packages(c("ggplot2", "enrichR"))
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichR)
library(ggplot2)
library(dplyr)

# -------------------------------
# STEP 1: Prepare gene list
# -------------------------------

# Assuming Upregulated is already loaded (columns: Gene, log2FoldChange, padj)
all_genes <- Upregulated$Gene
cat("📄 Loaded", length(all_genes), "upregulated genes\n")

# Use top N by significance for detailed annotation
top_genes <- head(Upregulated$Gene[order(Upregulated$padj)], 500)
cat("🥇 Using top", length(top_genes), "for detailed enrichment\n")

# -------------------------------
# STEP 2: GO Biological Process Enrichment
# -------------------------------

cat("\n[1] Performing Gene Ontology Enrichment...\n")
all_gene_ids <- bitr(all_genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

cat("Converted", nrow(all_gene_ids), "symbols to Entrez IDs\n")

go_results <- enrichGO(
  gene = all_gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 10,
  maxGSSize = 500
)

cat("GO terms found:", nrow(go_results@result), "\n")
if(nrow(go_results@result) > 0) {
  print(head(go_results@result[, c("Description","Count","pvalue")], 5))
}

# -------------------------------
# STEP 3: KEGG Pathway enrichment
# -------------------------------

cat("\n[2] Performing KEGG Pathway Enrichment...\n")
kegg_results <- enrichKEGG(
  gene = all_gene_ids$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

cat("KEGG pathways found:", nrow(kegg_results@result), "\n")
if(nrow(kegg_results@result) > 0) {
  print(head(kegg_results@result[, c("Description","Count","pvalue")], 5))
}

# -------------------------------
# STEP 4: Disease association (Enrichr)
# -------------------------------

cat("\n[3] Retrieving Disease Associations...\n")
databases <- c("DisGeNET", "GWAS_Catalog_2019", "OMIM_Expanded")
disease_results <- enrichr(all_genes, databases)

disgenet <- disease_results$DisGeNET
sig_dis <- disgenet[disgenet$Adjusted.P.value < 0.05, ]
cat("Significant diseases from DisGeNET:", nrow(sig_dis), "\n")

if(nrow(sig_dis) > 0) {
  print(head(sig_dis[, c("Term","Adjusted.P.value")], 5))
}

# -------------------------------
# STEP 5: Gene Prioritization
# -------------------------------

cat("\n[4] Scoring genes (Statistical + Knowledge)...\n")
gene_scores <- data.frame(
  Gene = all_genes,
  Stat_Score = -log10(Upregulated$padj[match(all_genes, Upregulated$Gene)]),
  Effect_Size = abs(Upregulated$log2FoldChange[match(all_genes, Upregulated$Gene)])
)

# Assign categories by p-value threshold
gene_scores$Signif_Category <- cut(
  Upregulated$padj[match(all_genes, Upregulated$Gene)],
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("Highly Significant","Very Significant","Significant","Marginal")
)

# Count GO & pathway occurrences for top genes
gene_scores$GO_Count <- 0
gene_scores$Path_Count <- 0

if(nrow(go_results@result) > 0) {
  go_map <- go_results@result$geneID
  for(i in seq_len(min(500, nrow(gene_scores)))) {
    gid <- all_gene_ids$ENTREZID[all_gene_ids$SYMBOL == gene_scores$Gene[i]]
    if(length(gid)) {
      gene_scores$GO_Count[i] <- sum(sapply(go_map, function(x) gid[1] %in% strsplit(x,"/")[[1]]))
    }
  }
}

if(nrow(kegg_results@result) > 0) {
  pw_map <- kegg_results@result$geneID
  for(i in seq_len(min(500, nrow(gene_scores)))) {
    gid <- all_gene_ids$ENTREZID[all_gene_ids$SYMBOL == gene_scores$Gene[i]]
    if(length(gid)) {
      gene_scores$Path_Count[i] <- sum(sapply(pw_map, function(x) gid[1] %in% strsplit(x,"/")[[1]]))
    }
  }
}

gene_scores$Knowledge_Score <- gene_scores$GO_Count + gene_scores$Path_Count

# Final weighted score
gene_scores$Final_Score <- scale(gene_scores$Stat_Score)[,1] * 0.6 +
                           scale(gene_scores$Knowledge_Score)[,1] * 0.4

gene_scores <- arrange(gene_scores, desc(Final_Score))

cat("Top 10 genes:\n")
print(head(gene_scores[, c("Gene","Stat_Score","Knowledge_Score","Final_Score","Signif_Category")], 10))

# -------------------------------
# STEP 6: Visualizations
# -------------------------------

cat("\n[5] Generating plots...\n")

if(nrow(go_results@result) > 0) {
  dotplot(go_results, showCategory=10) + ggtitle("Top Biological Processes")
}

ggplot(gene_scores[1:50,], aes(Stat_Score, Knowledge_Score,
                               size = Final_Score, color = Signif_Category,
                               label = Gene)) +
  geom_point(alpha=0.6) +
  geom_text(check_overlap = FALSE, size=2, vjust=-0.8) +
  theme_minimal() +
  labs(title="Top 50 Genes: Statistical vs Knowledge Evidence",
       x="-log10(padj)", y="Knowledge Score")

ggplot(gene_scores, aes(Signif_Category, fill=Signif_Category)) +
  geom_bar() +
  theme_minimal() +
  labs(title="Gene Significance Distribution",
       x="Category", y="Count") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# -------------------------------
# STEP 7: Save outputs
# -------------------------------

cat("\n[6] Saving results to CSV...\n")
write.csv(gene_scores, "All_Gene_Scores.csv", row.names=FALSE)
write.csv(gene_scores[1:100, ], "Top100_Genes.csv", row.names=FALSE)

if(nrow(go_results@result) > 0) write.csv(go_results@result, "GO_BP_Results.csv", row.names=FALSE)
if(nrow(kegg_results@result) > 0) write.csv(kegg_results@result, "KEGG_Pathways.csv", row.names=FALSE)

for(db in names(disease_results)) {
  sig_db <- subset(disease_results[[db]], Adjusted.P.value < 0.05)
  if(nrow(sig_db) > 0) {
    write.csv(sig_db, paste0("Disease_", gsub(" ","_",db), ".csv"), row.names=FALSE)
  }
}

cat("\n✅ DONE: Analysis complete.\n")
cat("Main output files:\n  All_Gene_Scores.csv, Top100_Genes.csv, GO_BP_Results.csv, KEGG_Pathways.csv, Disease_* .csv\n")
