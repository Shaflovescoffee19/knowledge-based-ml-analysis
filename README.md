# Knowledge-Based ML: DEGs involved in the Co-Morbidity of Patients who Undergo Maintenance Hemodialysis with Heartfailure Pipeline

A comprehensive R pipeline that integrates statistical evidence with biological knowledge to rank upregulated genes and interpret their functional impact.

## 🔬 Overview

This pipeline performs advanced analysis on upregulated genes by combining:
- Statistical significance from differential expression analysis
- Biological knowledge from Gene Ontology and pathway databases
- Disease association data
- Comprehensive visualization and reporting

## ✨ Features

- **Gene Ontology (GO) Enrichment**: Biological process analysis using clusterProfiler
- **KEGG Pathway Analysis**: Metabolic and signaling pathway enrichment
- **Disease Association**: Integration with DisGeNET, GWAS Catalog, and OMIM databases
- **Intelligent Gene Prioritization**: Combined statistical and knowledge-based scoring system
- **Rich Visualizations**: Dot plots, scatter plots, and distribution charts
- **Comprehensive Exports**: Multiple CSV outputs for further analysis

## 📋 Requirements

### System Requirements
- R >= 4.0.0
- Internet connection (for database queries)

### R Packages
Install all required packages by running:
```r
source("requirements.R")
```

**CRAN Packages:**
- `ggplot2` - Data visualization
- `dplyr` - Data manipulation  
- `enrichR` - Functional enrichment analysis

**Bioconductor Packages:**
- `clusterProfiler` - Statistical analysis and visualization of functional profiles
- `org.Hs.eg.db` - Human genome annotation database
- `enrichplot` - Visualization of functional enrichment results

## 🚀 Quick Start

### 1. Install Dependencies
```r
source("requirements.R")
```

### 2. Prepare Your Data
Your input data should be a dataframe named `Upregulated` with these columns:

| Column | Description | Example |
|--------|-------------|---------|
| `Gene` | Gene symbols (HGNC format) | "TP53", "BRCA1", "MYC" |
| `log2FoldChange` | Log2 fold change values | 2.34, -1.87, 3.12 |
| `padj` | Adjusted p-values | 0.001, 0.003, 0.0001 |

### 3. Load Your Data
```r
# From CSV file
Upregulated <- read.csv("data/your_upregulated_genes.csv")

# From Excel file
library(readxl)
Upregulated <- read_excel("data/your_upregulated_genes.xlsx")

# Verify data structure
str(Upregulated)
head(Upregulated)
```

### 4. Run the Pipeline
```r
source("src/upregulated_genes_pipeline.R")
```

## 📊 Output Files

The pipeline generates several output files:

| File | Description |
|------|-------------|
| `All_Gene_Scores.csv` | Complete gene rankings with all scores |
| `Top100_Genes.csv` | Top 100 prioritized genes |
| `GO_BP_Results.csv` | Gene Ontology biological processes |
| `KEGG_Pathways.csv` | KEGG pathway enrichment results |
| `Disease_DisGeNET.csv` | Disease associations from DisGeNET |
| `Disease_GWAS_Catalog_2019.csv` | GWAS catalog associations |
| `Disease_OMIM_Expanded.csv` | OMIM disease associations |

## 🧮 Methodology

### Gene Prioritization Algorithm

The pipeline uses a sophisticated weighted scoring system:

**Final Score = (Statistical Score × 0.6) + (Knowledge Score × 0.4)**

Where:
- **Statistical Score**: -log10(adjusted p-value)
- **Knowledge Score**: Number of GO terms + pathway associations

### Significance Categories

Genes are automatically categorized based on statistical significance:

| Category | P-value Range | Description |
|----------|---------------|-------------|
| Highly Significant | padj < 0.001 | Strong statistical evidence |
| Very Significant | 0.001 ≤ padj < 0.01 | Moderate statistical evidence |
| Significant | 0.01 ≤ padj < 0.05 | Standard significance threshold |
| Marginal | padj ≥ 0.05 | Below significance threshold |

### Enrichment Analysis

1. **GO Biological Process**: Identifies overrepresented biological processes
2. **KEGG Pathways**: Finds enriched metabolic and signaling pathways  
3. **Disease Associations**: Links genes to known diseases and phenotypes

## 📈 Visualizations

The pipeline automatically generates:

1. **GO Terms Dot Plot**: Top 10 enriched biological processes
2. **Gene Score Scatter Plot**: Statistical vs Knowledge evidence for top 50 genes
3. **Significance Distribution**: Bar chart showing gene categories

## 🗂️ Project Structure

```
upregulated-genes-pipeline/
├── README.md                    # This file
├── requirements.R               # Package installation script
├── src/
│   └── upregulated_genes_pipeline.R  # Main analysis script
├── data/
│   └── README.md               # Data format documentation
├── outputs/
│   └── .gitkeep               # Placeholder for output files
└── docs/
    └── analysis_guide.md      # Detailed analysis guide
```

## 💡 Tips for Best Results

1. **Gene Symbols**: Ensure gene symbols are current HGNC approved symbols
2. **Data Quality**: Remove rows with missing values in required columns
3. **Sample Size**: Pipeline works best with 100+ upregulated genes
4. **Internet Connection**: Required for database queries (GO, KEGG, disease databases)

## 🔧 Troubleshooting

### Common Issues

**"Package not found" errors:**
```r
# Reinstall missing packages
source("requirements.R")
```

**"No gene can be mapped" warning:**
- Check that gene symbols are in HGNC format
- Remove outdated or non-standard gene symbols

**Empty enrichment results:**
- Increase gene list size
- Check internet connection
- Verify gene symbols are correct

## 📝 Example Usage

```r
# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichR)
library(ggplot2)
library(dplyr)

# Load your data
Upregulated <- read.csv("data/my_upregulated_genes.csv")

# Check data structure
str(Upregulated)
cat("Loaded", nrow(Upregulated), "upregulated genes\n")

# Run the complete pipeline
source("src/upregulated_genes_pipeline.R")

# View top results
top_genes <- read.csv("Top100_Genes.csv")
head(top_genes)
```

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

If you encounter issues:

1. Check this README for common solutions
2. Search existing [Issues](../../issues)
3. Create a new issue with:
   - Your R version (`R.version`)
   - Package versions (`sessionInfo()`)
   - Complete error message
   - Sample of your input data structure

## 🔗 Related Resources

- [clusterProfiler Book](https://yulab-smu.top/biomedical-knowledge-mining-book/)
- [Gene Ontology Consortium](http://geneontology.org/)
- [KEGG Pathway Database](https://www.genome.jp/kegg/pathway.html)
- [DisGeNET](https://www.disgenet.org/)

## ✨ Acknowledgments

- Built with [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- Uses data from Gene Ontology, KEGG, DisGeNET, and other public databases
- Visualization powered by ggplot2

---

*Happy analyzing! 🧬📊*

# knowledge-based-ml-analysis

