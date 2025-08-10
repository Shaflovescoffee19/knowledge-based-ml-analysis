# Data Directory

## Input Data Requirements

Place your input data files in this directory. The pipeline expects a dataframe named `Upregulated` with the following structure:

### Required Columns:
- **Gene**: Gene symbols in HGNC format (e.g., "TP53", "BRCA1")
- **log2FoldChange**: Numeric values representing log2 fold changes
- **padj**: Adjusted p-values from differential expression analysis

### Example Data Format:Gene,log2FoldChange,padj
RNR2	135422.7179	9.634268263	0.375721748	25.64202984	5.19E-145	1.82E-141
RNR1	37429.05225	8.56269675	0.451740028	18.95492144	4.02E-80	1.99E-77
IGHA2	552.8423066	7.587750835	0.555462375	13.66024267	1.75E-42	1.39E-40
ALB	17.17384277	7.50870288	1.621213554	4.631532264	3.63E-06	1.53E-05
RPL23P8	42.04454494	7.257241308	0.559150343	12.97905188	1.61E-38	9.77E-37
RNA5S9	13.5458839	7.18991537	2.478679904	2.900703459	0.00372326	0.010050886
AMER2	12.29445076	7.048927954	2.432285577	2.898067571	0.003754697	0.010127964

