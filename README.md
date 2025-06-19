# 🔬 DESeq2 Shiny App

An interactive web application for performing differential gene expression analysis using **DESeq2**, built with **R Shiny**.

---

## 🚀 Features

- Upload **count matrix** and **design matrix**
- Run DESeq2 directly from the interface
- Preview and download **normalized counts**
- Quality control plots:
  - PCA
  - Elbow plot
  - Sample clustering heatmap
- Differential expression analysis:
  - Choose comparison conditions
  - Adjust thresholds for adjusted p-value and log2 fold change
  - Download DEG results
  - Visualizations: MA plot, p-value distribution, volcano plot

---

## 📦 Requirements

- R (≥ 4.0)

---

## ⚙️ Getting Started

1. Clone the repository:

```bash
git clone https://github.com/your-username/deseq2-shiny-app.git
cd deseq2-shiny-app
```
2. Launch the app in R or RStudio:

```
shiny::runApp()
```
## App Structure

```
deseq2-shiny-app/
├── global.R       # Loads and installs packages
├── ui.R           # Shiny UI components
├── server.R       # Shiny server logic
└── README.md      # Project documentation
```

## Input File Format

* Count Matrix: tab-delimited .tsv or .txt, genes in rows, samples in columns. Row names = gene IDs. Column names = sample IDs.
* Design Matrix: tab-delimited .tsv or .txt, must include at least a condition column with sample grouping.

Example:

### Count matrix

```
          sample1 sample2 sample3
GeneA        100     120     130
GeneB         50      60      65
```
### Design matrix

```
sample    condition
sample1   control
sample2   treated
sample3   treated
```
