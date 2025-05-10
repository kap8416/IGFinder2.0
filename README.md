# IGFinder2.0

**IGFinder2.0** is a Python-based tool to identify, analyze, and visualize **intronless genes (IGs)** from any Ensembl-supported genome using transcript-level annotations. It integrates gene structure, chromosomal analysis, statistical comparisons, and graphical outputs — with support for future biological enrichment and expression integration.

---

## 🔍 What does IGFinder do?

IGFinder automates the discovery and analysis of **intronless genes** (genes with only one exon) and compares them with **multi-exonic genes** (genes with two or more exons). It queries Ensembl's REST API, downloads gene models, filters and classifies them, and produces:

- A clean list of classified genes
- Statistical comparisons between intronless and multi-exonic gene groups
- Plots showing gene length and distribution by chromosome

---

## ✅ Key Features (Implemented)

### 1. 🔧  Modularization

- Structured with clean Python functions and modules.
- Uses `argparse` for command-line options and `logging` for progress tracking.
- Validates all required input files.

### 2. 🧬 Robust Detection of Intronless Genes

- Classifies genes as **intronless** if they contain only **one exon**.
- Excludes genes that appear intronless but contain **UTR-based introns** (requires a UTR database).
- Automatically excludes **mitochondrial chromosome (MT)** and non-standard contigs.

### 3. 💾 Reproducible Output

- Final results are exported in **tab-separated value (.tsv)** format for reuse and compatibility.

### 4. 📈 Integrated Gene Structure Analysis

- Compares **gene lengths** of intronless vs. multi-exonic genes.
- Maps genes across standard chromosomes (1–22).

### 5. 📊 Automatic Visualizations

- **Boxplots** and **violin plots** of gene size.
- **Density plots** of gene length distributions.
- **Chromosome barplots** showing gene counts.
- Plots are saved as `.png` images.

### 6. 📉 Built-in Statistical Comparisons

- Performs a **t-test** on gene lengths.
- Performs a **Chi-squared test** on chromosome distributions.
- Results are exported as `.txt` summaries.

---

## 🧪 How to Use

```bash
python3 IGFinder.py \
  --species homo_sapiens \
  --utr_db UTR_db.txt \
  --output genes_filtrados.tsv \
  --stats \
  --plots \
  --log IGFinder_log.txt
```

- `--species`: Ensembl species name (e.g., `homo_sapiens`, `mus_musculus`)
- `--utr_db`: TSV file with known UTR-based introns to exclude false IGs
- `--output`: Output table of filtered genes
- `--stats`: Enables statistical tests
- `--plots`: Enables plot generation
- `--log`: Path to output log file (optional)

---

## 📂 Output Files

- `genes_filtrados.tsv`: Filtered gene list with type, coordinates, and biotype.
- `IGFinder_stats.txt`: Summary of statistical results.
- `IGFinder_plots_*.png`: Visualizations of length and distribution.

---

| Feature                        | IGFinder1 ❌        | IGFinder2.0 ✅               |
|-------------------------------|---------------------|-----------------------------   |
| Modular design                | ❌ Monolithic script | ✅ `fetch`, `main`, `utils`  |
| Logging                       | ❌ Print only        | ✅ Logging with timestamps   |
| Visual outputs                | ❌ None              | ✅ PNG plots, statistics     |
| Reproducibility               | ⚠️ Manual            | ✅ Fully parameterized CLI   |
| Command-line support          | ❌ No                | ✅ With `argparse`           |
| Scalability                   | ❌ Limited           | ✅ Multi-species ready       |


## 🚧 Upcoming Features

### 7. 🧬 Gene Chromosome clustering-Enrichment Analysis (Sliding Window)

- Scan chromosomes with mobile windows (e.g., 1Mb).
- **Hypergeometric test** for detecting IG enrichment.
- Mark **enriched regions visually**:
  - Red dots = genes of interest
  - Purple bands = enriched windows



### 9. 🧬 Custom Background Support

- Load a user-defined file of **multi-exonic genes** as background.
- Enable tailored statistical comparisons.

### 10.🧬 Expression and GO enrichment Integration

- Connect to **Expression Atlas** or **GTEx** via APIs.
- Functional enrichment analysis


---

## 📬 Contact

Developed and maintained by **Katia Aviña Padilla**  
Test: Lesly Cárdenas López

📧 katia.avinap@cinvestav.mx  🔗 [ResearchGate](https://www.researchgate.net/profile/Katia-Avina-Padilla)
📧 lesly.cardenas.fcqb2022@uas.edu.mx


Feel free to contribute, request features, or collaborate.

