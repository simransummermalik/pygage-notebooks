# pyGAGE + Pathview Plus Real Dataset Walkthrough  
## Pathway-Level Shifts Across Colorectal Cancer Progression Using GSE50760

This project is a real-dataset walkthrough for connecting **RNA-seq expression data**, **pyGAGE pathway analysis**, and **Pathview Plus pathway visualization**.

The goal is not only to run a package successfully, but to document the full workflow that a new user would need to understand:

1. how to retrieve public RNA-seq data from GEO,  
2. how to inspect sample metadata before analysis,  
3. how to build a clean expression matrix from many individual FPKM files,  
4. how to run pyGAGE on a real biological comparison,  
5. how to troubleshoot unstable pathway statistics, and  
6. how to map gene-level changes onto KEGG pathways using Pathview Plus.

This walkthrough uses **GSE50760**, a colorectal cancer RNA-seq dataset with three biologically meaningful sample groups:

| Group | Sample count | Description |
|---|---:|---|
| Normal colon | 18 | Normal-looking surrounding colonic epithelium |
| Primary colorectal cancer | 18 | Primary colorectal cancer tissue |
| Liver metastasis | 18 | Metastatic colorectal cancer to the liver |

The first completed analysis in this walkthrough focuses on:

```text
normal colon vs primary colorectal cancer
```

The broader project direction is to compare pathway behavior across colorectal cancer progression:

```text
normal colon → primary colorectal cancer → liver metastasis
```



---

## Dataset deep dive: What is GSE50760?

**GSE50760** is a public human colorectal cancer RNA-seq dataset designed around disease progression and tumor heterogeneity. The study generated RNA-seq data from 54 samples collected from 18 colorectal cancer patients, with each patient contributing tissue from three biologically meaningful states:

```text
normal colon
primary colorectal cancer
liver metastasis
```

This structure makes the dataset especially useful for a pathway-analysis walkthrough because it is not just a simple “healthy vs disease” comparison. It allows the analysis to be framed as a progression problem:

```text
normal-looking colon tissue
    ↓
primary colorectal cancer
    ↓
metastatic colorectal cancer in liver
```

The GEO metadata confirms that the dataset contains three groups of 18 samples:

| Tissue code | Biological group | Count | Metadata wording |
|---:|---|---:|---|
| `.1` | Primary colorectal cancer | 18 | primary colorectal cancer |
| `.2` | Normal colon | 18 | normal-looking surrounding colonic epithelium |
| `.3` | Liver metastasis | 18 | metastatic colorectal cancer to the liver |

This design is powerful for a tutorial because each comparison can answer a different biological question:

| Comparison | Biological question |
|---|---|
| Normal colon vs primary colorectal cancer | What pathway programs shift during tumor formation? |
| Primary colorectal cancer vs liver metastasis | What pathway programs shift further during metastatic progression? |
| Normal colon vs liver metastasis | What is the full pathway-level difference between normal tissue and metastatic disease? |

For this first version, the notebook focuses on the first comparison:

```text
normal colon vs primary colorectal cancer
```

That choice is intentional. It provides a clean starting point before moving into the more complex progression comparisons.

---

## Why this dataset is a strong test case for pyGAGE + Pathview Plus

This dataset is useful for more than just testing whether the code runs.

It tests whether a real user can move through the entire analysis chain:

```text
GEO metadata
    ↓
sample group labeling
    ↓
expression matrix construction
    ↓
paired comparison setup
    ↓
gene set analysis
    ↓
pathway-level interpretation
    ↓
gene-level pathway visualization
```

That matters because real public RNA-seq datasets usually do not arrive in a perfect “ready for pathway analysis” format. In this project, several practical issues appeared:

1. The sample-level supplementary file column was empty, but the files existed at the series level.
2. The data were distributed as 54 separate compressed FPKM files instead of one expression matrix.
3. Repeated outer merges were memory-unstable.
4. Fold-style preparation created `NaN` and `inf` values because FPKM data contain zeros and near-zero values.
5. Pathview Plus rendered the pathway from gene symbols, but the gene-level expression values mapped correctly only after conversion to Entrez IDs.

These issues are exactly why this workflow is useful as documentation. It does not only show the successful path. It records the real points where a student, biologist, or first time user could get stuck.


---

## Why?

RNA-seq experiments often produce long gene level tables, but biological interpretation usually happens at the pathway level. A list of changed genes is useful, but researchers usually want to know:

- Which biological programs are shifting?
- Are changes related to proliferation, invasion, immune signaling, hypoxia, or apoptosis?
- Do the pathway-level signals make sense with the tissue comparison?
- Can the pathway-level result be connected back to gene-level changes visually?

This project uses **pyGAGE** to summarize gene level expression changes into pathway level signals, then uses **Pathview Plus** to map gene level changes back onto KEGG pathway diagrams.

The result is a workflow that connects:

```text
public RNA-seq data
    ↓
clean expression matrix
    ↓
pyGAGE pathway analysis
    ↓
gene-level change table
    ↓
Entrez ID mapping
    ↓
Pathview Plus KEGG visualization
```

---

## Research question

The first question in this walkthrough is:

> Can pyGAGE detect biologically interpretable pathway-level shifts between normal colon and primary colorectal cancer using a real public RNA-seq dataset?

The longer-term question is:

> Can pyGAGE and Pathview Plus be used together to distinguish pathway programs involved in tumor formation from those involved in metastatic progression?

This is important because tumor formation and metastasis are not the same biological event. A pathway that changes when normal tissue becomes a primary tumor may not be the same pathway that changes when a primary tumor becomes metastatic.

---

## Repository structure

Recommended folder structure:

```text
gse50760_pygage/
├── data/
│   ├── raw/
│   │   ├── GSE50760_RAW.tar
│   │   └── GSE50760_RAW/
│   │       ├── GSM1228184_AMC_2.1_FPKM.txt.gz
│   │       ├── GSM1228185_AMC_3.1_FPKM.txt.gz
│   │       └── ...
│   └── processed/
│       └── gse50760_fpkm_expression_matrix_clean.csv
├── results/
│   ├── normal_vs_primary_greater_results.csv
│   ├── normal_vs_primary_less_results.csv
│   ├── normal_vs_primary_stats_results.csv
│   ├── normal_vs_primary_gene_level_change_for_pathview.csv
│   ├── normal_vs_primary_pathview_simple_input.csv
│   └── pathview_outputs/
│       ├── hsa04110.normal_vs_primary_cell_cycle.png
│       └── hsa04110.normal_vs_primary_cell_cycle_entrez.png
├── notebooks/
│   └── 02_gse50760_real_dataset_pygage.ipynb
└── README.md
```

---

## Installation

This workflow was tested in Google Colab, but the same structure can be moved into VSCode or a local Python environment.

```python
!pip install -q GEOparse pandas numpy polars matplotlib seaborn scipy pygage pathview-plus mygene
```

Core packages used:

```python
from pathlib import Path
import urllib.request
import tarfile

import GEOparse
import pandas as pd
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from pygage.core import GAGEPreparation, GAGEAnalysis
```

Pathview Plus installs from PyPI as:

```text
pathview-plus
```

but imports as:

```python
import pathview
```

This distinction matters because:

```python
import pathview_plus
```

will fail even when the package is installed correctly.

---

## Step 1: Set up project folders

```python
base_dir = Path("gse50760_pygage")
raw_dir = base_dir / "data" / "raw"
processed_dir = base_dir / "data" / "processed"
results_dir = base_dir / "results"

for folder in [raw_dir, processed_dir, results_dir]:
    folder.mkdir(parents=True, exist_ok=True)

print("Project folders are ready")
print("Raw data folder:", raw_dir)
print("Processed data folder:", processed_dir)
print("Results folder:", results_dir)
```

This keeps raw files, processed matrices, and analysis outputs separate.

---

## Step 2: Load GEO metadata

```python
gse = GEOparse.get_GEO(
    geo="GSE50760",
    destdir=str(raw_dir),
    annotate_gpl=False
)

print("GEO series loaded")
print("Number of samples:", len(gse.gsms))
```

Expected result:

```text
Number of samples: 54
```

---

## Step 3: Build a sample metadata table

GEO stores each sample as a separate GSM object. This step turns those nested sample records into a readable sample metadata table.

```python
sample_rows = []

for gsm_id, gsm in gse.gsms.items():
    title = gsm.metadata.get("title", [""])[0]
    source_name = gsm.metadata.get("source_name_ch1", [""])[0]
    characteristics = gsm.metadata.get("characteristics_ch1", [])
    supplementary_files = gsm.metadata.get("supplementary_file", [])

    sample_rows.append({
        "gsm_id": gsm_id,
        "title": title,
        "source_name": source_name,
        "characteristics": " | ".join(characteristics),
        "supplementary_files": " | ".join(supplementary_files)
    })

metadata = pd.DataFrame(sample_rows)

pd.set_option("display.max_colwidth", 250)
display(metadata.head(10))
print("Metadata shape:", metadata.shape)
```

Expected shape:

```text
(54, 5)
```

This step matters because the sample labels define every downstream comparison.

---

## Step 4: Locate the processed expression files

The individual GSM records did not expose sample-level supplementary files. The processed files were instead stored at the GEO series level.

```python
print("Series metadata keys:")
print(gse.metadata.keys())

print("\nSeries supplementary files:")
print(gse.metadata.get("supplementary_file", "No series-level supplementary_file field found"))
```

The series-level supplementary file is:

```text
GSE50760_RAW.tar
```

Download it:

```python
series_files = gse.metadata.get("supplementary_file", [])

tar_url = series_files[0]
tar_path = raw_dir / "GSE50760_RAW.tar"

urllib.request.urlretrieve(tar_url, tar_path)

print("Downloaded file:")
print(tar_path)
print("File size in MB:", round(tar_path.stat().st_size / 1_000_000, 2))
```

Expected size:

```text
7.74 MB
```

---

## Step 5: Inspect the GEO archive

```python
with tarfile.open(tar_path, "r") as tar:
    members = tar.getmembers()

print("Number of files inside archive:", len(members))

for member in members[:30]:
    print(member.name, "|", round(member.size / 1_000_000, 3), "MB")
```

Expected result:

```text
Number of files inside archive: 54
```

Each file is a compressed FPKM text file:

```text
GSM1228184_AMC_2.1_FPKM.txt.gz
GSM1228185_AMC_3.1_FPKM.txt.gz
...
```

---

## Step 6: Extract and inspect one FPKM file

```python
extract_dir = raw_dir / "GSE50760_RAW"
extract_dir.mkdir(exist_ok=True)

with tarfile.open(tar_path, "r") as tar:
    tar.extractall(path=extract_dir)

extracted_files = sorted(list(extract_dir.glob("*.txt.gz")))

print("Number of extracted files:", len(extracted_files))
```

Inspect one file:

```python
first_file = extracted_files[0]

test_expr = pd.read_csv(first_file, sep="\t", compression="gzip")

print("Shape:", test_expr.shape)
display(test_expr.head())
print("Columns:", test_expr.columns.tolist())
```

Expected structure:

```text
Shape: (23505, 2)
Columns: ['genes', 'AMC_2.1_FPKM']
```

Each file contains:

```text
genes | one_sample_FPKM
```

---

## Step 7: Build a clean expression matrix

A first attempt using repeated outer merges caused Colab to crash, even with high RAM. The issue was not the dataset size. The issue was the merge strategy.

Repeated outer merges on gene names can expand unexpectedly if duplicated gene identifiers create many-to-many matches.

The safer strategy is:

1. check that gene order is consistent across files,
2. collect each FPKM column,
3. horizontally stack columns, and
4. collapse duplicated gene names afterward.

```python
sample_frames = []
gene_order = None
order_mismatch_files = []

for file in extracted_files:
    df = pl.read_csv(file, separator="\t")

    sample_name = file.name.replace("_FPKM.txt.gz", "")
    genes = df["genes"].cast(pl.Utf8).to_list()

    if gene_order is None:
        gene_order = genes
    else:
        if genes != gene_order:
            order_mismatch_files.append(file.name)

    fpkm_col = df.columns[1]
    sample_frame = df.select(
        pl.col(fpkm_col).cast(pl.Float64).alias(sample_name)
    )
    sample_frames.append(sample_frame)

if len(order_mismatch_files) > 0:
    print("Warning: gene order mismatches found")
    print(order_mismatch_files[:10])
else:
    print("Gene order is consistent across all files")

expression_matrix = pl.DataFrame({"gene_id": gene_order})

for sample_frame in sample_frames:
    expression_matrix = expression_matrix.hstack(sample_frame)

print("Combined expression matrix shape:", expression_matrix.shape)
display(expression_matrix.head())
```

Expected shape before duplicate cleanup:

```text
(23505, 55)
```

This means:

```text
23505 genes
1 gene_id column
54 sample columns
```

Collapse duplicated gene identifiers:

```python
duplicate_count = expression_matrix["gene_id"].is_duplicated().sum()
print("Duplicated gene IDs before cleanup:", duplicate_count)

expression_matrix_clean = (
    expression_matrix
    .group_by("gene_id")
    .agg(pl.exclude("gene_id").mean())
)

print("Clean expression matrix shape:", expression_matrix_clean.shape)
print("Duplicated gene IDs after cleanup:", expression_matrix_clean["gene_id"].is_duplicated().sum())
```

Expected clean shape:

```text
(23503, 55)
```

Save the cleaned matrix:

```python
processed_dir.mkdir(parents=True, exist_ok=True)

expression_matrix_path = processed_dir / "gse50760_fpkm_expression_matrix_clean.csv"
expression_matrix_clean.write_csv(expression_matrix_path)

print("Saved clean expression matrix to:")
print(expression_matrix_path)
```

---

## Step 8: Confirm sample groups

The 54 sample columns split evenly into three tissue codes:

```python
sample_columns = [col for col in expression_matrix_clean.columns if col != "gene_id"]

sample_column_rows = []

for col in sample_columns:
    gsm_id = col.split("_")[0]
    amc_part = col.replace(gsm_id + "_", "")
    tissue_code = amc_part.split(".")[-1]

    sample_column_rows.append({
        "sample_column": col,
        "gsm_id": gsm_id,
        "amc_sample": amc_part,
        "tissue_code": tissue_code
    })

sample_column_metadata = pd.DataFrame(sample_column_rows)

display(sample_column_metadata.head(25))
print(sample_column_metadata["tissue_code"].value_counts())
```

Expected:

```text
1    18
2    18
3    18
```

Merge with GEO metadata:

```python
sample_column_metadata = sample_column_metadata.merge(
    metadata[["gsm_id", "title", "source_name", "characteristics"]],
    on="gsm_id",
    how="left"
)

display(sample_column_metadata.head(54))

print("Source name counts:")
print(sample_column_metadata["source_name"].value_counts())
```

Confirmed groups:

```text
primary colorectal cancer    18
normal colon                 18
metastasized cancer          18
```

The tissue codes mean:

```text
.1 = primary colorectal cancer
.2 = normal colon
.3 = liver metastasis
```

---

## Step 9: Create matched normal vs primary comparison

Because the dataset has matched sample IDs, normal and primary samples are paired by AMC ID.

```python
def assign_biological_group(row):
    source = str(row["source_name"]).lower()
    characteristics = str(row["characteristics"]).lower()
    title = str(row["title"]).lower()
    text = source + " " + characteristics + " " + title

    if "normal colon" in text or "normal-looking" in text:
        return "normal_colon"
    elif "metastasized cancer" in text or "metastatic colorectal cancer to the liver" in text:
        return "liver_metastasis"
    elif "primary colorectal cancer" in text:
        return "primary_crc"
    else:
        return "unknown"

sample_column_metadata["group"] = sample_column_metadata.apply(assign_biological_group, axis=1)
```

Create matched pairs:

```python
sample_column_metadata["patient_id"] = (
    sample_column_metadata["amc_sample"]
    .str.replace(r"\.\d$", "", regex=True)
)

normal_meta = (
    sample_column_metadata[sample_column_metadata["group"] == "normal_colon"]
    .sort_values("patient_id")
)

primary_meta = (
    sample_column_metadata[sample_column_metadata["group"] == "primary_crc"]
    .sort_values("patient_id")
)

paired_metadata = normal_meta[["patient_id", "sample_column"]].merge(
    primary_meta[["patient_id", "sample_column"]],
    on="patient_id",
    suffixes=("_normal", "_primary")
)

display(paired_metadata)
print("Number of matched normal-primary pairs:", paired_metadata.shape[0])
```

Expected:

```text
18 matched normal-primary pairs
```

Build the expression table:

```python
normal_cols = paired_metadata["sample_column_normal"].tolist()
primary_cols = paired_metadata["sample_column_primary"].tolist()

normal_vs_primary_cols = ["gene_id"] + normal_cols + primary_cols

expr_normal_vs_primary = expression_matrix_clean.select(normal_vs_primary_cols)

print("Normal columns:", len(normal_cols))
print("Primary columns:", len(primary_cols))
print("Expression table shape:", expr_normal_vs_primary.shape)
```

Expected:

```text
Normal columns: 18
Primary columns: 18
Expression table shape: (23503, 37)
```

---

## Step 10: Define small CRC-relevant gene sets

For the first real pyGAGE run, this project uses small biologically motivated gene groups.

These are not meant to replace KEGG or GO. They serve as a controlled real-data checkpoint before scaling to larger pathway databases.

```python
crc_gene_sets = {
    "cell_cycle_proliferation": [
        "MKI67", "CDK1", "CDK2", "CCNB1", "CCND1", "PCNA", "TOP2A"
    ],
    "wnt_beta_catenin_signaling": [
        "APC", "CTNNB1", "MYC", "AXIN2", "TCF7", "LEF1", "LGR5"
    ],
    "ecm_invasion_metastasis": [
        "MMP2", "MMP9", "MMP14", "VIM", "SNAI1", "TWIST1", "ITGA5", "ITGB1"
    ],
    "apoptosis_balance": [
        "TP53", "CASP3", "BAX", "BCL2", "FAS", "BID"
    ],
    "angiogenesis_hypoxia": [
        "VEGFA", "HIF1A", "ANGPT2", "KDR", "FLT1"
    ],
    "immune_inflammation": [
        "IL6", "TNF", "CXCL8", "CCL2", "STAT3", "NFKB1", "RELA"
    ]
}
```

Check gene presence:

```python
available_genes = set(expression_matrix_clean["gene_id"].to_list())

gene_set_presence = []

for gene_set_name, genes in crc_gene_sets.items():
    found = [gene for gene in genes if gene in available_genes]
    missing = [gene for gene in genes if gene not in available_genes]

    gene_set_presence.append({
        "gene_set": gene_set_name,
        "genes_requested": len(genes),
        "genes_found": len(found),
        "genes_missing": len(missing),
        "found_genes": ", ".join(found),
        "missing_genes": ", ".join(missing)
    })

gene_set_presence_df = pd.DataFrame(gene_set_presence)
display(gene_set_presence_df)
```

---

## Step 11: Run pyGAGE

The first attempt used fold-style preparation and produced unstable values because FPKM data contains zeros and near-zero values. This led to `NaN` and `inf` values in the prepared table.

The corrected version uses:

```python
use_fold=False
```

This prepares expression differences rather than fold-style ratios.

```python
prep = GAGEPreparation()

prepared_normal_vs_primary_diff = prep.prepare_expression(
    expr_normal_vs_primary,
    ref_indices=list(range(1, 1 + len(normal_cols))),
    samp_indices=list(range(1 + len(normal_cols), 1 + len(normal_cols) + len(primary_cols))),
    comparison="paired",
    use_fold=False
)

prepared_normal_vs_primary_diff = prepared_normal_vs_primary_diff.with_columns(
    expr_normal_vs_primary["gene_id"]
)
```

Check the prepared table:

```python
prepared_diff_pd = prepared_normal_vs_primary_diff.to_pandas()

numeric_cols = [col for col in prepared_diff_pd.columns if col != "gene_id"]

nan_count = prepared_diff_pd[numeric_cols].isna().sum().sum()
inf_count = np.isinf(prepared_diff_pd[numeric_cols].to_numpy()).sum()

print("Total NaN values:", nan_count)
print("Total infinite values:", inf_count)
```

Expected:

```text
Total NaN values: 0
Total infinite values: 0
```

Run pyGAGE:

```python
gage = GAGEAnalysis()

results_normal_vs_primary_diff = gage.run_gage(
    prepared_normal_vs_primary_diff,
    crc_gene_sets,
    gene_col="gene_id",
    set_size_range=(2, 100),
    same_dir=True,
    test_method="t-test"
)

print("pyGAGE difference based run complete")
print(results_normal_vs_primary_diff.keys())

display(results_normal_vs_primary_diff["greater"].to_pandas())
display(results_normal_vs_primary_diff["less"].to_pandas())
display(results_normal_vs_primary_diff["stats"].to_pandas())
```

---

## pyGAGE result summary!

The difference-based pyGAGE run produced stable pathway-level results.

For the normal colon vs primary colorectal cancer comparison, the strongest positive pathway-level signals were:

| Gene set | Direction | q value |
|---|---|---:|
| cell_cycle_proliferation | higher in primary CRC | 0.032560 |
| ecm_invasion_metastasis | higher in primary CRC | 0.043520 |
| immune_inflammation | higher in primary CRC | 0.049742 |
| wnt_beta_catenin_signaling | positive, near threshold | 0.058651 |

Interpretation:

- **Cell cycle/proliferation** being higher in primary tumor fits tumor biology.
- **ECM/invasion/metastasis** may reflect extracellular matrix remodeling, invasion-related behavior, or tumor microenvironment changes.
- **Immune/inflammation** may reflect inflammatory signaling or differences in immune/stromal composition in bulk tissue.
- **Wnt beta catenin signaling** is biologically relevant to colorectal cancer, but the adjusted q value is slightly above 0.05 in this small custom gene set, so it should be described as suggestive rather than strongly significant.

Important caution:

This is a walkthrough-level analysis using small custom gene sets. It should be validated with broader KEGG or GO pathway collections before making stronger biological claims.

---

## Step 12: Save pyGAGE result tables

```python
greater_normal_primary = results_normal_vs_primary_diff["greater"].to_pandas()
less_normal_primary = results_normal_vs_primary_diff["less"].to_pandas()
stats_normal_primary = results_normal_vs_primary_diff["stats"].to_pandas()

greater_normal_primary["comparison"] = "normal_colon_vs_primary_crc"
less_normal_primary["comparison"] = "normal_colon_vs_primary_crc"
stats_normal_primary["comparison"] = "normal_colon_vs_primary_crc"

results_dir.mkdir(parents=True, exist_ok=True)

greater_path = results_dir / "normal_vs_primary_greater_results.csv"
less_path = results_dir / "normal_vs_primary_less_results.csv"
stats_path = results_dir / "normal_vs_primary_stats_results.csv"

greater_normal_primary.to_csv(greater_path, index=False)
less_normal_primary.to_csv(less_path, index=False)
stats_normal_primary.to_csv(stats_path, index=False)
```

---

## Step 13: Prepare Pathview Plus input

pyGAGE gives pathway-level results. Pathview Plus needs gene-level values.

For the first Pathview Plus input, this project calculates:

```text
primary colorectal cancer mean − normal colon mean
```

```python
gene_level_change = expression_matrix_clean.select(
    ["gene_id"] + normal_cols + primary_cols
).with_columns([
    pl.mean_horizontal(normal_cols).alias("normal_mean"),
    pl.mean_horizontal(primary_cols).alias("primary_mean")
]).with_columns([
    (pl.col("primary_mean") - pl.col("normal_mean")).alias("primary_minus_normal")
]).select([
    "gene_id",
    "normal_mean",
    "primary_mean",
    "primary_minus_normal"
])

display(gene_level_change.head())
print("Gene level change table shape:", gene_level_change.shape)
```

Save full and simple Pathview inputs:

```python
full_pathview_input_path = results_dir / "normal_vs_primary_gene_level_change_for_pathview.csv"
simple_pathview_input_path = results_dir / "normal_vs_primary_pathview_simple_input.csv"

gene_level_change.write_csv(full_pathview_input_path)

pathview_simple_input = gene_level_change.select([
    "gene_id",
    pl.col("primary_minus_normal").alias("value")
])

pathview_simple_input.write_csv(simple_pathview_input_path)
```

---

## Step 14: Install and inspect Pathview Plus

```python
!pip install -q pathview-plus
```

The package installs as:

```text
pathview-plus
```

but imports as:

```python
import pathview
```

Inspect the main function:

```python
import pathview
import inspect

inspect.signature(pathview.pathview)
```

The key arguments are:

```python
pathview.pathview(
    pathway_id: str,
    gene_data: Optional[pl.DataFrame] = None,
    species: str = "hsa",
    kegg_dir: str | Path = ".",
    kegg_native: bool = True,
    output_format: str = "png",
    gene_idtype: str = "ENTREZ",
    ...
)
```

The default gene identifier type is `ENTREZ`.

This matters because the original expression matrix uses gene symbols, but KEGG mapping works more reliably with Entrez IDs.

---

## Step 15: First Pathview Plus attempt using gene symbols

```python
gene_data_for_pathview = gene_level_change.select([
    pl.col("gene_id"),
    pl.col("primary_minus_normal").alias("value")
])

pathview_dir = results_dir / "pathview_outputs"
pathview_dir.mkdir(parents=True, exist_ok=True)

cell_cycle_result = pathview.pathview(
    pathway_id="hsa04110",
    gene_data=gene_data_for_pathview,
    species="hsa",
    kegg_dir=pathview_dir,
    kegg_native=True,
    output_format="png",
    gene_idtype="SYMBOL",
    out_suffix="normal_vs_primary_cell_cycle",
    map_symbol=True
)
```

This generated a KEGG pathway image, but the returned object showed:

```text
plot_data_gene: None
```

This suggested that the image rendered, but the gene values did not map properly onto the KEGG nodes.

---

## Step 16: Convert gene symbols to Entrez IDs

```python
!pip install -q mygene
```

```python
import mygene

mg = mygene.MyGeneInfo()

pathview_input_pd = pathview_simple_input.to_pandas()
pathview_input_pd = pathview_input_pd.dropna(subset=["gene_id", "value"])

gene_symbols = pathview_input_pd["gene_id"].astype(str).unique().tolist()

query_results = mg.querymany(
    gene_symbols,
    scopes="symbol",
    fields="entrezgene,symbol",
    species="human",
    as_dataframe=False
)

mapping_rows = []

for result in query_results:
    if "notfound" in result and result["notfound"]:
        continue

    if "entrezgene" in result:
        mapping_rows.append({
            "gene_id": result["query"],
            "symbol": result.get("symbol", result["query"]),
            "entrez_id": str(result["entrezgene"])
        })

gene_mapping = pd.DataFrame(mapping_rows)

display(gene_mapping.head())
print("Mapped genes:", gene_mapping.shape[0])
```

Merge Entrez IDs with expression-change values:

```python
pathview_entrez_pd = pathview_input_pd.merge(
    gene_mapping[["gene_id", "entrez_id"]],
    on="gene_id",
    how="inner"
)

pathview_entrez_pd = pathview_entrez_pd[["entrez_id", "value"]].dropna()

pathview_entrez_pd = (
    pathview_entrez_pd
    .groupby("entrez_id", as_index=False)
    .agg({"value": "mean"})
)

gene_data_entrez_for_pathview = pl.from_pandas(pathview_entrez_pd)

display(gene_data_entrez_for_pathview.head())
print(gene_data_entrez_for_pathview.shape)
```

---

## Step 17: Successful Pathview Plus visualization

The Entrez-based version successfully mapped gene values onto the KEGG cell cycle pathway.

```python
cell_cycle_result_entrez = pathview.pathview(
    pathway_id="hsa04110",
    gene_data=gene_data_entrez_for_pathview,
    species="hsa",
    kegg_dir=pathview_dir,
    kegg_native=True,
    output_format="png",
    gene_idtype="ENTREZ",
    out_suffix="normal_vs_primary_cell_cycle_entrez",
    map_symbol=True
)

print("Pathview Plus Entrez run complete")
print(cell_cycle_result_entrez)
```

Display the figure:

```python
from IPython.display import Image, display

cell_cycle_entrez_png = pathview_dir / "hsa04110.normal_vs_primary_cell_cycle_entrez.png"

print("Does Entrez output image exist?", cell_cycle_entrez_png.exists())
print("Output image path:", cell_cycle_entrez_png)

display(Image(filename=str(cell_cycle_entrez_png)))
```

This produced an expression-colored KEGG cell cycle pathway map.

---

## Key technical lessons

### 1. GEO metadata can be split across sample and series levels

The sample-level supplementary file column was empty, but the dataset had a series-level supplementary archive.

This is a useful documentation point because users may think data is missing when it is actually stored at the series level.

### 2. Repeated outer merges can crash on small-looking datasets

The archive was only around 7.74 MB, but repeated outer merges on gene names caused memory problems.

The safer workflow was to:

1. check gene order,
2. horizontally stack FPKM columns, and
3. collapse duplicated gene identifiers afterward.

### 3. FPKM fold values can create unstable pyGAGE preparation

Using fold-style preparation created `NaN` and `inf` values because FPKM data is zero-heavy.

Using:

```python
use_fold=False
```

made the prepared expression table stable.

### 4. pyGAGE and Pathview Plus may need different gene ID handling

pyGAGE worked with gene symbols for the custom gene sets.

Pathview Plus rendered the KEGG pathway with gene symbols, but the expression values did not map correctly until the input was converted to Entrez IDs.

For KEGG-based Pathview Plus visualization, Entrez IDs were more reliable.


---

## Results-driven biological observations

This analysis is a first-pass workflow result, not a final biological discovery. Still, the pathway-level output is biologically coherent and gives useful hypotheses for the next stage.

For the **normal colon vs primary colorectal cancer** comparison, the strongest positive pyGAGE signals were:

| Gene set | Direction | Adjusted q value | Interpretation |
|---|---|---:|---|
| cell_cycle_proliferation | Higher in primary CRC | 0.032560 | Suggests a strong tumor proliferation signal |
| ecm_invasion_metastasis | Higher in primary CRC | 0.043520 | Suggests extracellular matrix remodeling or invasion-related shifts |
| immune_inflammation | Higher in primary CRC | 0.049742 | Suggests inflammatory or microenvironment-associated differences |
| wnt_beta_catenin_signaling | Positive, near threshold | 0.058651 | Biologically relevant to CRC, but should be described as suggestive |

### Observation 1: Cell cycle/proliferation was the clearest primary tumor signal

The strongest significant result was the cell cycle/proliferation gene set. This is biologically sensible because primary tumor tissue should show stronger proliferative behavior than normal-looking surrounding colon tissue.

This result also made the **KEGG cell cycle pathway** a good first Pathview Plus candidate. The workflow successfully mapped primary-minus-normal gene-level values onto the KEGG cell cycle pathway after converting gene symbols to Entrez IDs.

### Observation 2: ECM/invasion signal appeared even in the primary tumor comparison

The ECM/invasion/metastasis gene set was also significant in primary colorectal cancer compared with normal colon.

This is interesting because the comparison is not primary tumor vs metastasis yet. Seeing an ECM/invasion-related signal already in primary tumor may suggest that extracellular matrix remodeling and invasion-associated programs are not only late metastatic features. They may already be visible in the primary tumor microenvironment.

That is not a final claim yet, but it gives a strong direction for the next comparison:

```text
primary colorectal cancer vs liver metastasis
```

If the ECM/invasion signal becomes even stronger in liver metastasis, that would support the idea of a progression-amplified pathway program.

### Observation 3: Immune/inflammation was significant, but needs careful interpretation

The immune/inflammation gene set was significant in the primary tumor comparison.

This is biologically plausible, but also needs caution. Bulk RNA-seq measures expression from mixed tissue. An immune/inflammatory signal could reflect:

- inflammatory signaling inside tumor cells,
- increased immune cell infiltration,
- stromal or microenvironmental composition,
- or a combination of these.

This is a good example of why pathway results should not be interpreted as pure pathway activation inside one cell type without additional evidence.

### Observation 4: Wnt signaling was positive but near threshold

Wnt/beta-catenin signaling is central to colorectal cancer biology, and the custom Wnt group showed a positive signal. However, the adjusted q value was slightly above 0.05.

The careful interpretation is:

```text
Wnt/beta-catenin signaling showed a suggestive positive shift, but this small custom gene set alone is not enough to claim strong enrichment.
```

This should be followed up with broader KEGG/GO pathway sets and Pathview Plus visualization of KEGG Wnt signaling.

---

## Candidate discoveries and hypotheses for future testing

The current workflow produces several testable hypotheses.

### Hypothesis 1: Cell cycle is an early tumor-formation program

Because cell cycle/proliferation was significant in primary colorectal cancer compared with normal colon, it may represent an early tumor-associated pathway program rather than a metastasis-specific signal.

Future test:

```text
Compare cell cycle signal across:
normal vs primary
primary vs liver metastasis
normal vs liver metastasis
```

If cell cycle is high in primary tumor and remains high in metastasis, it may represent a persistent cancer program.

### Hypothesis 2: ECM/invasion may be progression-amplified

The ECM/invasion/metastasis signal appeared in the primary tumor comparison. The next question is whether it increases further in liver metastasis.

Future test:

```text
If ECM/invasion is stronger in primary vs metastasis than in normal vs primary,
then it may represent a metastasis-associated or progression-amplified program.
```

### Hypothesis 3: Immune/inflammatory shifts may reflect tumor microenvironment changes

The immune/inflammation signal may be telling us about changes in the tissue environment, not just cancer-cell-intrinsic signaling.

Future test:

```text
Compare immune/inflammatory pathway scores with cell-type marker gene sets.
```

Potential follow-up gene sets:

- T cell markers
- macrophage markers
- neutrophil/inflammatory markers
- fibroblast or stromal markers

### Hypothesis 4: Pathway visualization requires identifier-aware preprocessing

The Pathview Plus result showed that symbol-based input could render the pathway, but Entrez-based input was needed for successful gene-value mapping.

This is a software/documentation discovery rather than a biology discovery:

```text
For KEGG-based Pathview Plus visualization, Entrez ID input is more reliable than gene symbols.
```

This is highly useful for a ReadTheDocs or vignette workflow because it prevents users from thinking their pathway values are mapped when they are only rendering the base pathway.

---

## Future directions

### 1. Complete all three biological comparisons

The next major step is to run the same pyGAGE workflow on:

```text
primary colorectal cancer vs liver metastasis
normal colon vs liver metastasis
```

This would allow pathway signals to be classified into progression categories:

| Category | Pattern |
|---|---|
| Early tumor-formation program | significant in normal vs primary |
| Metastasis-associated program | significant in primary vs liver metastasis |
| Persistent cancer program | significant in normal vs primary and normal vs metastasis |
| Progression-amplified program | stronger in metastasis than primary |
| Rewired program | changes direction between primary and metastasis |

### 2. Replace small custom gene sets with full pathway databases

The current custom gene sets are useful for debugging and interpretation. The next version should test larger pathway collections:

```text
KEGG pathways
Gene Ontology biological process sets
Reactome pathways
MSigDB Hallmark gene sets
```

This would make the analysis less dependent on hand-selected gene groups.

### 3. Add a pathway transition table

A future version should create a summary table like:

| Pathway | Normal vs Primary | Primary vs Metastasis | Normal vs Metastasis | Classification |
|---|---:|---:|---:|---|
| Cell cycle | positive | TBD | TBD | early or persistent tumor program |
| ECM/invasion | positive | TBD | TBD | possible progression program |
| Immune/inflammation | positive | TBD | TBD | microenvironment-associated |
| Wnt signaling | suggestive | TBD | TBD | CRC-relevant, needs validation |

This would turn the notebook from a single-comparison walkthrough into a progression-aware pathway analysis.

### 4. Generate a small Pathview Plus figure panel

The first successful Pathview Plus figure was KEGG cell cycle.

Future figure candidates:

| KEGG ID | Pathway |
|---|---|
| hsa04110 | Cell cycle |
| hsa04512 | ECM-receptor interaction |
| hsa04510 | Focal adhesion |
| hsa04310 | Wnt signaling pathway |
| hsa05210 | Colorectal cancer |

A polished result could include:

```text
Figure 1: pyGAGE workflow diagram
Figure 2: normal vs primary pathway result table
Figure 3: KEGG cell cycle Pathview Plus map
Figure 4: cross-stage pathway transition table
```

### 5. Add AWS or cloud reproducibility angle

Dr. White also mentioned that an AWS approach could be useful. A future version could include:

- running the notebook on AWS SageMaker Studio Lab or SageMaker notebooks,
- storing processed expression matrices in S3,
- packaging the workflow as a reproducible cloud-based tutorial,
- using an environment file for repeatable installation,
- and adding a lightweight command-line version for larger datasets.

### 6. Build ReadTheDocs-style pages from this workflow

Possible documentation pages:

```text
01_real_dataset_overview.md
02_download_geo_data.md
03_build_expression_matrix.md
04_run_pygage_paired_comparison.md
05_troubleshoot_nan_inf.md
06_prepare_pathview_input.md
07_convert_symbols_to_entrez.md
08_visualize_kegg_with_pathview_plus.md
```

The strongest part of this project is that it records the actual user experience, including the failure points and fixes.


## Current status

Completed:

- Loaded GSE50760 metadata from GEO
- Downloaded and inspected the series-level supplementary archive
- Extracted 54 FPKM files
- Built a clean expression matrix
- Confirmed sample groups
- Built matched normal vs primary comparison
- Ran pyGAGE successfully using difference-based preparation
- Identified positive pathway-level signals in primary CRC
- Created gene-level Pathview Plus input
- Installed and inspected Pathview Plus
- Rendered KEGG cell cycle pathway
- Fixed gene mapping using Entrez IDs
- Successfully generated an expression-colored KEGG cell cycle map

---

## Next steps

### 1. Run the second pyGAGE comparison

```text
primary colorectal cancer vs liver metastasis
```

This asks what changes further when colorectal cancer metastasizes to the liver.

### 2. Run the third comparison

```text
normal colon vs liver metastasis
```

This captures the full shift from normal tissue to metastatic disease.

### 3. Compare pathway behavior across disease stages

Create a table like:

| Pathway | Normal vs Primary | Primary vs Metastasis | Normal vs Metastasis | Interpretation |
|---|---:|---:|---:|---|
| Cell cycle | high | ? | ? | tumor formation or persistent tumor program |
| ECM/invasion | high | ? | ? | invasion or microenvironment shift |
| Wnt signaling | near threshold | ? | ? | CRC relevant, needs broader validation |

### 4. Add broader KEGG or GO gene sets

The custom gene sets are useful for a first walkthrough, but the next version should use larger pathway collections.

### 5. Generate more Pathview Plus maps

Candidate pathways:

```text
hsa04110  Cell cycle
hsa04512  ECM receptor interaction
hsa04510  Focal adhesion
hsa04310  Wnt signaling pathway
hsa05210  Colorectal cancer
```

### 6. Turn this into ReadTheDocs material

Potential documentation pages:

```text
Getting started with pyGAGE on a real RNA-seq dataset
Understanding GEO sample metadata
Building an expression matrix from FPKM files
Running pyGAGE with paired samples
Troubleshooting NaN and inf values
Preparing gene-level input for Pathview Plus
Mapping gene symbols to Entrez IDs
Visualizing KEGG pathways with Pathview Plus
```

---

## Interpretation caution


The current gene sets are intentionally small and hand curated. They are useful for teaching and debugging the pyGAGE to Pathview workflow, but stronger biological conclusions would require:

- broader pathway databases,
- multiple testing across complete pathway collections,
- careful normalization choices,
- validation of gene identifiers,
- consideration of bulk tissue composition,
- and possibly comparison with established differential expression pipelines.

The strength of this project is that it documents the real analysis path clearly, including the parts that did not work on the first attempt.