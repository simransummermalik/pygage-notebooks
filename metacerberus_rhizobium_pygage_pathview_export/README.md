# MetaCerberus Rhizobium Output to pyGAGE and Pathview Plus

## Project overview

This repository contains a workflow test using existing MetaCerberus rhizobium output as input for pyGAGE and Pathview Plus.

The goal was to inspect the MetaCerberus output structure, identify which files are appropriate for downstream analysis, run pyGAGE on the KEGG functional count table, and visualize the Rhizobium versus Non-Rhizobium signal on KEGG nitrogen metabolism using Pathview Plus.

This is not an RNA-seq expression analysis. The MetaCerberus count table used here is mostly binary, so the results are interpreted as differences in functional presence or functional representation across genomes.

The main comparison is:

```text
Rhizobium genomes versus Non-Rhizobium genomes
```

The main pathway focus is:

```text
KEGG nitrogen metabolism, ko00910
```

## Repository source

The workflow starts from the rhizobium test output already present in the MetaCerberus repository.

Relevant output folder:

```text
results/rhizobium/23-06-01_rhizobium/step_10-visualizeData/combined/
```

Relevant files inspected:

```text
FOAM_Loading_Matrix.tsv
KEGG_Loading_Matrix.tsv
counts_FOAM.tsv
counts_KEGG.tsv
stats.tsv
pathview/
```

The main analysis file used for pyGAGE was:

```text
counts_KEGG.tsv
```

The main reason for using `counts_KEGG.tsv` is that it contains KEGG Orthology IDs as rows and sample columns as the matrix values.

Example structure:

```text
ID       GIC31_complete   JQGI01_HUD   RW2_complete   S17   SH31   gsA_SM152B   gsB_3841   gsC_SM41   gsD_SM51   gsE_BIHB1217
K00001   0                0            0              0     0      1            1          1          1          1
K00002   0                0            0              0     0      1            0          0          0          0
K00003   0                0            0              0     0      1            1          1          1          1
```

## Output file inspection

The first file inspected was:

```text
FOAM_Loading_Matrix.tsv
```

That file had columns like:

```text
ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10
```

Because the columns are principal component columns, this table appears to be a loading or visualization summary table rather than the feature by sample matrix needed for pyGAGE.

The workflow then moved to the count tables:

```text
counts_FOAM.tsv
counts_KEGG.tsv
```

The KEGG count table was selected for the main pyGAGE and Pathview Plus workflow because it uses KEGG Orthology identifiers.

## Sample groups

The rhizobium test data defines two groups: Rhizobium and Non-Rhizobium.

The samples used in this count matrix were:

```python
rhizobium_cols = [
    "gsA_SM152B",
    "gsB_3841",
    "gsC_SM41",
    "gsD_SM51",
    "gsE_BIHB1217"
]

non_rhizobium_cols = [
    "JQGI01_HUD",
    "RW2_complete",
    "SH31",
    "GIC31_complete",
    "S17"
]
```

One Rhizobium sample listed in the repository documentation, `gsF_Vaf-108`, was not present in this specific `counts_KEGG.tsv` output, so it was not included in the analysis.

The final comparison used five Rhizobium genomes and five Non-Rhizobium genomes.

## Data interpretation

The values in `counts_KEGG.tsv` are mostly binary.

```text
0 = KEGG function not detected
1 = KEGG function detected
```

Because of this, the analysis should be described as a functional presence or absence comparison.

Correct interpretation language:

```text
more represented
more frequently present
functional presence
functional potential
KEGG function representation
```

Language to avoid:

```text
upregulated
downregulated
differential expression
gene expression increase
transcriptional change
```

This distinction is important because the input is not gene expression data. It is functional annotation output across genomes.

## Workflow summary

The analysis workflow was:

```text
MetaCerberus rhizobium output
to output file inspection
to KEGG count table selection
to Rhizobium versus Non-Rhizobium grouping
to pyGAGE input construction
to nitrogen related KO sanity check
to official KEGG nitrogen metabolism pyGAGE test
to Pathview Plus input construction
to KEGG nitrogen metabolism pathway visualization
```

## pyGAGE input construction

The pyGAGE input table was built from `counts_KEGG.tsv`.

The table was reordered so that the `ID` column came first, followed by the Non-Rhizobium reference samples and then the Rhizobium comparison samples.

```python
kegg_pygage_input = counts_kegg[
    ["ID"] + non_rhizobium_cols + rhizobium_cols
].copy()

print("pyGAGE input shape:", kegg_pygage_input.shape)
print(kegg_pygage_input.columns.tolist())

display(kegg_pygage_input.head(10))
```

The resulting table had:

```text
4340 rows
11 columns
1 KEGG Orthology ID column
10 genome/sample columns
```

## pyGAGE preparation

The table was converted to Polars before passing it into pyGAGE.

```python
import polars as pl

kegg_pygage_pl = pl.from_pandas(kegg_pygage_input)
```

Because the input values were mostly binary, the comparison was prepared using difference based preparation rather than fold based preparation.

```python
from pygage.core import GAGEPreparation, GAGEAnalysis

prep = GAGEPreparation()

prepared_rhizobium = prep.prepare_expression(
    kegg_pygage_pl,
    ref_indices=list(range(1, 1 + len(non_rhizobium_cols))),
    samp_indices=list(range(
        1 + len(non_rhizobium_cols),
        1 + len(non_rhizobium_cols) + len(rhizobium_cols)
    )),
    comparison="unpaired",
    use_fold=False
)

prepared_rhizobium = prepared_rhizobium.with_columns(
    kegg_pygage_pl["ID"]
)
```

The choice of `use_fold=False` was important because binary data contains many zeros. Fold based ratios would not be the best first approach for this type of matrix.

The prepared table was checked for missing and infinite values.

```python
prepared_pd = prepared_rhizobium.to_pandas()

numeric_cols = [col for col in prepared_pd.columns if col != "ID"]

nan_count = prepared_pd[numeric_cols].isna().sum().sum()
inf_count = np.isinf(prepared_pd[numeric_cols].to_numpy()).sum()

print("Total NaN values:", nan_count)
print("Total infinite values:", inf_count)
```

Result:

```text
Total NaN values: 0
Total infinite values: 0
```

## Nitrogen related KO sanity check

Before running the official KEGG nitrogen metabolism pathway, I first ran a small targeted sanity check using nitrogenase related KEGG Orthology sets.

This was used to confirm that the workflow produced a biologically sensible signal before scaling to a larger pathway set.

```python
nitrogen_ko_sets = {
    "core_nitrogenase_complex": [
        "K02586",
        "K02588",
        "K02591"
    ],
    "nitrogenase_cofactor_and_assembly": [
        "K02584",
        "K02585",
        "K02587",
        "K02589",
        "K02590",
        "K02592"
    ],
    "nitrogen_metabolism_context": [
        "K00360",
        "K00366",
        "K00367",
        "K00368",
        "K02567",
        "K02568",
        "K02575"
    ]
}
```

Before running pyGAGE, the notebook checked which requested KO IDs were actually present in the input matrix.

```python
available_kos = set(kegg_pygage_input["ID"].astype(str))

ko_presence_rows = []

for set_name, kos in nitrogen_ko_sets.items():
    found = [ko for ko in kos if ko in available_kos]
    missing = [ko for ko in kos if ko not in available_kos]

    ko_presence_rows.append({
        "ko_set": set_name,
        "requested": len(kos),
        "found": len(found),
        "missing": len(missing),
        "found_kos": ", ".join(found),
        "missing_kos": ", ".join(missing)
    })

ko_presence_df = pd.DataFrame(ko_presence_rows)
display(ko_presence_df)
```

## pyGAGE sanity check results

The targeted nitrogenase related KO sets produced a positive Rhizobium signal.

```python
gage = GAGEAnalysis()

results_rhizobium = gage.run_gage(
    prepared_rhizobium,
    nitrogen_ko_sets,
    gene_col="ID",
    set_size_range=(2, 100),
    same_dir=True,
    test_method="t-test"
)
```

Key output:

```text
nitrogenase_cofactor_and_assembly
set_size = 4
stat_mean = 1.000000
p_greater = 0.00000
q_greater = 0.00000

core_nitrogenase_complex
set_size = 3
stat_mean = 0.933333
p_greater = 0.004620
q_greater = 0.004620
```

Interpretation:

```text
Nitrogenase related KEGG functions are more represented in the Rhizobium genomes than in the Non-Rhizobium genomes in this test output.
```

This result is consistent with the expected biology of Rhizobium as a nitrogen fixation associated group.

## Official KEGG nitrogen metabolism pathway

After the targeted sanity check, the workflow tested the official KEGG nitrogen metabolism pathway:

```text
ko00910
```

The KO list for the pathway was pulled from KEGG.

```python
import requests

def get_kegg_pathway_kos(pathway_id):
    url = f"https://rest.kegg.jp/link/ko/{pathway_id}"
    response = requests.get(url)
    response.raise_for_status()

    kos = []
    for line in response.text.strip().split("\n"):
        if not line:
            continue
        _, ko_entry = line.split("\t")
        kos.append(ko_entry.replace("ko:", ""))

    return sorted(set(kos))

nitrogen_metabolism_kos = get_kegg_pathway_kos("ko00910")
```

Only KOs present in the MetaCerberus KEGG count table were used for the pyGAGE pathway test.

```python
available_kos = set(kegg_pygage_input["ID"].astype(str))

found_nitrogen_pathway_kos = [
    ko for ko in nitrogen_metabolism_kos
    if ko in available_kos
]

missing_nitrogen_pathway_kos = [
    ko for ko in nitrogen_metabolism_kos
    if ko not in available_kos
]
```

The pathway level pyGAGE run used:

```python
kegg_pathway_sets = {
    "KEGG_nitrogen_metabolism_ko00910": found_nitrogen_pathway_kos
}

results_kegg_nitrogen = gage.run_gage(
    prepared_rhizobium,
    kegg_pathway_sets,
    gene_col="ID",
    set_size_range=(2, 500),
    same_dir=True,
    test_method="t-test"
)
```

Result:

```text
KEGG_nitrogen_metabolism_ko00910
set_size = 24
stat_mean = 0.466667
p_greater = 0.042952
q_greater = 0.042952
```

Interpretation:

```text
The KEGG nitrogen metabolism pathway is more represented in the Rhizobium genomes compared with the Non-Rhizobium genomes in this MetaCerberus output.
```

## Pathview Plus input construction

Pathview Plus needs one value per KO.

For this workflow, the mapped value was:

```text
Rhizobium mean presence minus Non-Rhizobium mean presence
```

This value ranges from `-1` to `1`.

```text
positive value = more frequently present in Rhizobium
negative value = more frequently present in Non-Rhizobium
value near zero = similarly represented across groups
```

Input construction:

```python
pathview_input = counts_kegg[
    ["ID"] + non_rhizobium_cols + rhizobium_cols
].copy()

pathview_input["non_rhizobium_mean"] = pathview_input[
    non_rhizobium_cols
].mean(axis=1)

pathview_input["rhizobium_mean"] = pathview_input[
    rhizobium_cols
].mean(axis=1)

pathview_input["rhizobium_minus_non_rhizobium"] = (
    pathview_input["rhizobium_mean"]
    - pathview_input["non_rhizobium_mean"]
)

pathview_simple_input = pathview_input[
    ["ID", "rhizobium_minus_non_rhizobium"]
].copy()

pathview_simple_input = pathview_simple_input.rename(columns={
    "ID": "gene_id",
    "rhizobium_minus_non_rhizobium": "value"
})
```

The input was saved as:

```text
pathview_input_rhizobium_vs_non_rhizobium.csv
```
## Interpretation of the pyGAGE results

The pyGAGE results show that nitrogen related KEGG functions are more represented in the Rhizobium genomes than in the Non-Rhizobium genomes in this MetaCerberus test output.

This matters because the strongest targeted signal came from nitrogenase related KO sets. The `core_nitrogenase_complex` set includes KEGG Orthology IDs associated with the nitrogenase enzyme complex, which is central to biological nitrogen fixation. The `nitrogenase_cofactor_and_assembly` set also showed a strong Rhizobium signal, which makes biological sense because nitrogen fixation does not only require the structural nitrogenase genes. It also requires proteins involved in assembling and supporting the nitrogenase machinery.

The targeted pyGAGE sanity check produced:

```text
core_nitrogenase_complex
stat_mean = 0.933333
q_greater = 0.004620

nitrogenase_cofactor_and_assembly
stat_mean = 1.000000
q_greater = 0.000000
```

Because the input matrix is binary, a high positive `stat_mean` means these KO sets are more consistently present in the Rhizobium genomes compared with the Non-Rhizobium genomes. It does not mean the genes are more highly expressed. It means the functional capacity is more represented across that group of genomes.

The broader KEGG nitrogen metabolism pathway also showed a positive Rhizobium signal:

```text
KEGG_nitrogen_metabolism_ko00910
set_size = 24
stat_mean = 0.466667
p_greater = 0.042952
q_greater = 0.042952
```

This result is weaker than the targeted nitrogenase sanity check, which makes sense. The full KEGG nitrogen metabolism pathway contains many functions beyond nitrogen fixation itself, including nitrate reduction, nitrite reduction, denitrification, assimilation, and related nitrogen transformations. Not every part of nitrogen metabolism should be expected to separate Rhizobium from Non-Rhizobium in the same direction. A broader pathway will usually dilute the very specific nitrogenase signal because it includes multiple nitrogen related processes, not only nitrogen fixation.

So the interpretation is not simply:

```text
Rhizobium has more nitrogen metabolism.
```

A more careful interpretation is:

```text
The Rhizobium genomes show stronger representation of nitrogenase related functions, and this targeted signal is still detectable when expanded to the broader KEGG nitrogen metabolism pathway.
```

That distinction is important because the targeted nitrogenase KO sets are more directly tied to nitrogen fixation, while KEGG nitrogen metabolism is a larger pathway category that includes several nitrogen cycling processes.

## Pathview Plus run

The Pathview Plus run used:

```text
pathway_id = ko00910
species = ko
gene_idtype = KEGG
```

```python
import pathview
import polars as pl

gene_data_for_pathview = pl.from_pandas(pathview_simple_input)

pathview_plus_dir = combined_dir / "pathview_plus_outputs"
pathview_plus_dir.mkdir(exist_ok=True)

nitrogen_pathview_result = pathview.pathview(
    pathway_id="ko00910",
    gene_data=gene_data_for_pathview,
    species="ko",
    kegg_dir=pathview_plus_dir,
    kegg_native=True,
    output_format="png",
    gene_idtype="KEGG",
    out_suffix="rhizobium_vs_non_rhizobium_nitrogen_metabolism",
    map_symbol=False
)
```

Pathview Plus produced:

```text
ko00910.png
ko00910.rhizobium_vs_non_rhizobium_nitrogen_metabolism.png
```

The first file is the base nitrogen metabolism map.

The second file is the mapped Pathview Plus output using the Rhizobium versus Non-Rhizobium comparison.

## Final results

### Targeted nitrogenase KO sets

| KO set | Set size | Direction | q value |
|---|---:|---|---:|
| nitrogenase cofactor and assembly | 4 | more represented in Rhizobium | 0.00000 |
| core nitrogenase complex | 3 | more represented in Rhizobium | 0.00462 |

### Official KEGG nitrogen metabolism pathway

| Pathway | Set size | Direction | p greater | q greater |
|---|---:|---|---:|---:|
| KEGG nitrogen metabolism, ko00910 | 24 | more represented in Rhizobium | 0.042952 | 0.042952 |

## Main interpretation

The workflow successfully connected MetaCerberus rhizobium KEGG output to pyGAGE and Pathview Plus.

The strongest result is that nitrogenase related KO sets and the broader KEGG nitrogen metabolism pathway show greater functional representation in the Rhizobium genomes compared with the Non-Rhizobium genomes.

Because the input data is binary, this should be interpreted as a difference in functional presence or functional potential, not expression.

## Output files

Recommended files to keep:

```text
notebooks/MetaCerberus_Rhizobium_FOAM_Output_to_pyGAGE_Pathview_Plus.ipynb

data/counts_KEGG.tsv
data/counts_FOAM.tsv
data/stats.tsv

results/pathview_input_rhizobium_vs_non_rhizobium.csv
results/pathview_nitrogen_metabolism_rhizobium_vs_non_rhizobium.png
results/ko00910_all_blue_nitrogen_metabolism.png
```

## repository structure

```text
metacerberus-rhizobium-pygage-pathview/
│
├── README.md
│
├── notebooks/
│   └── MetaCerberus_Rhizobium_FOAM_Output_to_pyGAGE_Pathview_Plus.ipynb
│
├── data/
│   ├── counts_KEGG.tsv
│   ├── counts_FOAM.tsv
│   └── stats.tsv
│
├── results/
│   ├── pathview_input_rhizobium_vs_non_rhizobium.csv
│   ├── pathview_nitrogen_metabolism_rhizobium_vs_non_rhizobium.png
│   └── ko00910_all_blue_nitrogen_metabolism.png
│
└── notes/
    └── workflow_notes.md
```

## Troubleshooting notes

### FOAM loading matrix

`FOAM_Loading_Matrix.tsv` is not the best direct pyGAGE input because it contains principal component columns.

Use `counts_KEGG.tsv` for the KEGG based pyGAGE workflow.

### Binary values

The count matrix is mostly 0 and 1.

This means the analysis compares functional presence and absence across groups.

### Difference based preparation

`use_fold=False` was used because fold based ratios are not ideal for binary matrices with many zeros.

### Pathview Plus mapping

Pathview Plus may report that many input IDs are unmapped.

This is expected because the full KEGG count table contains thousands of KOs, while a single KEGG pathway only contains a subset of those KOs.

### Correct Pathview file

The mapped comparison image is:

```text
ko00910.rhizobium_vs_non_rhizobium_nitrogen_metabolism.png
```

The base pathway image is:

```text
ko00910.png
```

## Limitations

This notebook is a workflow test and not a complete biological study.

Key limitations:

1. The dataset is small.
2. The values are binary rather than continuous abundance values.
3. The analysis compares functional presence and absence, not expression.
4. The first KO sets were manually selected as a sanity check.
5. The official KEGG nitrogen metabolism test is stronger, but still limited by which KOs are present in the MetaCerberus output.
6. Additional pathways should be tested before making broader claims.
7. The Pathview Plus output should be compared with the existing MetaCerberus Pathview output or the R Pathview workflow.

## Next steps

Possible next steps:

1. Compare the Pathview Plus figure to the existing MetaCerberus Pathview output.
2. Run additional KEGG pathway sets through pyGAGE.
3. Test whether FOAM based outputs are better suited for some functional categories than KEGG based outputs.
4. Write a cleaner vignette showing how to move from MetaCerberus output to pyGAGE and Pathview Plus.
5. Add a small pre-analysis checker that identifies whether a MetaCerberus output table is a loading matrix, count matrix, or pathway output.
6. Document when to use `counts_KEGG.tsv` versus loading matrices.
7. Explore how this workflow could be integrated into a GUI or notebook template.

## Final interpretation

This workflow shows that MetaCerberus rhizobium KEGG output can be passed into pyGAGE and Pathview Plus in a biologically interpretable way.

The strongest signal came from nitrogenase related KO sets, especially the core nitrogenase complex and nitrogenase cofactor or assembly functions. This is the expected direction because Rhizobium is associated with nitrogen fixation.

The broader KEGG nitrogen metabolism pathway also showed greater representation in the Rhizobium genomes, although the signal was less specific than the targeted nitrogenase sets. That makes sense because KEGG nitrogen metabolism includes many nitrogen cycling functions, not only nitrogen fixation.

The main biological interpretation is:

```text
Rhizobium genomes in this MetaCerberus output show stronger functional representation of nitrogen fixation related KEGG orthology groups, and this signal is still visible at the broader KEGG nitrogen metabolism pathway level.
```

The main workflow interpretation is:

```text
counts_KEGG.tsv is a usable bridge between MetaCerberus output, pyGAGE pathway analysis, and Pathview Plus pathway visualization.
```

The most important caution is that this is presence and absence data. The result describes functional potential across genomes, not gene expression or transcriptional activity.

So the cleanest final conclusion is:

```text
This notebook validates a Python based handoff from MetaCerberus KEGG functional output into pyGAGE and Pathview Plus, using Rhizobium versus Non-Rhizobium nitrogen metabolism as a test case.
```
