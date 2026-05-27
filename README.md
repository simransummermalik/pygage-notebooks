# File Structure for Navigation

Use this README as a quick guide to the project layout and current work.

## Most Recent Work

[MetaCerberus Rhizobium pyGAGE + Pathview Plus export](metacerberus_rhizobium_pygage_pathview_export/README.md)

## Clickable Navigation

### MetaCerberus Rhizobium Workflow

- [Workflow README](metacerberus_rhizobium_pygage_pathview_export/README.md)
- [Analysis notebook](metacerberus_rhizobium_pygage_pathview_export/MetaCerberus_Rhizobium_FOAM_Output_to_pyGAGE_Pathview_Plus.ipynb)
- [KEGG counts input](metacerberus_rhizobium_pygage_pathview_export/data/counts_KEGG.tsv)
- [FOAM counts input](metacerberus_rhizobium_pygage_pathview_export/data/counts_FOAM.tsv)
- [Stats input](metacerberus_rhizobium_pygage_pathview_export/data/stats.tsv)
- [Pathview input CSV](metacerberus_rhizobium_pygage_pathview_export/results/pathview_input_rhizobium_vs_non_rhizobium.csv)
- [Nitrogen metabolism Pathview image](metacerberus_rhizobium_pygage_pathview_export/results/pathview_nitrogen_metabolism_rhizobium_vs_non_rhizobium.png)
- [Original nitrogen metabolism output](metacerberus_rhizobium_pygage_pathview_export/results/nitrogenMetabolismOutputOne.png)

### GSE50760 Real Dataset Workflow

- [GSE50760 walkthrough README](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/README.md)
- [Requirements](pygage_real_dataset_walkthrough/requirements.txt)
- [Analysis notebook](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/02_gse50760_real_dataset_pygage.ipynb)
- [Processed expression matrix](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/data/processed/gse50760_fpkm_expression_matrix_clean.csv)
- [pyGAGE stats results](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/results/normal_vs_primary_stats_results.csv)
- [Successful Pathview outputs list](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/results/successful_pathview_outputs.csv)
- [Cell cycle Pathview image](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/results/pathview_outputs/hsa04110.normal_vs_primary_cell_cycle_entrez.png)
- [ECM receptor Pathview image](pygage_real_dataset_walkthrough/gse50760_pygage_work_so_far/results/pathview_outputs/hsa04512.normal_vs_primary_ecm_receptor_entrez.png)

### Toy Dataset Attempt

- [Toy pyGAGE notebook](toy-dataset-initial-attempt/tryingoutpygage.ipynb)

## Repository Structure

```text
pygage-walkthrough/
|-- README.md
|-- metacerberus_rhizobium_pygage_pathview_export/
|   |-- README.md
|   |-- MetaCerberus_Rhizobium_FOAM_Output_to_pyGAGE_Pathview_Plus.ipynb
|   |-- data/
|   |   |-- counts_KEGG.tsv
|   |   |-- counts_FOAM.tsv
|   |   `-- stats.tsv
|   `-- results/
|       |-- pathview_input_rhizobium_vs_non_rhizobium.csv
|       |-- pathview_nitrogen_metabolism_rhizobium_vs_non_rhizobium.png
|       |-- pathview_nitrogen_metabolism.png
|       `-- nitrogenMetabolismOutputOne.png
|-- pygage_real_dataset_walkthrough/
|   |-- requirements.txt
|   `-- gse50760_pygage_work_so_far/
|       |-- README.md
|       |-- 02_gse50760_real_dataset_pygage.ipynb
|       |-- data/
|       |   |-- raw/
|       |   |   |-- GSE50760_RAW.tar
|       |   |   |-- GSE50760_family.soft.gz
|       |   |   `-- GSE50760_RAW/
|       |   |       `-- GSM*_FPKM.txt.gz
|       |   `-- processed/
|       |       `-- gse50760_fpkm_expression_matrix_clean.csv
|       `-- results/
|           |-- normal_vs_primary_greater_results.csv
|           |-- normal_vs_primary_less_results.csv
|           |-- normal_vs_primary_stats_results.csv
|           |-- normal_vs_primary_gene_level_change_for_pathview.csv
|           |-- normal_vs_primary_pathview_simple_input.csv
|           |-- successful_pathview_outputs.csv
|           `-- pathview_outputs/
|               |-- hsa04110.normal_vs_primary_cell_cycle_entrez.png
|               `-- hsa04512.normal_vs_primary_ecm_receptor_entrez.png
`-- toy-dataset-initial-attempt/
    `-- tryingoutpygage.ipynb
```
