# File Structure for Navigation

Use this README as a quick guide to the project layout and current work.

## Most Recent Work

[PyGAGE organism-specific egSymb mapping test](changes-test/README.md)

This is the current test folder for organism-specific Entrez ID and gene-symbol mapping. It includes mouse (`mmu`) outputs, a `pmav` NCBI fallback test, a notebook, and copies of the changed PyGAGE files for review.

## Clickable Navigation

### PyGAGE Organism Mapping Test

- [Test README](changes-test/README.md)
- [Test notebook](changes-test/test_organism_specific_egsymb_mapping.ipynb)
- [Mouse rich mapping output](changes-test/outputs/mmu_gene_mapping.tsv)
- [Mouse egSymb-compatible output](changes-test/outputs/mmu_egSymb.tsv)
- [pmav rich mapping test output](changes-test/outputs/pmav_gene_mapping.test25.tsv)
- [pmav egSymb-compatible test output](changes-test/outputs/pmav_egSymb.test25.tsv)
- [Copied changed-file notes](changes-test/files%20changed/README.md)
- [Copied gene ID utility file](changes-test/files%20changed/lib/gene_id_utils.py)
- [Copied command-line script](changes-test/files%20changed/bin/pygage-build_egsymb_mapping.py)
- [Copied pyproject file](changes-test/files%20changed/pyproject.toml)
- [Copied PyGAGE README](changes-test/files%20changed/pygage_README.md)

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
|-- .gitignore
|-- changes-test/
|   |-- README.md
|   |-- test_organism_specific_egsymb_mapping.ipynb
|   |-- outputs/
|   |   |-- mmu_gene_mapping.tsv
|   |   |-- mmu_egSymb.tsv
|   |   |-- pmav_gene_mapping.test25.tsv
|   |   `-- pmav_egSymb.test25.tsv
|   `-- files changed/
|       |-- README.md
|       |-- pyproject.toml
|       |-- pygage_README.md
|       |-- bin/
|       |   `-- pygage-build_egsymb_mapping.py
|       `-- lib/
|           `-- gene_id_utils.py
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

The local `pygage-dev/` clone is intentionally ignored by Git and is not part of this repository navigation.
