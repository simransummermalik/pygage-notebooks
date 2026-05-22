# File Structure for Navigation

Use this README as a quick guide to the project layout.

The most recent work is here:

[GSE50760 pyGAGE + Pathview Plus walkthrough](gse50760_pygage_work_so_far/README.md)

```text
pygage_real_dataset_walkthrough/
|-- README.md
|-- requirements.txt
`-- gse50760_pygage_work_so_far/
    |-- README.md
    |-- README_pygage_pathview_gse50760_EXPANDED.md
    |-- 02_gse50760_real_dataset_pygage.ipynb
    |-- data/
    |   |-- raw/
    |   |   |-- GSE50760_RAW.tar
    |   |   |-- GSE50760_family.soft.gz
    |   |   `-- GSE50760_RAW/
    |   |       `-- GSM*_FPKM.txt.gz
    |   `-- processed/
    |       `-- gse50760_fpkm_expression_matrix_clean.csv
    `-- results/
        |-- normal_vs_primary_greater_results.csv
        |-- normal_vs_primary_less_results.csv
        |-- normal_vs_primary_stats_results.csv
        |-- normal_vs_primary_gene_level_change_for_pathview.csv
        |-- normal_vs_primary_pathview_simple_input.csv
        |-- successful_pathview_outputs.csv
        `-- pathview_outputs/
            |-- hsa04110.normal_vs_primary_cell_cycle.png
            |-- hsa04110.normal_vs_primary_cell_cycle_entrez.png
            |-- hsa04512.normal_vs_primary_ecm_receptor_entrez.png
            |-- hsa04110.xml
            `-- hsa04512.xml
```
