# Files Changed for Organism-Specific Mapping Test

This folder contains copies of the local PyGAGE files changed for the organism-specific `egSymb` mapping test.

These are copied here so the changes can be reviewed without committing the whole `pygage-dev/` clone.

## Copied Files

```text
files changed/
|-- README.md
|-- pyproject.toml
|-- pygage_README.md
|-- bin/
|   `-- pygage-build_egsymb_mapping.py
`-- lib/
    `-- gene_id_utils.py
```

## 1. `lib/gene_id_utils.py`

Original purpose:

```text
PyGAGE already used this file for Entrez ID <-> gene symbol conversion.
```

Main change:

```python
class KEGGOrganismGeneMapper:
```

Why it was added:

```text
GeneIDConverter can convert IDs only after a mapping file already exists.
KEGGOrganismGeneMapper creates organism-specific mapping files first.
```

Important methods added:

```python
list_organisms()
validate_organism_code()
fetch_kegg_gene_list()
fetch_entrez_mapping()
build_mapping()
build_ncbi_gene_mapping()
write_mapping_files()
```

What each part does:

- `list_organisms()` pulls KEGG's organism-code table.
- `validate_organism_code()` checks codes such as `mmu`, `hsa`, or `pmav`.
- `fetch_kegg_gene_list()` pulls KEGG gene records from `list/<organism_code>`.
- `fetch_entrez_mapping()` pulls KEGG-to-NCBI GeneID mappings from `conv/ncbi-geneid/<organism_code>`.
- `build_mapping()` builds the rich KEGG-based mapping table.
- `build_ncbi_gene_mapping()` is the fallback for cases like `pmav`, where KEGG recognizes the organism but does not return a KEGG gene list.
- `write_mapping_files()` writes both the rich mapping and the PyGAGE-compatible `egSymb` file.

The rich output keeps:

```text
organism_code
kegg_gene_id
kegg_gene_number
entrez_id
symbol
description
has_entrez_id
source
```

The PyGAGE-compatible output keeps:

```text
entrez_id
symbol
```

## 2. `bin/pygage-build_egsymb_mapping.py`

This is a new command-line script.

Purpose:

```text
Let users build an organism-specific mapping file from the terminal or from a notebook.
```

Mouse example:

```powershell
python bin\pygage-build_egsymb_mapping.py mmu `
  --output ..\changes-test\outputs\mmu_gene_mapping.tsv `
  --egsymb-output ..\changes-test\outputs\mmu_egSymb.tsv `
  --timeout 60
```

`pmav` example using the NCBI fallback:

```powershell
python bin\pygage-build_egsymb_mapping.py pmav `
  --fallback ncbi `
  --output ..\changes-test\outputs\pmav_gene_mapping.tsv `
  --egsymb-output ..\changes-test\outputs\pmav_egSymb.tsv `
  --timeout 60
```

Why fallback exists:

```text
KEGG recognizes pmav, but KEGG does not return records for list/pmav.
NCBI Gene does have Peromyscus maniculatus bairdii gene records.
```

## 3. `pyproject.toml`

Change:

```text
bin/pygage-build_egsymb_mapping.py
```

was added to the `script-files` list.

Why:

```text
If this change were packaged later, the new script would install with the other PyGAGE command-line tools.
```

## 4. `pygage_README.md`

This is a copy of the local PyGAGE README after adding a short usage section.

Added documentation explains:

- how to build organism-specific mappings,
- what output files are written,
- why rich mappings preserve KEGG IDs and descriptions,
- why some organisms may need fallback behavior.

## What Was Tested

Mouse (`mmu`) worked directly through KEGG:

```text
Trp53 -> 22059
Brca1 -> 12189
Egfr  -> 13649
```

`pmav` worked through NCBI fallback:

```text
Oaz3  -> 102915534
Snrpn -> 102910316
Epas1 -> 102906331
```

## Git Note

The real local PyGAGE clone is still ignored by the outer repo:

```text
pygage-dev/
```

This folder only stores review copies of the changed files and notes.
