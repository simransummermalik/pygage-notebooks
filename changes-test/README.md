# Organism-Specific egSymb Mapping Test

This folder is a test area for the organism-specific gene mapping idea Andra described.

It does **not** mean the official PyGAGE GitHub repository or the `pip install pygage` package has been changed. This is local experimental work that tests whether the approach works before asking anyone to review or merge it.

## Why This Test Exists

PyGAGE already has an `egSymb.tsv`-style mapping file for converting between Entrez Gene IDs and gene symbols.

The existing mapping format is simple:

```text
entrez_id    symbol
```

Andra asked whether this could be made more flexible for different KEGG organism codes, for example mouse (`mmu`) or another organism code.

The main idea tested here is:

1. take a KEGG organism code,
2. pull KEGG gene records for that organism,
3. pull KEGG-to-NCBI GeneID conversion when KEGG provides it,
4. parse a useful symbol-like field,
5. write a rich mapping table,
6. also write a two-column `egSymb.tsv`-compatible table that PyGAGE can already use.

## Important Colab / `pip install pygage` Note

In Colab, this still works with the normal PyGAGE dependency style:

```python
!pip install pygage
```

But the new mapping generator is **not part of the released PyPI package yet**.

So right now there are two different cases:

### Current released package

```python
!pip install pygage
```

This installs the official released PyGAGE package. It will have the current PyGAGE tools, but it will **not** include this new organism-specific mapping script unless the change is merged and released later.

### Local test version

This test notebook imports the local experimental code from:

```text
pygage-dev/lib/gene_id_utils.py
```

That is why the notebook adds this folder to Python's import path:

```python
sys.path.insert(0, str(pygage_dev / "lib"))
```

This lets the notebook test the new mapping code before it exists in the public `pip install pygage` package.

If this change eventually gets pushed to a branch or merged upstream, then Colab could install it from GitHub instead, for example:

```python
!pip install git+https://github.com/<user-or-org>/pygage.git@<branch-name>
```

After an official release, plain `pip install pygage` would be enough.

## Notebook

Main notebook:

```text
changes-test/test_organism_specific_egsymb_mapping.ipynb
```

The notebook tests:

- dependency setup,
- local import of the experimental mapping class,
- mouse mapping generation using KEGG organism code `mmu`,
- output file inspection,
- conversion with PyGAGE's existing `GeneIDConverter`,
- and the `pmav` NCBI fallback path.

## Setup Cell Added to Fix the Current Error

The notebook error was:

```text
ModuleNotFoundError: No module named 'polars'
```

That happened because `gene_id_utils.py` imports `polars`, but the notebook kernel did not have `polars` installed.

The fix is this setup cell near the top of the notebook:

```python
%pip install -q polars requests
```

After running that cell, restart the kernel if VS Code asks you to, then run the notebook from the top.

## Experimental Code Being Tested

The main experimental code is in:

```text
pygage-dev/lib/gene_id_utils.py
```

The added class is:

```python
class KEGGOrganismGeneMapper:
    """Build organism-specific KEGG gene to Entrez/symbol mapping tables."""
```

### Why This Class Was Added

The original `GeneIDConverter` can already convert IDs if a mapping file exists.

The missing piece was a way to **generate** a mapping file for a specific organism code.

So the test adds a mapper that can build a mapping file first, then pass that generated file into the existing converter.

## Main Code Snippets

### 1. Validate KEGG organism codes

Location:

```text
pygage-dev/lib/gene_id_utils.py
```

Purpose:

Checks whether a code like `mmu` exists in KEGG before trying to build a mapping.

```python
def validate_organism_code(self, organism_code: str) -> Dict[str, str]:
    """Validate an organism code against KEGG and return its metadata."""
    organisms = self.list_organisms()
    match = organisms.filter(pl.col("organism_code") == organism_code)

    if match.height == 1:
        return match.row(0, named=True)
```

### 2. Fetch KEGG gene records

Purpose:

Pulls gene records from:

```text
https://rest.kegg.jp/list/<organism_code>
```

Example for mouse:

```text
https://rest.kegg.jp/list/mmu
```

Code shape:

```python
def fetch_kegg_gene_list(self, organism_code: str) -> Dict[str, Dict[str, str]]:
    """Fetch KEGG gene descriptions for an organism."""
    text = self._get_text(f"list/{organism_code}")
```

The parser keeps:

```text
organism_code
kegg_gene_id
kegg_gene_number
symbol
description
```

### 3. Fetch Entrez mappings

Purpose:

Pulls KEGG-to-NCBI GeneID mappings from:

```text
https://rest.kegg.jp/conv/ncbi-geneid/<organism_code>
```

Example for mouse:

```text
https://rest.kegg.jp/conv/ncbi-geneid/mmu
```

Code shape:

```python
def fetch_entrez_mapping(self, organism_code: str) -> Dict[str, str]:
    """Fetch KEGG gene to NCBI GeneID mapping for an organism."""
    try:
        text = self._get_text(f"conv/ncbi-geneid/{organism_code}")
    except RuntimeError:
        return {}
```

The `return {}` fallback is intentional. Some organisms have KEGG gene records but do not provide NCBI GeneID conversion. In that case, the rich KEGG mapping can still be useful, but the strict `egSymb` output may be empty or partial.

### 4. Fall back to NCBI Gene when KEGG has no gene list

Purpose:

Supports cases like `pmav`, where KEGG recognizes the organism code but does not return gene records from:

```text
https://rest.kegg.jp/list/pmav
```

The fallback uses the KEGG organism name to query NCBI Gene.

Code shape:

```python
def build_ncbi_gene_mapping(
    self,
    organism_code: str,
    organism_name: Optional[str] = None,
    retmax: int = 100000,
    batch_size: int = 500,
    validate: bool = True,
) -> pl.DataFrame:
```

It writes rows with:

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

For NCBI fallback rows, the KEGG-specific ID fields are empty because the mapping came from NCBI Gene, not from KEGG gene records.

### 5. Write both output formats

Purpose:

Writes one rich file and one PyGAGE-compatible file.

Code shape:

```python
def write_mapping_files(
    self,
    organism_code: str,
    output: Path,
    eg_symb_output: Optional[Path] = None,
    validate: bool = True,
) -> Dict[str, Path]:
```

Rich output columns:

```text
organism_code
kegg_gene_id
kegg_gene_number
entrez_id
symbol
description
has_entrez_id
```

PyGAGE-compatible output columns:

```text
entrez_id
symbol
```

This is the key compatibility point. PyGAGE's existing `GeneIDConverter` already expects:

```text
entrez_id    symbol
```

So the generated `mmu_egSymb.tsv` can be used with the existing converter.

## Command-Line Script Being Tested

Script:

```text
pygage-dev/bin/pygage-build_egsymb_mapping.py
```

Purpose:

Makes the mapping generator runnable from the terminal or notebook.

Mouse test command:

```powershell
python pygage-dev\bin\pygage-build_egsymb_mapping.py mmu `
  --output changes-test\outputs\mmu_gene_mapping.tsv `
  --egsymb-output changes-test\outputs\mmu_egSymb.tsv `
  --timeout 60
```

Expected output files:

```text
changes-test/outputs/mmu_gene_mapping.tsv
changes-test/outputs/mmu_egSymb.tsv
```

`pmav` fallback test command:

```powershell
python pygage-dev\bin\pygage-build_egsymb_mapping.py pmav `
  --fallback ncbi `
  --ncbi-retmax 25 `
  --output changes-test\outputs\pmav_gene_mapping.test25.tsv `
  --egsymb-output changes-test\outputs\pmav_egSymb.test25.tsv `
  --timeout 60
```

The `--ncbi-retmax 25` option is only for a quick test. For a full `pmav` mapping, remove that small limit or set it higher.

## Mouse Test Result

The `mmu` test worked.

Verified conversions:

```text
Trp53 -> 22059
Brca1 -> 12189
Egfr  -> 13649

22059 -> Trp53
12189 -> Brca1
13649 -> Egfr
```

This means the generated mouse `egSymb` file works with PyGAGE's existing `GeneIDConverter`.

## `pmav` Test Result

`pmav` is recognized by KEGG as:

```text
Peromyscus maniculatus bairdii
```

KEGG currently returns an error for:

```text
https://rest.kegg.jp/list/pmav
```

So KEGG alone cannot build the mapping.

The new fallback can still build an `egSymb`-style file for `pmav` from NCBI Gene:

```text
entrez_id    symbol
102915534    Oaz3
102910316    Snrpn
102906331    Epas1
102908193    Hoxd13
```

This is important because the future-proof version should not assume every KEGG organism code has complete gene-list and Entrez-conversion support. The practical behavior is:

1. try KEGG first,
2. if KEGG has no gene list, optionally fall back to NCBI Gene,
3. write a clear source column in the rich mapping so users know where the mapping came from.

## Git Safety Note

The local `pygage-dev/` clone is ignored by the outer walkthrough repo through:

```text
.gitignore
```

with:

```text
pygage-dev/
```

That means pushing this walkthrough repo does **not** push the local PyGAGE clone.

