#!/usr/bin/env python3
"""
Gene ID Conversion Utilities

Functions for converting between Entrez Gene IDs and official gene symbols
for human genes using the egSymb mapping data.

"""

import polars as pl
import argparse
import re
import time
import requests
from pathlib import Path
from typing import Dict, Optional, Union, List


class GeneIDConverter:
    """Class for converting between Entrez Gene IDs and gene symbols."""
    
    def __init__(self, mapping_file: Optional[Path] = None):
        """
        Initialize the converter with a mapping file.
        
        Args:
            mapping_file: Path to mapping file with columns [entrez_id, symbol]
                         If None, will look for default 'egSymb.csv' or 'egSymb.tsv'
        """
        self.mapping_df = None
        if mapping_file is None:
            mapping_file = Path(__file__).parent / "data" / "egSymb.tsv"
        self.load_mapping(mapping_file)
    
    def load_mapping(self, mapping_file: Path):
        """
        Load gene ID to symbol mapping.
        
        Args:
            mapping_file: Path to mapping file (CSV or TSV)
        """

        if mapping_file.suffix == '.csv':
            self.mapping_df = pl.read_csv(mapping_file, infer_schema=False)
        elif mapping_file.suffix in ['.tsv', '.txt']:
            self.mapping_df = pl.read_csv(mapping_file, infer_schema=False, separator='\t')
        else:
            raise ValueError(f"Unsupported file format: {mapping_file.suffix}")
        
        # Ensure columns are named correctly
        if self.mapping_df.shape[1] >= 2:
            self.mapping_df = self.mapping_df.select([
                pl.col(self.mapping_df.columns[0]).alias('entrez_id'),
                pl.col(self.mapping_df.columns[1]).alias('symbol')
            ])
    
    def eg2sym(self, entrez_ids: Union[List[str], List[int], pl.Series]) -> List[Optional[str]]:
        """
        Convert Entrez Gene IDs to official gene symbols.
        
        Args:
            entrez_ids: List or Series of Entrez Gene IDs
            
        Returns:
            List of gene symbols (None for missing IDs)
        """

        if self.mapping_df is None:
            raise ValueError("Mapping data not loaded. Call load_mapping() first.")
        
        # Convert to list if needed
        if isinstance(entrez_ids, pl.Series):
            entrez_ids = entrez_ids.to_list()
        
        # Convert to strings for consistent matching
        entrez_ids = [str(x) for x in entrez_ids]
        
        # Create lookup DataFrame
        lookup_df = pl.DataFrame({'entrez_id': entrez_ids})
        
        # Join with mapping
        result = lookup_df.join(
            self.mapping_df.select(['entrez_id', 'symbol']),
            on='entrez_id',
            how='left'
        )
        result = result.select([
                pl.col(result.columns[0]).alias('input'),
                pl.col(result.columns[1]).alias('output')
            ])
        
        return result
    
    def sym2eg(self, symbols: Union[List[str], pl.Series]) -> List[Optional[str]]:
        """
        Convert official gene symbols to Entrez Gene IDs.
        
        Args:
            symbols: List or Series of gene symbols
            
        Returns:
            List of Entrez Gene IDs (None for missing symbols)
        """
        if self.mapping_df is None:
            raise ValueError("Mapping data not loaded. Call load_mapping() first.")
        
        # Convert to list if needed
        if isinstance(symbols, pl.Series):
            symbols = symbols.to_list()
        
        # Create lookup DataFrame
        lookup_df = pl.DataFrame({'symbol': symbols})
        
        # Join with mapping
        result = lookup_df.join(
            self.mapping_df.select(['symbol', 'entrez_id']),
            on='symbol',
            how='left'
        )
        result = result.select([
                pl.col(result.columns[0]).alias('input'),
                pl.col(result.columns[1]).alias('output')
            ])
        
        return result


class KEGGOrganismGeneMapper:
    """Build organism-specific KEGG gene to Entrez/symbol mapping tables."""

    KEGG_REST_BASE = "https://rest.kegg.jp"

    def __init__(self, timeout: int = 30, retries: int = 3, sleep_seconds: float = 0.5):
        self.timeout = timeout
        self.retries = retries
        self.sleep_seconds = sleep_seconds

    def _get_text(self, endpoint: str, base_url: str = KEGG_REST_BASE, params: Optional[Dict] = None) -> str:
        """Fetch a text endpoint with small retry handling."""
        url = f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"
        last_error = None

        for attempt in range(1, self.retries + 1):
            try:
                response = requests.get(url, params=params, timeout=self.timeout)
                if response.status_code == 200:
                    return response.text
                last_error = RuntimeError(
                    f"KEGG request failed for {url}: HTTP {response.status_code}"
                )
            except requests.RequestException as exc:
                last_error = exc

            if attempt < self.retries:
                time.sleep(self.sleep_seconds * attempt)

        raise RuntimeError(f"Could not fetch KEGG endpoint after {self.retries} attempts: {url}") from last_error

    def _get_json(self, url: str, params: Optional[Dict] = None) -> Dict:
        """Fetch a JSON endpoint with small retry handling."""
        last_error = None

        for attempt in range(1, self.retries + 1):
            try:
                response = requests.get(url, params=params, timeout=self.timeout)
                if response.status_code == 200:
                    return response.json()
                last_error = RuntimeError(
                    f"Request failed for {url}: HTTP {response.status_code}"
                )
            except (requests.RequestException, ValueError) as exc:
                last_error = exc

            if attempt < self.retries:
                time.sleep(self.sleep_seconds * attempt)

        raise RuntimeError(f"Could not fetch JSON endpoint after {self.retries} attempts: {url}") from last_error

    def list_organisms(self) -> pl.DataFrame:
        """Return KEGG organism metadata."""
        text = self._get_text("list/organism")
        rows = []

        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                rows.append({
                    "taxon_id": parts[0],
                    "organism_code": parts[1],
                    "organism_name": parts[2],
                    "lineage": parts[3],
                })

        return pl.DataFrame(rows)

    def validate_organism_code(self, organism_code: str) -> Dict[str, str]:
        """Validate an organism code against KEGG and return its metadata."""
        organisms = self.list_organisms()
        match = organisms.filter(pl.col("organism_code") == organism_code)

        if match.height == 1:
            return match.row(0, named=True)

        suggestions = (
            organisms
            .filter(
                pl.col("organism_code").str.contains(organism_code, literal=True)
                | pl.col("organism_name").str.to_lowercase().str.contains(organism_code.lower(), literal=True)
            )
            .select(["organism_code", "organism_name"])
            .head(10)
        )

        suggestion_text = ""
        if suggestions.height > 0:
            suggestion_text = "\nPossible matches:\n" + "\n".join(
                f"  {row['organism_code']}\t{row['organism_name']}"
                for row in suggestions.iter_rows(named=True)
            )

        raise ValueError(f"Unknown KEGG organism code: {organism_code}{suggestion_text}")

    def fetch_kegg_gene_list(self, organism_code: str) -> Dict[str, Dict[str, str]]:
        """Fetch KEGG gene descriptions for an organism."""
        try:
            text = self._get_text(f"list/{organism_code}")
        except RuntimeError as exc:
            raise RuntimeError(
                f"KEGG organism code '{organism_code}' was recognized, but KEGG did not "
                f"return a gene list for endpoint list/{organism_code}. This can happen "
                "for organisms without KEGG gene records. Try another organism code or "
                "use a local annotation file if Entrez/symbol mappings are needed."
            ) from exc

        genes = {}

        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            kegg_gene_id = parts[0]
            description = "\t".join(parts[1:])
            symbol_source = parts[-1]
            gene_number = kegg_gene_id.split(":", 1)[1] if ":" in kegg_gene_id else kegg_gene_id
            symbol = self._parse_symbol(symbol_source, gene_number)

            genes[kegg_gene_id] = {
                "organism_code": organism_code,
                "kegg_gene_id": kegg_gene_id,
                "kegg_gene_number": gene_number,
                "symbol": symbol,
                "description": description,
            }

        return genes

    def fetch_entrez_mapping(self, organism_code: str) -> Dict[str, str]:
        """Fetch KEGG gene to NCBI GeneID mapping for an organism."""
        try:
            text = self._get_text(f"conv/ncbi-geneid/{organism_code}")
        except RuntimeError:
            return {}

        mapping = {}

        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) != 2:
                continue

            kegg_gene_id = None
            entrez_id = None
            for part in parts:
                if part.startswith(f"{organism_code}:"):
                    kegg_gene_id = part
                elif part.startswith("ncbi-geneid:"):
                    entrez_id = part.replace("ncbi-geneid:", "", 1)

            if kegg_gene_id and entrez_id:
                mapping[kegg_gene_id] = entrez_id

        return mapping

    def build_mapping(self, organism_code: str, validate: bool = True) -> pl.DataFrame:
        """
        Build a KEGG organism gene mapping table.

        The returned table preserves KEGG identifiers even when an Entrez mapping
        is unavailable, which is common for some non-model organisms.
        """
        if validate:
            self.validate_organism_code(organism_code)

        genes = self.fetch_kegg_gene_list(organism_code)
        entrez_by_kegg_id = self.fetch_entrez_mapping(organism_code)

        rows = []
        for kegg_gene_id in sorted(genes):
            row = genes[kegg_gene_id].copy()
            row["entrez_id"] = entrez_by_kegg_id.get(kegg_gene_id)
            row["has_entrez_id"] = row["entrez_id"] is not None
            rows.append(row)

        return pl.DataFrame(rows).select([
            "organism_code",
            "kegg_gene_id",
            "kegg_gene_number",
            "entrez_id",
            "symbol",
            "description",
            "has_entrez_id",
        ])

    def build_ncbi_gene_mapping(
        self,
        organism_code: str,
        organism_name: Optional[str] = None,
        retmax: int = 100000,
        batch_size: int = 500,
        validate: bool = True,
    ) -> pl.DataFrame:
        """
        Build an Entrez/symbol mapping from NCBI Gene.

        This fallback is useful when KEGG recognizes an organism code but does
        not expose gene records through list/<organism_code>.
        """
        if organism_name is None:
            if validate:
                organism = self.validate_organism_code(organism_code)
                organism_name = self._normalize_ncbi_organism_name(organism["organism_name"])
            else:
                raise ValueError("organism_name is required when validate=False")

        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        search = self._get_json(search_url, params={
            "db": "gene",
            "term": f'"{organism_name}"[Organism]',
            "retmode": "json",
            "retmax": retmax,
        })

        ids = search.get("esearchresult", {}).get("idlist", [])
        rows = []

        for start in range(0, len(ids), batch_size):
            batch_ids = ids[start:start + batch_size]
            summary = self._get_json(summary_url, params={
                "db": "gene",
                "id": ",".join(batch_ids),
                "retmode": "json",
            })
            result = summary.get("result", {})

            for gene_id in batch_ids:
                item = result.get(gene_id, {})
                symbol = item.get("name") or item.get("nomenclaturesymbol") or gene_id
                rows.append({
                    "organism_code": organism_code,
                    "kegg_gene_id": None,
                    "kegg_gene_number": None,
                    "entrez_id": gene_id,
                    "symbol": symbol,
                    "description": item.get("description", ""),
                    "has_entrez_id": True,
                    "source": "NCBI Gene",
                })

        if not rows:
            raise RuntimeError(f"No NCBI Gene records found for organism name: {organism_name}")

        return pl.DataFrame(rows).select([
            "organism_code",
            "kegg_gene_id",
            "kegg_gene_number",
            "entrez_id",
            "symbol",
            "description",
            "has_entrez_id",
            "source",
        ])

    def write_mapping_files(
        self,
        organism_code: str,
        output: Path,
        eg_symb_output: Optional[Path] = None,
        validate: bool = True,
        fallback: str = "none",
        ncbi_retmax: int = 100000,
    ) -> Dict[str, Path]:
        """Write a rich mapping TSV and an egSymb-compatible TSV."""
        try:
            mapping = self.build_mapping(organism_code, validate=validate)
        except RuntimeError:
            if fallback != "ncbi":
                raise
            mapping = self.build_ncbi_gene_mapping(
                organism_code=organism_code,
                retmax=ncbi_retmax,
                validate=validate,
            )
        output.parent.mkdir(parents=True, exist_ok=True)
        mapping.write_csv(output, separator="\t")

        if eg_symb_output is None:
            eg_symb_output = output.with_name(f"{output.stem}.egSymb.tsv")

        eg_symb_output.parent.mkdir(parents=True, exist_ok=True)
        eg_symb = (
            mapping
            .filter(pl.col("entrez_id").is_not_null())
            .select(["entrez_id", "symbol"])
            .unique(subset=["entrez_id", "symbol"], maintain_order=True)
        )
        eg_symb.write_csv(eg_symb_output, separator="\t", quote_style="never")

        return {
            "mapping": output,
            "egSymb": eg_symb_output,
        }

    @staticmethod
    def _parse_symbol(description: str, fallback: str) -> str:
        """
        Parse the most useful symbol-like field from a KEGG gene description.

        KEGG descriptions are not perfectly uniform across organisms. The first
        semicolon-delimited field usually contains gene names or a locus tag.
        """
        head = description.split(";", 1)[0].strip()
        if not head:
            return fallback

        first_symbol = head.split(",", 1)[0].strip()
        return first_symbol or fallback

    @staticmethod
    def _normalize_ncbi_organism_name(organism_name: str) -> str:
        """Remove KEGG common-name parentheticals before NCBI organism search."""
        return re.sub(r"\s*\([^)]*\)\s*$", "", organism_name).strip()
