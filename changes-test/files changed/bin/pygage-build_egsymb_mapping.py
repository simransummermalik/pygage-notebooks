#!/usr/bin/env python3
"""
Build organism-specific KEGG gene mapping files for PyGAGE.

This creates:
1. a rich TSV with KEGG IDs, Entrez IDs, symbols, and descriptions
2. an egSymb-compatible TSV with only entrez_id and symbol
"""

import argparse
import sys
from pathlib import Path

repo_lib = Path(__file__).resolve().parents[1] / "lib"

try:
    from pygage.gene_id_utils import KEGGOrganismGeneMapper
except ImportError:
    if repo_lib.exists():
        sys.path.insert(0, str(repo_lib))
    from gene_id_utils import KEGGOrganismGeneMapper


def main():
    parser = argparse.ArgumentParser(
        description="Build an organism-specific egSymb-style mapping from KEGG."
    )
    parser.add_argument(
        "organism_code",
        help="KEGG organism code, for example hsa, mmu, eco, or pmav.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Rich mapping TSV output path. Default: <organism_code>_gene_mapping.tsv",
    )
    parser.add_argument(
        "--egsymb-output",
        type=Path,
        help="egSymb-compatible TSV output path. Default: <output stem>.egSymb.tsv",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip organism-code validation against KEGG's organism list.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Request timeout in seconds. Default: 30.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of KEGG request attempts. Default: 3.",
    )
    parser.add_argument(
        "--fallback",
        choices=["none", "ncbi"],
        default="none",
        help="Fallback source if KEGG has no gene list. Default: none.",
    )
    parser.add_argument(
        "--ncbi-retmax",
        type=int,
        default=100000,
        help="Maximum NCBI Gene records to retrieve when --fallback ncbi is used. Default: 100000.",
    )

    args = parser.parse_args()

    output = args.output
    if output is None:
        output = Path(f"{args.organism_code}_gene_mapping.tsv")

    mapper = KEGGOrganismGeneMapper(timeout=args.timeout, retries=args.retries)
    try:
        written = mapper.write_mapping_files(
            organism_code=args.organism_code,
            output=output,
            eg_symb_output=args.egsymb_output,
            validate=not args.no_validate,
            fallback=args.fallback,
            ncbi_retmax=args.ncbi_retmax,
        )
    except (RuntimeError, ValueError) as exc:
        parser.exit(status=1, message=f"Error: {exc}\n")

    print(f"Rich mapping written to: {written['mapping']}")
    print(f"egSymb-compatible mapping written to: {written['egSymb']}")


if __name__ == "__main__":
    main()
