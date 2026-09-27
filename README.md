# TaxonWeave

**TaxonWeave** is a Python framework for provenance-aware reconciliation of taxonomic concepts and molecular sequence records.

Taxonomic names change through time, but molecular records deposited in public databases retain the names and annotations associated with their original submissions. As a result, sequence data relevant to a currently accepted species may be distributed across accepted names, historical synonyms, and other nomenclatural combinations.

TaxonWeave connects taxonomic information from the **World Register of Marine Species (WoRMS)** with nucleotide records from **GenBank** to make those relationships explicit and reproducible.

TaxonWeave is currently under active development.

**Current version:** 0.1.0

## What TaxonWeave does

Given a scientific name, TaxonWeave:

1. resolves the supplied name through WoRMS;
2. identifies the currently accepted taxonomic concept and AphiaID;
3. retrieves nomenclatural synonyms associated with that concept;
4. searches GenBank using the accepted name and WoRMS synonyms;
5. tracks which taxonomic search name retrieved each GenBank record;
6. deduplicates records retrieved through multiple names;
7. preserves the taxonomic name originally represented on each GenBank record;
8. reconciles deposited names against the current WoRMS concept;
9. extracts molecular, geographic, collection, and specimen-associated metadata;
10. conservatively normalizes commonly used molecular-marker annotations; and
11. produces a structured summary that can be exported as CSV files.

The objective is not to replace taxonomic judgment. TaxonWeave provides a reproducible evidence-reconciliation layer connecting changing taxonomic concepts with molecular records.

## Search provenance and deposited names

A central design principle of TaxonWeave is the distinction between **search provenance** and **original taxonomic assertion**.

For each GenBank record, TaxonWeave distinguishes:

- `found_via_name` — the taxonomic name used in a GenBank search that retrieved the record;
- `deposited_name` — the organism name represented on the GenBank record itself;
- `accepted_name` — the currently accepted taxonomic concept resolved through WoRMS.

These fields should not be treated as interchangeable.

For example, a historical synonym may retrieve a GenBank record even when that record is currently represented in GenBank under the accepted species name. TaxonWeave preserves both pieces of information rather than interpreting the search term as the deposited identification.

This distinction allows synonym-mediated retrieval to be documented without rewriting the original database assertion.

## Taxonomic reconciliation

TaxonWeave currently classifies deposited GenBank names relative to the WoRMS concept as:

- `ACCEPTED_NAME`
- `WORMS_SYNONYM`
- `UNRESOLVED_NAME`
- `MISSING_NAME`

Associated reconciliation states distinguish records that match the accepted name, records that can be reconciled nomenclaturally through WoRMS, and records requiring further investigation.

A WoRMS synonym is treated as a **nomenclatural assertion**, not as independent biological evidence that two sampled organisms are conspecific.

## Molecular-marker normalization

GenBank records contain heterogeneous gene and product annotations. TaxonWeave preserves the raw annotations while also applying conservative normalization rules for commonly encountered markers.

Examples include:

- COI, COII, and COIII
- CytB
- ATP6 and ATP8
- ND1–ND6 and ND4L
- 12S and 16S mitochondrial rRNA
- 18S and 28S nuclear rRNA
- histone H3

Normalization status is explicitly recorded as:

- `NORMALIZED`
- `RECOGNIZED`
- `AMBIGUOUS`
- `UNRESOLVED`

TaxonWeave deliberately avoids making unsupported inferences. For example, a generic annotation such as `large subunit ribosomal RNA` is retained as ambiguous rather than automatically being interpreted as 28S.

## Installation

TaxonWeave currently requires Python 3.9 or later.

Clone the repository:

```bash
git clone https://github.com/parasiteguy/taxonweave.git
cd taxonweave
```

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install TaxonWeave in editable mode:

```bash
pip install -e .
```

## Command-line usage

A species can be queried directly from the command line:

```bash
taxonweave query "Ficopomatus enigmaticus"
```

TaxonWeave will resolve the taxonomic concept through WoRMS, search GenBank across the accepted name and associated synonyms, reconcile the resulting records, and print a summary report.

An obsolete or historical name can also be supplied:

```bash
taxonweave query "Mercierella enigmatica"
```

TaxonWeave first resolves the supplied name through WoRMS before constructing the molecular search.

### Exporting results

Results can be exported as CSV files:

```bash
taxonweave query "Ficopomatus enigmaticus" --export
```

A specific output directory can be supplied:

```bash
taxonweave query "Ficopomatus enigmaticus" \
    --export \
    --output ficopomatus_results
```

## Python usage

TaxonWeave can also be used directly from Python:

```python
from taxonweave import query_species

report = query_species("Ficopomatus enigmaticus")

report.summary()
```

The report can then be exported:

```python
report.export("ficopomatus_results")
```

## Exported files

A TaxonWeave report currently contains five CSV files:

| File | Contents |
|---|---|
| `taxonomy.csv` | Accepted WoRMS taxonomic concept and classification |
| `synonyms.csv` | WoRMS nomenclatural synonyms associated with the concept |
| `genbank_searches.csv` | GenBank search names and retrieval counts |
| `sequences.csv` | Reconciled GenBank records and associated molecular/specimen metadata |
| `summary.csv` | Summary statistics for the complete query |

The `sequences.csv` file preserves both raw database annotations and TaxonWeave-derived reconciliation fields so that downstream analyses can distinguish source data from interpreted or normalized information.

## Example: *Ficopomatus enigmaticus*

The invasive serpulid polychaete *Ficopomatus enigmaticus* provides a useful example of the TaxonWeave workflow.

WoRMS recognizes historical names including *Mercierella enigmatica* and *Phycopomatus enigmaticus*. Searching these names independently can retrieve molecular records associated with the same contemporary taxonomic concept.

TaxonWeave resolves the nomenclatural relationships, records which search names recover each sequence, deduplicates overlapping GenBank results, and preserves the taxonomic identification represented on each individual record.

This makes it possible to distinguish the history of the **name used to discover a record** from the taxonomic assertion associated with the **record itself**.

## Current scope

TaxonWeave v0.1.0 currently focuses on:

**WoRMS → taxonomic concept and nomenclatural history**

**GenBank → molecular records and associated metadata**

Occurrence-data integration is not part of the core v0.1.0 workflow.

TaxonWeave currently relies on live external database services. Results may therefore change as WoRMS and GenBank records are added, revised, or reannotated.

Taxonomic reconciliation should be interpreted as evidence organization rather than automated species delimitation or taxonomic revision.

## Development

Run the automated test suite with:

```bash
pytest -v
```

The test suite currently covers package integrity, taxonomic reconciliation logic, molecular-marker normalization, and command-line argument parsing.

Live database queries are conceptually distinct from the deterministic unit tests because external API content can change independently of TaxonWeave.

## Roadmap

Planned development includes:

- batch processing of multiple taxonomic names;
- expanded provenance reporting;
- additional reconciliation diagnostics;
- improved handling of specimen and voucher metadata;
- stable command-line workflows;
- programmatic report access for downstream biodiversity-informatics analyses; and
- a web interface built on the same TaxonWeave scientific engine.

## Citation

TaxonWeave is under active development. Formal citation information will be added with the first archived software release.

## License

TaxonWeave is released under the MIT License.

Copyright (c) 2026 Andrew Davinack