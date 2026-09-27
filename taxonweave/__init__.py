"""
TaxonWeave

Provenance-aware reconciliation of taxonomic concepts
and molecular sequence records.
"""

from .query import query_species

__version__ = "0.1.0"

__all__ = [
    "query_species",
]