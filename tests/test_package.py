"""
Basic tests for the TaxonWeave package.
"""

import taxonweave


def test_package_imports():
    """TaxonWeave should import successfully."""
    assert taxonweave is not None


def test_version():
    """TaxonWeave should expose its package version."""
    assert taxonweave.__version__ == "0.2.0"


def test_query_species_is_available():
    """query_species should be exposed through the public API."""
    assert callable(taxonweave.query_species)