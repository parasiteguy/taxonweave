"""
Tests for TaxonWeave taxonomic reconciliation logic.

These tests do not contact WoRMS, GenBank, or any other
external service.
"""

import pandas as pd

from taxonweave.query import (
    normalize_taxon_name,
    reconcile_genbank_names,
)


# ============================================================
# TAXON NAME NORMALIZATION
# ============================================================

def test_normalize_taxon_name_lowercases():
    assert (
        normalize_taxon_name("Alitta succinea")
        == "alitta succinea"
    )


def test_normalize_taxon_name_removes_extra_whitespace():
    assert (
        normalize_taxon_name("  Alitta   succinea  ")
        == "alitta succinea"
    )


def test_normalize_taxon_name_handles_none():
    assert normalize_taxon_name(None) is None


# ============================================================
# SYNTHETIC TEST DATA
# ============================================================

def make_test_records():
    """
    Construct synthetic GenBank-like records representing
    the four possible TaxonWeave reconciliation outcomes.
    """

    return pd.DataFrame(
        {
            "accession": [
                "TEST001",
                "TEST002",
                "TEST003",
                "TEST004",
            ],
            "deposited_name": [
                "Alitta succinea",
                "Neanthes succinea",
                "Nereis mysteryensis",
                None,
            ],
        }
    )


def make_test_synonyms():
    """
    Construct synthetic WoRMS-style synonym records.
    """

    return [
        {
            "scientificname": "Neanthes succinea"
        },
        {
            "scientificname": "Nereis succinea"
        },
    ]


# ============================================================
# RECONCILIATION RELATIONSHIPS
# ============================================================

def test_reconciliation_relationships():

    records = make_test_records()
    synonyms = make_test_synonyms()

    reconciled = reconcile_genbank_names(
        genbank=records,
        accepted_name="Alitta succinea",
        synonyms=synonyms,
    )

    results = dict(
        zip(
            reconciled["accession"],
            reconciled["name_relationship"],
        )
    )

    assert results["TEST001"] == "ACCEPTED_NAME"

    assert results["TEST002"] == "WORMS_SYNONYM"

    assert results["TEST003"] == "UNRESOLVED_NAME"

    assert results["TEST004"] == "MISSING_NAME"


# ============================================================
# RECONCILIATION STATUS
# ============================================================

def test_reconciliation_status():

    records = make_test_records()
    synonyms = make_test_synonyms()

    reconciled = reconcile_genbank_names(
        genbank=records,
        accepted_name="Alitta succinea",
        synonyms=synonyms,
    )

    results = dict(
        zip(
            reconciled["accession"],
            reconciled["reconciliation_status"],
        )
    )

    assert results["TEST001"] == "RECONCILED"

    assert (
        results["TEST002"]
        == "RECONCILED_NOMENCLATURALLY"
    )

    assert results["TEST003"] == "UNRESOLVED"

    assert results["TEST004"] == "UNRESOLVED"


# ============================================================
# CONSERVATIVE STRING MATCHING
# ============================================================

def test_reconciliation_ignores_case_and_extra_whitespace():

    records = pd.DataFrame(
        {
            "accession": [
                "TEST005",
                "TEST006",
            ],
            "deposited_name": [
                "  ALITTA   SUCCINEA  ",
                "  NEANTHES   SUCCINEA  ",
            ],
        }
    )

    reconciled = reconcile_genbank_names(
        genbank=records,
        accepted_name="Alitta succinea",
        synonyms=make_test_synonyms(),
    )

    relationships = list(
        reconciled["name_relationship"]
    )

    assert relationships == [
        "ACCEPTED_NAME",
        "WORMS_SYNONYM",
    ]


# ============================================================
# EMPTY GENBANK DATASET
# ============================================================

def test_empty_genbank_dataset():

    records = pd.DataFrame()

    reconciled = reconcile_genbank_names(
        genbank=records,
        accepted_name="Alitta succinea",
        synonyms=make_test_synonyms(),
    )

    assert reconciled.empty

    assert (
        "name_relationship"
        in reconciled.columns
    )

    assert (
        "reconciliation_status"
        in reconciled.columns
    )