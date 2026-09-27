"""
Tests for TaxonWeave molecular marker normalization.

These tests exercise the annotation-normalization logic without
contacting GenBank, NCBI, WoRMS, or any other external service.
"""

import pytest

from taxonweave.genbank import normalize_marker


# ============================================================
# STANDARD MITOCHONDRIAL MARKERS
# ============================================================

@pytest.mark.parametrize(
    "annotation, expected",
    [
        ("COI", "COI"),
        ("CO1", "COI"),
        ("cox1", "COI"),
        ("cytochrome c oxidase subunit I", "COI"),
        ("cytochrome c oxidase subunit 1", "COI"),

        ("COII", "COII"),
        ("cox2", "COII"),
        ("cytochrome c oxidase subunit II", "COII"),

        ("COIII", "COIII"),
        ("cox3", "COIII"),
        ("cytochrome c oxidase subunit III", "COIII"),

        ("cytb", "CytB"),
        ("cyt b", "CytB"),
        ("cob", "CytB"),
        ("cytochrome b", "CytB"),

        ("ATP6", "ATP6"),
        ("ATP synthase subunit 6", "ATP6"),

        ("ATP8", "ATP8"),
        ("ATP synthase subunit 8", "ATP8"),

        ("ND1", "ND1"),
        ("ND2", "ND2"),
        ("ND3", "ND3"),
        ("ND4", "ND4"),
        ("ND4L", "ND4L"),
        ("ND5", "ND5"),
        ("ND6", "ND6"),
    ],
)
def test_mitochondrial_protein_coding_markers(
    annotation,
    expected,
):
    normalized, category, status = normalize_marker(
        annotation
    )

    assert normalized == expected
    assert category == "mitochondrial_protein_coding"
    assert status == "NORMALIZED"


# ============================================================
# RIBOSOMAL RNA MARKERS
# ============================================================

@pytest.mark.parametrize(
    "annotation, expected, category",
    [
        (
            "12S ribosomal RNA",
            "12S",
            "mitochondrial_rRNA",
        ),
        (
            "16S ribosomal RNA",
            "16S",
            "mitochondrial_rRNA",
        ),
        (
            "18S ribosomal RNA",
            "18S",
            "nuclear_rRNA",
        ),
        (
            "28S ribosomal RNA",
            "28S",
            "nuclear_rRNA",
        ),
    ],
)
def test_specific_ribosomal_markers(
    annotation,
    expected,
    category,
):
    normalized, observed_category, status = (
        normalize_marker(annotation)
    )

    assert normalized == expected
    assert observed_category == category
    assert status == "NORMALIZED"


# ============================================================
# HISTONE H3
# ============================================================

def test_histone_h3():

    normalized, category, status = normalize_marker(
        "histone H3"
    )

    assert normalized == "H3"
    assert category == "nuclear_protein_coding"
    assert status == "NORMALIZED"


# ============================================================
# MITOCHONDRIAL tRNA
# ============================================================

@pytest.mark.parametrize(
    "annotation",
    [
        "trnN",
        "trnA",
        "trnN(aac)",
    ],
)
def test_mitochondrial_trna(annotation):

    normalized, category, status = normalize_marker(
        annotation
    )

    assert normalized == annotation
    assert category == "mitochondrial_tRNA"
    assert status == "RECOGNIZED"


# ============================================================
# OTHER RECOGNIZED GENES
# ============================================================

@pytest.mark.parametrize(
    "annotation",
    [
        "extracellular globin",
        "hemoglobin linker chain",
        "hemerythrin",
        "actin",
        "tubulin",
        "heat shock protein 70",
    ],
)
def test_other_recognized_genes(annotation):

    normalized, category, status = normalize_marker(
        annotation
    )

    assert normalized == annotation

    assert (
        category
        == "nuclear_or_nonmitochondrial_gene"
    )

    assert status == "RECOGNIZED"


# ============================================================
# AMBIGUOUS RIBOSOMAL ANNOTATIONS
# ============================================================

@pytest.mark.parametrize(
    "annotation",
    [
        "large subunit ribosomal RNA",
        "large subunit rRNA",
        "small subunit ribosomal RNA",
        "small subunit rRNA",
    ],
)
def test_generic_ribosomal_annotations_are_ambiguous(
    annotation,
):

    normalized, category, status = normalize_marker(
        annotation
    )

    assert normalized is None
    assert category == "ribosomal_RNA"
    assert status == "AMBIGUOUS"


# ============================================================
# UNRESOLVED ANNOTATIONS
# ============================================================

@pytest.mark.parametrize(
    "annotation",
    [
        None,
        "UNVERIFIED: organism mitochondrion sequence",
        "mystery sequence",
    ],
)
def test_unresolved_annotations(annotation):

    normalized, category, status = normalize_marker(
        annotation
    )

    assert normalized is None
    assert category == "unknown"
    assert status == "UNRESOLVED"


# ============================================================
# CASE INSENSITIVITY
# ============================================================

def test_marker_normalization_is_case_insensitive():

    normalized, category, status = normalize_marker(
        "CYTOCHROME C OXIDASE SUBUNIT I"
    )

    assert normalized == "COI"
    assert status == "NORMALIZED"