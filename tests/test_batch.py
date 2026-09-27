"""
Tests for TaxonWeave batch-query support.
"""

import pandas as pd
import pytest

import taxonweave.batch as batch


class FakeReport:
    """
    Minimal stand-in for TaxonWeaveReport.
    """

    def __init__(
        self,
        query_name,
        accepted_name,
    ):
        self.query_name = query_name
        self.accepted_name = accepted_name

    def summary_dict(self):
        return {
            "query_name":
                self.query_name,

            "accepted_name":
                self.accepted_name,

            "aphia_id":
                12345,

            "genbank_records":
                10,
        }

    def export(
        self,
        directory_name=None,
    ):
        directory_name.mkdir(
            parents=True,
            exist_ok=True,
        )

        (
            directory_name
            / "summary.csv"
        ).write_text(
            "test\n",
            encoding="utf-8",
        )


def test_safe_directory_name():
    assert (
        batch.safe_directory_name(
            "Ficopomatus enigmaticus"
        )
        ==
        "Ficopomatus_enigmaticus"
    )


def test_read_batch_file(tmp_path):

    input_file = (
        tmp_path
        / "species.csv"
    )

    pd.DataFrame({
        "scientific_name": [
            "Alitta succinea",
            "Ficopomatus enigmaticus",
        ]
    }).to_csv(
        input_file,
        index=False,
    )

    names = batch.read_batch_file(
        input_file
    )

    assert names == [
        "Alitta succinea",
        "Ficopomatus enigmaticus",
    ]


def test_read_batch_requires_scientific_name_column(
    tmp_path
):

    input_file = (
        tmp_path
        / "species.csv"
    )

    pd.DataFrame({
        "species": [
            "Alitta succinea"
        ]
    }).to_csv(
        input_file,
        index=False,
    )

    with pytest.raises(
        ValueError,
        match="scientific_name",
    ):
        batch.read_batch_file(
            input_file
        )


def test_batch_continues_after_failure(
    tmp_path,
    monkeypatch,
):

    input_file = (
        tmp_path
        / "species.csv"
    )

    output_directory = (
        tmp_path
        / "results"
    )

    pd.DataFrame({
        "scientific_name": [
            "Alitta succinea",
            "Invalid taxon",
            "Ficopomatus enigmaticus",
        ]
    }).to_csv(
        input_file,
        index=False,
    )

    def fake_query_species(
        scientific_name
    ):

        if scientific_name == "Invalid taxon":
            raise ValueError(
                "WoRMS could not resolve taxon."
            )

        return FakeReport(
            query_name=scientific_name,
            accepted_name=scientific_name,
        )

    monkeypatch.setattr(
        batch,
        "query_species",
        fake_query_species,
    )

    summary = batch.run_batch(
        input_file=input_file,
        output_directory=output_directory,
    )

    assert len(summary) == 3

    assert list(
        summary["batch_status"]
    ) == [
        "SUCCESS",
        "FAILED",
        "SUCCESS",
    ]

    assert (
        output_directory
        / "batch_summary.csv"
    ).exists()

    assert (
        output_directory
        / "Alitta_succinea"
        / "summary.csv"
    ).exists()

    assert (
        output_directory
        / "Ficopomatus_enigmaticus"
        / "summary.csv"
    ).exists()

    failed_row = summary[
        summary["query_name"]
        == "Invalid taxon"
    ].iloc[0]

    assert (
        "ValueError"
        in failed_row["batch_error"]
    )
def test_percentage():

    assert (
        batch.percentage(
            50,
            100,
        )
        == 50.0
    )

    assert (
        batch.percentage(
            1,
            3,
        )
        == 33.3
    )

    assert (
        batch.percentage(
            0,
            0,
        )
        == 0.0
    )


def test_add_comparative_metrics():

    summary = {
        "genbank_records": 200,
        "accepted_name_records": 180,
        "worms_synonym_records": 10,
        "unresolved_name_records": 10,
        "geographic_records": 150,
        "coordinate_records": 100,
        "collection_date_records": 120,
        "voucher_records": 50,
    }

    result = (
        batch.add_comparative_metrics(
            summary
        )
    )

    assert (
        result["accepted_name_pct"]
        == 90.0
    )

    assert (
        result["synonym_name_pct"]
        == 5.0
    )

    assert (
        result["unresolved_name_pct"]
        == 5.0
    )

    assert (
        result["geographic_metadata_pct"]
        == 75.0
    )

    assert (
        result["coordinate_metadata_pct"]
        == 50.0
    )

    assert (
        result["collection_date_pct"]
        == 60.0
    )

    assert (
        result["voucher_metadata_pct"]
        == 25.0
    )