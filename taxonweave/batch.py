"""
Batch-query support for TaxonWeave.

This module provides an orchestration layer around query_species().
It does not implement independent taxonomic or molecular
reconciliation logic.

Batch summaries additionally contain derived percentage metrics
that facilitate comparison of taxonomic reconciliation and metadata
completeness across taxa.
"""

from pathlib import Path
import re

import pandas as pd

from taxonweave.query import query_species


def safe_directory_name(name):
    """
    Convert a taxonomic name into a filesystem-safe directory name.
    """

    text = str(name).strip()

    text = re.sub(
        r"[^\w.-]+",
        "_",
        text,
    )

    text = text.strip("._")

    return text or "taxon"


def read_batch_file(file_path):
    """
    Read scientific names from a CSV file.

    The input file must contain a column named 'scientific_name'.

    Returns
    -------
    list of str
        Scientific names in input order, with blank values removed.
        Duplicate names are retained because each input row represents
        a submitted query.
    """

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Batch input file not found: {file_path}"
        )

    dataframe = pd.read_csv(file_path)

    if "scientific_name" not in dataframe.columns:
        raise ValueError(
            "Batch input CSV must contain a "
            "'scientific_name' column."
        )

    names = []

    for value in dataframe["scientific_name"]:

        if pd.isna(value):
            continue

        name = str(value).strip()

        if name:
            names.append(name)

    if not names:
        raise ValueError(
            "Batch input CSV contains no scientific names."
        )

    return names


def percentage(count, total):
    """
    Calculate a percentage for a batch summary metric.

    Returns 0.0 when no GenBank records were retrieved.
    """

    if not total:
        return 0.0

    return round(
        (count / total) * 100,
        1,
    )


def add_comparative_metrics(summary):
    """
    Add derived percentage metrics to a TaxonWeave summary.

    These metrics are calculated from the raw counts produced by
    TaxonWeaveReport.summary_dict(). They do not alter the underlying
    taxonomic or molecular reconciliation results.
    """

    summary = summary.copy()

    total = summary.get(
        "genbank_records",
        0,
    )

    summary["accepted_name_pct"] = percentage(
        summary.get(
            "accepted_name_records",
            0,
        ),
        total,
    )

    summary["synonym_name_pct"] = percentage(
        summary.get(
            "worms_synonym_records",
            0,
        ),
        total,
    )

    summary["unresolved_name_pct"] = percentage(
        summary.get(
            "unresolved_name_records",
            0,
        ),
        total,
    )

    summary["geographic_metadata_pct"] = percentage(
        summary.get(
            "geographic_records",
            0,
        ),
        total,
    )

    summary["coordinate_metadata_pct"] = percentage(
        summary.get(
            "coordinate_records",
            0,
        ),
        total,
    )

    summary["collection_date_pct"] = percentage(
        summary.get(
            "collection_date_records",
            0,
        ),
        total,
    )

    summary["voucher_metadata_pct"] = percentage(
        summary.get(
            "voucher_records",
            0,
        ),
        total,
    )

    return summary


def run_batch(
    input_file,
    output_directory="taxonweave_batch",
):
    """
    Run TaxonWeave for every scientific name in a CSV file.

    Each successful taxon receives its own TaxonWeave report
    directory. A combined batch_summary.csv records the outcome
    of every submitted query.

    Failure of one taxon does not terminate the batch.

    Parameters
    ----------
    input_file : str or pathlib.Path
        CSV file containing a 'scientific_name' column.

    output_directory : str or pathlib.Path
        Directory in which batch results will be written.

    Returns
    -------
    pandas.DataFrame
        Combined batch summary.
    """

    names = read_batch_file(
        input_file
    )

    output_directory = Path(
        output_directory
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_rows = []

    total = len(names)

    print()
    print("=" * 70)
    print("TAXONWEAVE BATCH QUERY")
    print("=" * 70)
    print(f"Input taxa: {total}")
    print(
        "Output directory:",
        output_directory.resolve(),
    )

    for index, scientific_name in enumerate(
        names,
        start=1,
    ):

        print()
        print("-" * 70)
        print(
            f"[{index}/{total}] "
            f"{scientific_name}"
        )
        print("-" * 70)

        try:

            report = query_species(
                scientific_name
            )

            directory_base = (
                safe_directory_name(
                    report.accepted_name
                    or scientific_name
                )
            )

            taxon_directory = (
                output_directory
                / directory_base
            )

            # Avoid overwriting a previous result if the same
            # accepted concept occurs more than once in the input.
            if taxon_directory.exists():

                suffix = 2

                while (
                    output_directory
                    / f"{directory_base}_{suffix}"
                ).exists():

                    suffix += 1

                taxon_directory = (
                    output_directory
                    / f"{directory_base}_{suffix}"
                )

            report.export(
                directory_name=taxon_directory
            )

            row = report.summary_dict()

            row = add_comparative_metrics(
                row
            )

            row["batch_status"] = "SUCCESS"
            row["batch_error"] = None
            row["output_directory"] = (
                taxon_directory.name
            )

            summary_rows.append(
                row
            )

            print(
                f"Batch status: SUCCESS — "
                f"{scientific_name}"
            )

        except Exception as error:

            summary_rows.append({
                "query_name":
                    scientific_name,

                "accepted_name":
                    None,

                "aphia_id":
                    None,

                "batch_status":
                    "FAILED",

                "batch_error":
                    f"{type(error).__name__}: "
                    f"{error}",

                "output_directory":
                    None,
            })

            print(
                f"Batch status: FAILED — "
                f"{scientific_name}"
            )

            print(
                f"Reason: "
                f"{type(error).__name__}: "
                f"{error}"
            )

    summary = pd.DataFrame(
        summary_rows
    )

    summary_file = (
        output_directory
        / "batch_summary.csv"
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    successful = int(
        (
            summary["batch_status"]
            == "SUCCESS"
        ).sum()
    )

    failed = int(
        (
            summary["batch_status"]
            == "FAILED"
        ).sum()
    )

    print()
    print("=" * 70)
    print("BATCH COMPLETE")
    print("=" * 70)
    print(f"Submitted: {total}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(
        "Batch summary:",
        summary_file.resolve(),
    )
    print()

    return summary