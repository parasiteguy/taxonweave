"""
Command-line interface for TaxonWeave.
"""

import argparse

from taxonweave.batch import run_batch
from taxonweave.query import query_species


def build_parser():
    """
    Construct the TaxonWeave command-line parser.
    """

    parser = argparse.ArgumentParser(
        prog="taxonweave",
        description=(
            "Provenance-aware reconciliation of taxonomic "
            "concepts and molecular sequence records."
        ),
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # --------------------------------------------------------
    # QUERY COMMAND
    # --------------------------------------------------------

    query_parser = subparsers.add_parser(
        "query",
        help=(
            "Reconcile a taxonomic concept "
            "with GenBank records."
        ),
    )

    query_parser.add_argument(
        "scientific_name",
        help='Scientific name, e.g. "Alitta succinea".',
    )

    query_parser.add_argument(
        "--export",
        action="store_true",
        help="Export the TaxonWeave report as CSV files.",
    )

    query_parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional output directory. "
            "Used only with --export."
        ),
    )

    # --------------------------------------------------------
    # BATCH COMMAND
    # --------------------------------------------------------

    batch_parser = subparsers.add_parser(
        "batch",
        help=(
            "Run TaxonWeave for multiple scientific "
            "names from a CSV file."
        ),
    )

    batch_parser.add_argument(
        "input_file",
        help=(
            "CSV file containing a "
            "'scientific_name' column."
        ),
    )

    batch_parser.add_argument(
        "--output",
        default="taxonweave_batch",
        help=(
            "Directory for batch results. "
            "Default: taxonweave_batch"
        ),
    )

    return parser


def main():
    """
    Run the TaxonWeave command-line interface.
    """

    parser = build_parser()

    args = parser.parse_args()

    if args.command == "query":

        report = query_species(
            args.scientific_name
        )

        report.summary()

        if args.export:

            report.export(
                directory_name=args.output
            )

    elif args.command == "batch":

        run_batch(
            input_file=args.input_file,
            output_directory=args.output,
        )


if __name__ == "__main__":
    main()