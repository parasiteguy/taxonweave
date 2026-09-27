"""
Tests for the TaxonWeave command-line interface.
"""

from taxonweave.cli import build_parser


def test_cli_parser_builds():
    parser = build_parser()
    assert parser is not None


def test_query_command_parses_scientific_name():

    parser = build_parser()

    args = parser.parse_args(
        [
            "query",
            "Alitta succinea",
        ]
    )

    assert args.command == "query"
    assert (
        args.scientific_name
        == "Alitta succinea"
    )
    assert args.export is False
    assert args.output is None


def test_query_command_parses_export():

    parser = build_parser()

    args = parser.parse_args(
        [
            "query",
            "Alitta succinea",
            "--export",
        ]
    )

    assert args.command == "query"
    assert args.export is True


def test_query_command_parses_output_directory():

    parser = build_parser()

    args = parser.parse_args(
        [
            "query",
            "Alitta succinea",
            "--export",
            "--output",
            "alitta_results",
        ]
    )

    assert (
        args.output
        == "alitta_results"
    )


def test_batch_command_parses_input_file():

    parser = build_parser()

    args = parser.parse_args(
        [
            "batch",
            "species.csv",
        ]
    )

    assert args.command == "batch"
    assert (
        args.input_file
        == "species.csv"
    )

    assert (
        args.output
        == "taxonweave_batch"
    )


def test_batch_command_parses_output_directory():

    parser = build_parser()

    args = parser.parse_args(
        [
            "batch",
            "species.csv",
            "--output",
            "polychaete_results",
        ]
    )

    assert args.command == "batch"

    assert (
        args.output
        == "polychaete_results"
    )