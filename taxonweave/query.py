"""
Integrated WoRMS-GenBank taxonomic reconciliation for TaxonWeave.

The central purpose of this module is to connect a current WoRMS
taxonomic concept with molecular records retrieved from GenBank
while preserving the distinction between:

1. the name entered by the user,
2. the current accepted WoRMS name,
3. WoRMS synonyms,
4. the taxonomic name represented in each GenBank record, and
5. the taxonomic names through which each GenBank record was found.

A WoRMS synonym relationship is treated as a nomenclatural
assertion. It is not interpreted as independent evidence that a
particular sequence is biologically conspecific with the current
accepted taxon.
"""

from pathlib import Path

import pandas as pd

from taxonweave.worms import (
    resolve_name,
    get_accepted_record,
    get_classification,
    flatten_classification,
    get_synonyms
)

from taxonweave.genbank import (
    search_multiple_names,
    combine_search_results,
    fetch_genbank_records
)


# ============================================================
# NAME NORMALIZATION
# ============================================================

def normalize_taxon_name(name):
    """
    Normalize a taxon name for conservative string comparison.

    This does not perform taxonomic inference.
    """

    if name is None:
        return None

    return " ".join(
        str(name)
        .strip()
        .lower()
        .split()
    )


# ============================================================
# TAXONOMIC RECONCILIATION
# ============================================================

def reconcile_genbank_names(
    genbank,
    accepted_name,
    synonyms
):
    """
    Compare GenBank deposited names against the current WoRMS
    accepted name and WoRMS synonym set.

    Adds
    ----
    name_relationship

    reconciliation_status

    Categories
    ----------
    ACCEPTED_NAME
        GenBank deposited name exactly matches the current accepted
        WoRMS scientific name.

    WORMS_SYNONYM
        GenBank deposited name exactly matches a name currently
        returned by WoRMS as a synonym of the accepted taxon.

    UNRESOLVED_NAME
        GenBank contains a deposited organism name, but that name
        does not exactly match the accepted name or supplied WoRMS
        synonym set.

    MISSING_NAME
        No deposited organism name was available.

    Important
    ---------
    WORMS_SYNONYM means only that WoRMS currently asserts the
    nomenclatural relationship. It does not independently validate
    the biological identity of the sequence.
    """

    if genbank.empty:

        genbank[
            "name_relationship"
        ] = pd.Series(dtype="object")

        genbank[
            "reconciliation_status"
        ] = pd.Series(dtype="object")

        return genbank


    accepted_normalized = (
        normalize_taxon_name(
            accepted_name
        )
    )


    synonym_lookup = {}


    for synonym in synonyms:

        synonym_name = synonym.get(
            "scientificname"
        )

        normalized = normalize_taxon_name(
            synonym_name
        )

        if normalized:

            synonym_lookup[
                normalized
            ] = synonym_name


    relationships = []
    statuses = []


    for _, row in genbank.iterrows():

        deposited_name = row.get(
            "deposited_name"
        )

        deposited_normalized = (
            normalize_taxon_name(
                deposited_name
            )
        )


        # ----------------------------------------------------
        # MISSING
        # ----------------------------------------------------

        if not deposited_normalized:

            relationships.append(
                "MISSING_NAME"
            )

            statuses.append(
                "UNRESOLVED"
            )

            continue


        # ----------------------------------------------------
        # ACCEPTED NAME
        # ----------------------------------------------------

        if (
            deposited_normalized
            ==
            accepted_normalized
        ):

            relationships.append(
                "ACCEPTED_NAME"
            )

            statuses.append(
                "RECONCILED"
            )

            continue


        # ----------------------------------------------------
        # WORMS SYNONYM
        # ----------------------------------------------------

        if (
            deposited_normalized
            in synonym_lookup
        ):

            relationships.append(
                "WORMS_SYNONYM"
            )

            statuses.append(
                "RECONCILED_NOMENCLATURALLY"
            )

            continue


        # ----------------------------------------------------
        # OTHERWISE UNRESOLVED
        # ----------------------------------------------------

        relationships.append(
            "UNRESOLVED_NAME"
        )

        statuses.append(
            "UNRESOLVED"
        )


    genbank = genbank.copy()


    genbank[
        "name_relationship"
    ] = relationships


    genbank[
        "reconciliation_status"
    ] = statuses


    return genbank


# ============================================================
# REPORT OBJECT
# ============================================================

class TaxonWeaveReport:
    """
    Integrated WoRMS-GenBank reconciliation report.
    """

    def __init__(
        self,
        query_name,
        query_record,
        accepted_record,
        taxonomy,
        synonyms,
        search_results,
        genbank
    ):

        self.query_name = query_name

        self.query_record = query_record

        self.accepted_record = (
            accepted_record
        )

        self.accepted_name = (

            accepted_record.get(
                "scientificname"
            )

            if accepted_record

            else None
        )

        self.aphia_id = (

            accepted_record.get(
                "AphiaID"
            )

            if accepted_record

            else None
        )

        self.taxonomy = taxonomy

        self.synonyms = synonyms

        self.search_results = (
            search_results
        )

        self.genbank = genbank


    # ========================================================
    # SUMMARY DICTIONARY
    # ========================================================

    def summary_dict(self):

        total_records = len(
            self.genbank
        )


        names_with_hits = sum(

            1

            for ids
            in self.search_results.values()

            if ids
        )


        # ----------------------------------------------------
        # TAXONOMIC RECONCILIATION COUNTS
        # ----------------------------------------------------

        reconciliation_counts = {}


        if (
            not self.genbank.empty
            and
            "name_relationship"
            in self.genbank.columns
        ):

            reconciliation_counts = (

                self.genbank[
                    "name_relationship"
                ]
                .value_counts(
                    dropna=False
                )
                .to_dict()
            )


        # ----------------------------------------------------
        # MARKER STATUS COUNTS
        # ----------------------------------------------------

        marker_status_counts = {}


        if (
            not self.genbank.empty
            and
            "normalization_status"
            in self.genbank.columns
        ):

            for value in (
                self.genbank[
                    "normalization_status"
                ]
                .dropna()
            ):

                statuses = [

                    item.strip()

                    for item
                    in str(value).split(";")

                    if item.strip()
                ]


                for status in statuses:

                    marker_status_counts[
                        status
                    ] = (

                        marker_status_counts.get(
                            status,
                            0
                        )
                        + 1
                    )


        # ----------------------------------------------------
        # METADATA
        # ----------------------------------------------------

        def count_present(column):

            if (
                self.genbank.empty
                or
                column
                not in self.genbank.columns
            ):

                return 0

            return int(
                self.genbank[
                    column
                ]
                .notna()
                .sum()
            )


        voucher_records = 0


        if (
            not self.genbank.empty
            and
            "voucher_present"
            in self.genbank.columns
        ):

            voucher_records = int(

                self.genbank[
                    "voucher_present"
                ]
                .fillna(False)
                .astype(bool)
                .sum()
            )


        deposited_names = 0


        if (
            not self.genbank.empty
            and
            "deposited_name"
            in self.genbank.columns
        ):

            deposited_names = (

                self.genbank[
                    "deposited_name"
                ]
                .dropna()
                .nunique()
            )


        return {

            "query_name":
                self.query_name,

            "accepted_name":
                self.accepted_name,

            "aphia_id":
                self.aphia_id,

            "worms_synonyms":
                len(
                    self.synonyms
                ),

            "genbank_search_names":
                len(
                    self.search_results
                ),

            "search_names_with_hits":
                names_with_hits,

            "genbank_records":
                total_records,

            "deposited_taxon_names":
                deposited_names,

            "accepted_name_records":
                reconciliation_counts.get(
                    "ACCEPTED_NAME",
                    0
                ),

            "worms_synonym_records":
                reconciliation_counts.get(
                    "WORMS_SYNONYM",
                    0
                ),

            "unresolved_name_records":
                reconciliation_counts.get(
                    "UNRESOLVED_NAME",
                    0
                ),

            "missing_name_records":
                reconciliation_counts.get(
                    "MISSING_NAME",
                    0
                ),

            "marker_normalized_records":
                marker_status_counts.get(
                    "NORMALIZED",
                    0
                ),

            "marker_recognized_records":
                marker_status_counts.get(
                    "RECOGNIZED",
                    0
                ),

            "marker_ambiguous_records":
                marker_status_counts.get(
                    "AMBIGUOUS",
                    0
                ),

            "marker_unresolved_records":
                marker_status_counts.get(
                    "UNRESOLVED",
                    0
                ),

            "geographic_records":
                count_present(
                    "geo_loc_name"
                ),

            "coordinate_records":
                count_present(
                    "lat_lon"
                ),

            "collection_date_records":
                count_present(
                    "collection_date"
                ),

            "voucher_records":
                voucher_records
        }


    # ========================================================
    # HUMAN-READABLE SUMMARY
    # ========================================================

    def summary(self):

        stats = self.summary_dict()


        print("\n")
        print("=" * 70)
        print("TAXONWEAVE TAXONOMIC-MOLECULAR RECONCILIATION REPORT")
        print("=" * 70)


        # ====================================================
        # QUERY
        # ====================================================

        print("\nQUERY\n")

        print(
            "Name entered:",
            self.query_name
        )


        if self.query_record:

            print(
                "WoRMS matched name:",
                self.query_record.get(
                    "scientificname"
                )
            )

            print(
                "WoRMS matched status:",
                self.query_record.get(
                    "status"
                )
            )

            print(
                "Matched AphiaID:",
                self.query_record.get(
                    "AphiaID"
                )
            )


        # ====================================================
        # CURRENT CONCEPT
        # ====================================================

        print(
            "\nCURRENT WORMS CONCEPT\n"
        )

        print(
            "Accepted name:",
            self.accepted_name
        )

        print(
            "Authority:",
            self.accepted_record.get(
                "authority"
            )
        )

        print(
            "Accepted AphiaID:",
            self.aphia_id
        )


        # ====================================================
        # CLASSIFICATION
        # ====================================================

        print("\nCLASSIFICATION\n")

        for taxon in self.taxonomy:

            rank = taxon.get(
                "rank"
            )

            name = taxon.get(
                "scientificname"
            )

            if rank and name:

                print(
                    f"{rank}: {name}"
                )


        # ====================================================
        # SYNONYMS
        # ====================================================

        print("\nWORMS SYNONYMS\n")


        if not self.synonyms:

            print(
                "No WoRMS synonyms returned."
            )


        for synonym in self.synonyms:

            name = synonym.get(
                "scientificname"
            )

            authority = synonym.get(
                "authority"
            )

            aphia_id = synonym.get(
                "AphiaID"
            )


            text = f"- {name}"


            if authority:

                text += (
                    f" {authority}"
                )


            if aphia_id:

                text += (
                    f" [AphiaID {aphia_id}]"
                )


            print(text)


        # ====================================================
        # SEARCH PROVENANCE
        # ====================================================

        print(
            "\nGENBANK SEARCH PROVENANCE\n"
        )


        for name, ids in self.search_results.items():

            print(
                f"{name}: "
                f"{len(ids):,} records"
            )


        print(
            "\nSearch names producing records:",
            f"{stats['search_names_with_hits']}"
            f"/{stats['genbank_search_names']}"
        )


        # ====================================================
        # MOLECULAR EVIDENCE
        # ====================================================

        print(
            "\nGENBANK MOLECULAR EVIDENCE\n"
        )

        print(
            "Unique nucleotide records:",
            f"{stats['genbank_records']:,}"
        )

        print(
            "Distinct deposited taxon names:",
            f"{stats['deposited_taxon_names']:,}"
        )


        # ====================================================
        # DEPOSITED NAMES
        # ====================================================

        print(
            "\nDEPOSITED TAXON IDENTIFICATIONS\n"
        )


        if (
            not self.genbank.empty
            and
            "deposited_name"
            in self.genbank.columns
        ):

            print(

                self.genbank[
                    "deposited_name"
                ]
                .value_counts(
                    dropna=False
                )
                .to_string()
            )


        # ====================================================
        # TAXONOMIC RECONCILIATION
        # ====================================================

        print(
            "\nTAXONOMIC RECONCILIATION\n"
        )

        print(
            "Accepted-name records:",
            f"{stats['accepted_name_records']:,}"
        )

        print(
            "WoRMS-synonym records:",
            f"{stats['worms_synonym_records']:,}"
        )

        print(
            "Unresolved deposited names:",
            f"{stats['unresolved_name_records']:,}"
        )

        print(
            "Missing deposited names:",
            f"{stats['missing_name_records']:,}"
        )


        if (
            not self.genbank.empty
            and
            "name_relationship"
            in self.genbank.columns
        ):

            unresolved = self.genbank[

                self.genbank[
                    "name_relationship"
                ].isin(
                    [
                        "UNRESOLVED_NAME",
                        "MISSING_NAME"
                    ]
                )
            ]


            if not unresolved.empty:

                print(
                    "\nRecords requiring taxonomic review:"
                )

                columns = [

                    column

                    for column in [
                        "accession",
                        "deposited_name",
                        "found_via_name",
                        "name_relationship"
                    ]

                    if column
                    in unresolved.columns
                ]


                print()

                print(

                    unresolved[
                        columns
                    ]
                    .head(30)
                    .to_string(
                        index=False
                    )
                )


        # ====================================================
        # MARKERS
        # ====================================================

        print(
            "\nMOLECULAR MARKERS\n"
        )


        if (
            not self.genbank.empty
            and
            "marker_normalized"
            in self.genbank.columns
        ):

            print(

                self.genbank[
                    "marker_normalized"
                ]
                .value_counts(
                    dropna=False
                )
                .to_string()
            )


        # ====================================================
        # NORMALIZATION STATUS
        # ====================================================

        print(
            "\nMOLECULAR ANNOTATION STATUS\n"
        )

        print(
            "Records containing normalized annotations:",
            f"{stats['marker_normalized_records']:,}"
        )

        print(
            "Records containing recognized annotations:",
            f"{stats['marker_recognized_records']:,}"
        )

        print(
            "Records containing ambiguous annotations:",
            f"{stats['marker_ambiguous_records']:,}"
        )

        print(
            "Records containing unresolved annotations:",
            f"{stats['marker_unresolved_records']:,}"
        )


        # ====================================================
        # AMBIGUOUS / UNRESOLVED ANNOTATIONS
        # ====================================================

        print(
            "\nANNOTATIONS REQUIRING REVIEW\n"
        )


        if (
            not self.genbank.empty
            and
            "normalization_status"
            in self.genbank.columns
        ):

            review_mask = (

                self.genbank[
                    "normalization_status"
                ]
                .fillna("")
                .str.contains(
                    "AMBIGUOUS|UNRESOLVED",
                    regex=True
                )
            )


            review_records = (
                self.genbank[
                    review_mask
                ]
            )


            if review_records.empty:

                print(
                    "No ambiguous or unresolved annotations."
                )

            else:

                columns = [

                    column

                    for column in [
                        "accession",
                        "marker_raw",
                        "marker_normalized",
                        "marker_category",
                        "normalization_status"
                    ]

                    if column
                    in review_records.columns
                ]


                print(

                    review_records[
                        columns
                    ]
                    .head(40)
                    .to_string(
                        index=False
                    )
                )


        # ====================================================
        # METADATA
        # ====================================================

        print(
            "\nGENBANK METADATA COMPLETENESS\n"
        )

        print(
            "Geographic metadata:",
            f"{stats['geographic_records']:,}"
        )

        print(
            "Coordinates:",
            f"{stats['coordinate_records']:,}"
        )

        print(
            "Collection dates:",
            f"{stats['collection_date_records']:,}"
        )

        print(
            "specimen_voucher metadata:",
            f"{stats['voucher_records']:,}"
        )


        # ====================================================
        # VOUCHERS
        # ====================================================

        print(
            "\nSPECIMEN VOUCHER METADATA\n"
        )


        if (
            not self.genbank.empty
            and
            "voucher_present"
            in self.genbank.columns
        ):

            vouchers = self.genbank[

                self.genbank[
                    "voucher_present"
                ]
                .fillna(False)
                .astype(bool)
            ]


            print(
                "Records containing specimen_voucher:",
                f"{len(vouchers):,}"
            )


            if not vouchers.empty:

                columns = [

                    column

                    for column in [
                        "accession",
                        "deposited_name",
                        "name_relationship",
                        "marker_normalized",
                        "geo_loc_name",
                        "voucher_raw"
                    ]

                    if column
                    in vouchers.columns
                ]


                print()

                print(

                    vouchers[
                        columns
                    ]
                    .head(30)
                    .to_string(
                        index=False
                    )
                )


        print("\n")
        print("=" * 70)
        print("END TAXONWEAVE REPORT")
        print("=" * 70)
        print("\n")


    # ========================================================
    # EXPORT
    # ========================================================

    def export(
        self,
        directory_name=None
    ):

        if directory_name is None:

            safe_name = (
                self.accepted_name
                .replace(
                    " ",
                    "_"
                )
            )

            directory_name = (
                f"{safe_name}_report"
            )


        output_directory = Path(
            directory_name
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True
        )


        # ----------------------------------------------------
        # TAXONOMY
        # ----------------------------------------------------

        pd.DataFrame(
            self.taxonomy
        ).to_csv(

            output_directory
            / "taxonomy.csv",

            index=False
        )


        # ----------------------------------------------------
        # SYNONYMS
        # ----------------------------------------------------

        pd.DataFrame(
            self.synonyms
        ).to_csv(

            output_directory
            / "synonyms.csv",

            index=False
        )


        # ----------------------------------------------------
        # SEARCH PROVENANCE
        # ----------------------------------------------------

        search_rows = []


        for search_name, ids in self.search_results.items():

            search_rows.append({

                "search_name":
                    search_name,

                "records_retrieved":
                    len(ids),

                "retrieved_records":
                    bool(ids)
            })


        pd.DataFrame(
            search_rows
        ).to_csv(

            output_directory
            / "genbank_searches.csv",

            index=False
        )


        # ----------------------------------------------------
        # SEQUENCES / RECONCILIATION
        # ----------------------------------------------------

        self.genbank.to_csv(

            output_directory
            / "sequences.csv",

            index=False
        )


        # ----------------------------------------------------
        # SUMMARY
        # ----------------------------------------------------

        pd.DataFrame(
            [
                self.summary_dict()
            ]
        ).to_csv(

            output_directory
            / "summary.csv",

            index=False
        )


        print(
            "\nTaxonWeave report exported to:"
        )

        print(
            output_directory.resolve()
        )


# ============================================================
# MAIN QUERY
# ============================================================

def query_species(scientific_name):
    """
    Run the complete WoRMS-GenBank reconciliation workflow.
    """

    print("\n")
    print("=" * 65)
    print("TAXONWEAVE")
    print("=" * 65)

    print(
        "\nQuery:",
        scientific_name
    )


    # ========================================================
    # 1. WORMS
    # ========================================================

    print(
        "\n[1/2] Resolving taxonomic concept through WoRMS..."
    )


    query_record = resolve_name(
        scientific_name
    )


    if query_record is None:

        raise ValueError(
            f"WoRMS could not resolve "
            f"'{scientific_name}'."
        )


    accepted_record = get_accepted_record(
        query_record
    )


    if accepted_record is None:

        raise ValueError(
            "An accepted WoRMS concept "
            "could not be resolved."
        )


    accepted_name = accepted_record.get(
        "scientificname"
    )

    aphia_id = accepted_record.get(
        "AphiaID"
    )


    print(
        "      Query matched:",
        query_record.get(
            "scientificname"
        )
    )

    print(
        "      Query status:",
        query_record.get(
            "status"
        )
    )

    print(
        "      Accepted taxon:",
        accepted_name
    )


    classification_raw = (
        get_classification(
            aphia_id
        )
    )

    taxonomy = flatten_classification(
        classification_raw
    )


    synonyms = get_synonyms(
        aphia_id
    )


    synonym_names = []


    for synonym in synonyms:

        name = synonym.get(
            "scientificname"
        )

        if (
            name
            and
            name not in synonym_names
        ):

            synonym_names.append(
                name
            )


    print(
        "      WoRMS synonyms:",
        len(synonym_names)
    )


    # ========================================================
    # 2. GENBANK
    # ========================================================

    print(
        "\n[2/2] Reconciling with GenBank..."
    )


    search_names = [
        accepted_name
    ]


    for synonym_name in synonym_names:

        if (
            synonym_name
            not in search_names
        ):

            search_names.append(
                synonym_name
            )


    print(
        "      Taxonomic names to search:",
        len(search_names)
    )

    print()


    search_results = search_multiple_names(
        search_names
    )


    (
        unique_ids,
        search_name_map

    ) = combine_search_results(
        search_results
    )


    print(
        "\n      Unique nucleotide records:",
        f"{len(unique_ids):,}"
    )


    if unique_ids:

        genbank = fetch_genbank_records(

            ids=unique_ids,

            search_name_map=(
                search_name_map
            ),

            accepted_name=(
                accepted_name
            ),

            batch_size=100
        )

    else:

        genbank = pd.DataFrame()


    print(
        "      GenBank records downloaded:",
        f"{len(genbank):,}"
    )


    # ========================================================
    # TAXONOMIC RECONCILIATION
    # ========================================================

    print(
        "      Reconciling deposited names "
        "against WoRMS..."
    )


    genbank = reconcile_genbank_names(

        genbank=genbank,

        accepted_name=(
            accepted_name
        ),

        synonyms=(
            synonyms
        )
    )


    print(
        "      Taxonomic reconciliation complete."
    )


    # ========================================================
    # REPORT
    # ========================================================

    report = TaxonWeaveReport(

        query_name=(
            scientific_name
        ),

        query_record=(
            query_record
        ),

        accepted_record=(
            accepted_record
        ),

        taxonomy=(
            taxonomy
        ),

        synonyms=(
            synonyms
        ),

        search_results=(
            search_results
        ),

        genbank=(
            genbank
        )
    )


    print(
        "\nTaxonomic-molecular reconciliation complete."
    )


    return report