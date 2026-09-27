"""
Test the GenBank component of TaxonWeave.

This script:

1. Resolves a taxonomic name through WoRMS.
2. Retrieves the currently accepted name.
3. Retrieves WoRMS synonyms.
4. Searches GenBank independently using the accepted name
   and all WoRMS synonyms.
5. Preserves search-name provenance.
6. Combines duplicate NCBI records.
7. Downloads each unique GenBank record once.
8. Extracts taxonomic, molecular, geographic, collection,
   voucher, and publication metadata.
9. Summarizes marker normalization status.
10. Exports the standardized GenBank evidence table.
"""

from taxonweave.worms import (
    resolve_name,
    get_accepted_record,
    get_synonyms
)

from taxonweave.genbank import (
    configure_entrez,
    search_multiple_names,
    combine_search_results,
    fetch_genbank_records
)


# ============================================================
# NCBI CONFIGURATION
# ============================================================

# Replace with your real email address.

ENTREZ_EMAIL = "YOUR_EMAIL_HERE"

configure_entrez(
    email=ENTREZ_EMAIL
)


# ============================================================
# USER QUERY
# ============================================================

species = "Mercierella enigmatica"


# ============================================================
# STEP 1: WORMS RESOLUTION
# ============================================================

worms_record = resolve_name(
    species
)

if worms_record is None:

    raise ValueError(
        "Species could not be resolved through WoRMS."
    )


accepted = get_accepted_record(
    worms_record
)

accepted_name = accepted.get(
    "scientificname"
)

accepted_aphia_id = accepted.get(
    "AphiaID"
)


print("\nWORMS RESOLUTION\n")

print(
    "User query:",
    species
)

print(
    "WoRMS matched name:",
    worms_record.get(
        "scientificname"
    )
)

print(
    "WoRMS matched status:",
    worms_record.get(
        "status"
    )
)

print(
    "Accepted name:",
    accepted_name
)

print(
    "Accepted AphiaID:",
    accepted_aphia_id
)


# ============================================================
# STEP 2: WORMS SYNONYMS
# ============================================================

synonym_records = get_synonyms(
    accepted_aphia_id
)


synonym_names = []

for synonym in synonym_records:

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


print("\nWORMS SYNONYMS\n")

if synonym_names:

    for name in synonym_names:

        print(
            "-",
            name
        )

else:

    print(
        "No synonyms found."
    )


# ============================================================
# STEP 3: BUILD GENBANK SEARCH NAMES
# ============================================================

search_names = [
    accepted_name
]


for synonym in synonym_names:

    if synonym not in search_names:

        search_names.append(
            synonym
        )


print("\nGENBANK SEARCH NAMES\n")

for name in search_names:

    print(
        "-",
        name
    )


# ============================================================
# STEP 4: SEARCH GENBANK
# ============================================================

print("\nSEARCHING GENBANK\n")


search_results = search_multiple_names(
    search_names
)


unique_ids, search_name_map = (
    combine_search_results(
        search_results
    )
)


print("\nGENBANK SEARCH SUMMARY\n")


for name, ids in search_results.items():

    print(
        f"{name}: "
        f"{len(ids):,}"
    )


print(
    "\nUnique nucleotide records:",
    f"{len(unique_ids):,}"
)


if not unique_ids:

    raise ValueError(
        "No GenBank nucleotide records were found."
    )


# ============================================================
# STEP 5: DOWNLOAD + PARSE UNIQUE RECORDS
# ============================================================

print("\nDOWNLOADING GENBANK RECORDS\n")


genbank_records = fetch_genbank_records(
    ids=unique_ids,
    search_name_map=search_name_map,
    accepted_name=accepted_name,
    batch_size=100
)


print("\nGENBANK DATASET RETRIEVED\n")


print(
    "Unique records downloaded:",
    f"{len(genbank_records):,}"
)


# ============================================================
# STEP 6: METADATA COMPLETENESS
# ============================================================

print("\nMETADATA COMPLETENESS\n")


fields = {

    "Raw marker annotation":
        "marker_raw",

    "Normalized marker":
        "marker_normalized",

    "Geographic location":
        "geo_loc_name",

    "Legacy country":
        "country",

    "Locality":
        "locality",

    "Coordinates":
        "lat_lon",

    "Collection date":
        "collection_date",

    "Collector":
        "collected_by",

    "Identifier":
        "identified_by",

    "Voucher metadata":
        "voucher_raw"
}


for label, column in fields.items():

    count = (
        genbank_records[
            column
        ]
        .notna()
        .sum()
    )

    percent = (
        count
        / len(genbank_records)
        * 100
    )

    print(
        f"{label}: "
        f"{count:,} "
        f"({percent:.1f}%)"
    )


# ============================================================
# STEP 7: DEPOSITED TAXON NAMES
# ============================================================

print("\nDEPOSITED TAXON NAMES\n")


print(
    genbank_records[
        "deposited_name"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# STEP 8: SEARCH PROVENANCE
# ============================================================

print("\nSEARCH PROVENANCE\n")


print(
    genbank_records[
        "found_via_name"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# STEP 9: NORMALIZED MARKERS
# ============================================================

print("\nNORMALIZED MARKERS\n")


normalized_markers = (
    genbank_records[
        "marker_normalized"
    ]
    .value_counts(
        dropna=False
    )
)


print(
    normalized_markers.to_string()
)


# ============================================================
# STEP 10: MOLECULAR ANNOTATION STATUS
# ============================================================

print("\nMOLECULAR ANNOTATION STATUS\n")


status_counts = {
    "NORMALIZED": 0,
    "RECOGNIZED": 0,
    "AMBIGUOUS": 0,
    "UNRESOLVED": 0
}


for value in (
    genbank_records[
        "normalization_status"
    ]
    .dropna()
):

    statuses = [
        item.strip()
        for item in str(value).split(";")
    ]

    for status in statuses:

        if status in status_counts:

            status_counts[
                status
            ] += 1


for status, count in status_counts.items():

    print(
        f"{status}: "
        f"{count:,}"
    )


# ============================================================
# STEP 11: ANNOTATIONS REQUIRING REVIEW
# ============================================================

print("\nANNOTATIONS REQUIRING REVIEW\n")


review_mask = (
    genbank_records[
        "normalization_status"
    ]
    .fillna("")
    .str.contains(
        "AMBIGUOUS|UNRESOLVED",
        regex=True
    )
)


review_records = (
    genbank_records[
        review_mask
    ]
)


if len(review_records) == 0:

    print(
        "No ambiguous or unresolved "
        "marker annotations."
    )

else:

    review_columns = [
        "accession",
        "deposited_name",
        "marker_raw",
        "marker_normalized",
        "marker_category",
        "normalization_status"
    ]

    print(
        review_records[
            review_columns
        ]
        .to_string(
            index=False
        )
    )


# ============================================================
# STEP 12: GEOGRAPHIC SUMMARY
# ============================================================

print("\nGEOGRAPHIC METADATA\n")


geo_summary = (
    genbank_records[
        "geo_loc_name"
    ]
    .value_counts(
        dropna=False
    )
    .head(30)
)


print(
    geo_summary.to_string()
)


# ============================================================
# STEP 13: SPECIMEN VOUCHER METADATA
# ============================================================

print("\nSPECIMEN VOUCHER METADATA\n")


voucher_records = (
    genbank_records[
        genbank_records[
            "voucher_present"
        ] == True
    ]
)


print(
    "Records containing specimen_voucher:",
    f"{len(voucher_records):,}"
)


if len(voucher_records) > 0:

    columns = [
        "accession",
        "deposited_name",
        "marker_normalized",
        "geo_loc_name",
        "voucher_raw"
    ]

    print()

    print(
        voucher_records[
            columns
        ]
        .head(40)
        .to_string(
            index=False
        )
    )


# ============================================================
# STEP 14: SAVE STANDARDIZED DATA
# ============================================================

output_file = (
    "genbank_test_records.csv"
)


genbank_records.to_csv(
    output_file,
    index=False
)


print("\nOUTPUT\n")


print(
    f"Saved "
    f"{len(genbank_records):,} "
    f"standardized GenBank records to:"
)

print(
    output_file
)


print(
    "\nTaxonWeave GenBank test complete.\n"
)