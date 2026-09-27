"""
Test the WoRMS → GBIF component of TaxonWeave.

This script:

1. Takes a user-supplied taxonomic name.
2. Resolves it through WoRMS.
3. Determines the currently accepted name.
4. Matches the accepted taxon against GBIF.
5. Counts all GBIF occurrence records.
6. Downloads all available occurrence records.
7. Generates basic data summaries.
8. Saves the occurrence dataset as a CSV file.
"""

from taxonweave.worms import (
    resolve_name,
    get_accepted_record
)

from taxonweave.gbif import (
    resolve_gbif_taxon,
    get_occurrence_count,
    get_occurrences
)


# ==================================================
# USER QUERY
# ==================================================

species = "Mercierella enigmatica"


# ==================================================
# STEP 1: RESOLVE NAME THROUGH WORMS
# ==================================================

worms_record = resolve_name(species)

if worms_record is None:

    raise ValueError(
        "Species could not be resolved through WoRMS."
    )


accepted = get_accepted_record(worms_record)

accepted_name = accepted.get("scientificname")

accepted_aphia_id = accepted.get("AphiaID")


print("\nWORMS RESOLUTION\n")

print("User query:", species)

print(
    "WoRMS matched name:",
    worms_record.get("scientificname")
)

print(
    "WoRMS matched status:",
    worms_record.get("status")
)

print(
    "Accepted name:",
    accepted_name
)

print(
    "Accepted WoRMS AphiaID:",
    accepted_aphia_id
)


# ==================================================
# STEP 2: RESOLVE ACCEPTED NAME THROUGH GBIF
# ==================================================

gbif_taxon = resolve_gbif_taxon(
    accepted_name
)

taxon_key = gbif_taxon.get("usageKey")


if taxon_key is None:

    raise ValueError(
        "GBIF could not resolve the accepted taxon."
    )


print("\nGBIF TAXON MATCH\n")

print(
    "Scientific name:",
    gbif_taxon.get("scientificName")
)

print(
    "Canonical name:",
    gbif_taxon.get("canonicalName")
)

print(
    "Rank:",
    gbif_taxon.get("rank")
)

print(
    "Status:",
    gbif_taxon.get("status")
)

print(
    "GBIF taxon key:",
    taxon_key
)

print(
    "Match type:",
    gbif_taxon.get("matchType")
)

print(
    "Confidence:",
    gbif_taxon.get("confidence")
)


# ==================================================
# STEP 3: COUNT GBIF OCCURRENCES
# ==================================================

count = get_occurrence_count(
    taxon_key
)


print("\nGBIF OCCURRENCES\n")

print(
    "Total records reported by GBIF:",
    f"{count:,}"
)


# ==================================================
# STEP 4: DOWNLOAD ALL GBIF OCCURRENCES
# ==================================================

print("\nDOWNLOADING GBIF DATA\n")

occurrences = get_occurrences(
    taxon_key
)


print("\nGBIF DATASET RETRIEVED\n")

print(
    "Records downloaded:",
    f"{len(occurrences):,}"
)


# ==================================================
# STEP 5: DATA COMPLETENESS SUMMARY
# ==================================================

print("\nDATA SUMMARY\n")


records_with_coordinates = (
    occurrences["latitude"]
    .notna()
    .sum()
)


records_with_identification = (
    occurrences["identified_by"]
    .notna()
    .sum()
)


records_with_catalog_numbers = (
    occurrences["catalog_number"]
    .notna()
    .sum()
)


countries_represented = (
    occurrences["country"]
    .nunique()
)


print(
    "Records with coordinates:",
    f"{records_with_coordinates:,}"
)

print(
    "Records with identification information:",
    f"{records_with_identification:,}"
)

print(
    "Records with catalog numbers:",
    f"{records_with_catalog_numbers:,}"
)

print(
    "Countries represented:",
    f"{countries_represented:,}"
)


# ==================================================
# STEP 6: BASIS OF RECORD SUMMARY
# ==================================================

print("\nBASIS OF RECORD\n")


basis_summary = (
    occurrences["basis_of_record"]
    .value_counts(
        dropna=False
    )
)


print(
    basis_summary.to_string()
)


# ==================================================
# STEP 7: COUNTRY SUMMARY
# ==================================================

print("\nTOP 20 COUNTRIES\n")


country_summary = (
    occurrences["country"]
    .value_counts()
    .head(20)
)


print(
    country_summary.to_string()
)


# ==================================================
# STEP 8: IDENTIFICATION METADATA
# ==================================================

print("\nIDENTIFICATION METADATA\n")


if records_with_identification > 0:

    identifiers = (
        occurrences["identified_by"]
        .dropna()
        .value_counts()
        .head(20)
    )

    print("Top identifiers:\n")

    print(
        identifiers.to_string()
    )

else:

    print(
        "No identification metadata available."
    )


# ==================================================
# STEP 9: SAVE COMPLETE DATASET
# ==================================================

output_file = "gbif_test_occurrences.csv"


occurrences.to_csv(
    output_file,
    index=False
)


print("\nOUTPUT\n")

print(
    f"Saved {len(occurrences):,} records to:"
)

print(
    output_file
)


print("\nTaxonWeave GBIF test complete.\n")