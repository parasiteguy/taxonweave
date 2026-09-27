"""
TaxonWeave WoRMS-GenBank integration test.

This test uses Alitta succinea as a taxonomically more complex
case than Ficopomatus enigmaticus.

TaxonWeave will:

1. Resolve Alitta succinea through WoRMS.
2. Retrieve its current accepted concept.
3. Retrieve classification and synonyms.
4. Search GenBank independently under the accepted name
   and every WoRMS synonym.
5. Combine duplicate NCBI records.
6. Preserve the search-name provenance of every record.
7. Download each unique nucleotide record once.
8. Summarize deposited identifications, markers, geography,
   collection metadata, and voucher metadata.
9. Export the complete evidence package.
"""

from taxonweave.genbank import (
    configure_entrez
)

from taxonweave.query import (
    query_species
)


# ============================================================
# NCBI CONFIGURATION
# ============================================================

# Replace with your real email address.

ENTREZ_EMAIL = "davinack_drew@wheatoncollege.edu"


configure_entrez(
    email=ENTREZ_EMAIL
)


# ============================================================
# TEST TAXON
# ============================================================

species = "Alitta succinea"


# ============================================================
# RUN QUERY
# ============================================================

report = query_species(
    species
)


# ============================================================
# DISPLAY REPORT
# ============================================================

report.summary()


# ============================================================
# EXPORT REPORT
# ============================================================

report.export(
    "Alitta_succinea_report"
)