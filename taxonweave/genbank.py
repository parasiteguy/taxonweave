"""
GenBank interface for TaxonWeave.

This module searches NCBI Nucleotide using accepted and historical
taxonomic names, downloads each unique GenBank record once, and
extracts molecular and provenance metadata.

Design principles
-----------------
1. Preserve the taxonomic name represented in the GenBank record.
2. Preserve every search name through which a record was retrieved.
3. Standardize gene/marker terminology only when the annotation
   provides sufficient evidence.
4. Distinguish normalized, recognized, ambiguous, and unresolved
   molecular annotations.
5. Preserve raw annotations alongside standardized interpretations.
6. Treat specimen_voucher as metadata, not as proof that a verified
   museum voucher exists.
7. Leave WoRMS taxonomic reconciliation to the integrated query layer.
"""

from Bio import Entrez
from Bio import SeqIO

import pandas as pd
import re
import time


# ============================================================
# ENTREZ CONFIGURATION
# ============================================================

def configure_entrez(email, api_key=None):
    """
    Configure Biopython Entrez.
    """

    Entrez.email = email

    if api_key:
        Entrez.api_key = api_key


# ============================================================
# GENBANK SEARCH
# ============================================================

def search_genbank(scientific_name):
    """
    Search NCBI Nucleotide using an organism query.
    """

    print(
        f"Searching GenBank for: {scientific_name}"
    )

    term = f'"{scientific_name}"[Organism]'

    handle = Entrez.esearch(
        db="nucleotide",
        term=term,
        retmax=100000
    )

    result = Entrez.read(handle)

    handle.close()

    ids = result.get(
        "IdList",
        []
    )

    print(
        f"  Found {len(ids):,} nucleotide records."
    )

    return ids


def search_multiple_names(names):
    """
    Search GenBank independently for multiple taxonomic names.

    Results remain separate so search-name provenance is preserved.
    """

    results = {}

    for name in names:

        if not name:
            continue

        if name in results:
            continue

        results[name] = search_genbank(
            name
        )

        time.sleep(0.1)

    return results


# ============================================================
# COMBINE SEARCH RESULTS
# ============================================================

def combine_search_results(search_results):
    """
    Combine NCBI search results while retaining every search name
    through which each nucleotide record was retrieved.
    """

    unique_ids = []
    seen = set()

    search_name_map = {}

    for search_name, ids in search_results.items():

        for ncbi_id in ids:

            ncbi_id = str(
                ncbi_id
            )

            if ncbi_id not in search_name_map:

                search_name_map[
                    ncbi_id
                ] = []

            if (
                search_name
                not in
                search_name_map[ncbi_id]
            ):

                search_name_map[
                    ncbi_id
                ].append(
                    search_name
                )

            if ncbi_id not in seen:

                seen.add(
                    ncbi_id
                )

                unique_ids.append(
                    ncbi_id
                )

    return (
        unique_ids,
        search_name_map
    )


# ============================================================
# BASIC HELPERS
# ============================================================

def first_value(qualifiers, key):
    """
    Return the first value of a GenBank qualifier.
    """

    value = qualifiers.get(key)

    if not value:
        return None

    if isinstance(value, list):

        if not value:
            return None

        return value[0]

    return value


def join_values(qualifiers, key):
    """
    Join all values of a GenBank qualifier.
    """

    value = qualifiers.get(key)

    if not value:
        return None

    if isinstance(value, list):

        return "; ".join(
            str(item)
            for item in value
        )

    return str(value)


def clean_text(value):
    """
    Normalize whitespace without changing substantive content.
    """

    if value is None:
        return None

    value = re.sub(
        r"\s+",
        " ",
        str(value)
    ).strip()

    if not value:
        return None

    return value


# ============================================================
# FEATURE EXTRACTION
# ============================================================

def get_source_feature(record):
    """
    Return the source feature.
    """

    for feature in record.features:

        if feature.type == "source":
            return feature

    return None


def get_gene_features(record):
    """
    Return features potentially containing molecular annotations.
    """

    relevant_types = {
        "gene",
        "CDS",
        "rRNA",
        "tRNA",
        "misc_RNA",
        "misc_feature"
    }

    return [
        feature
        for feature in record.features
        if feature.type in relevant_types
    ]


def extract_feature_annotation(feature):
    """
    Extract the most informative annotation from one feature.
    """

    qualifiers = feature.qualifiers

    for key in [
        "gene",
        "product",
        "standard_name",
        "note"
    ]:

        value = first_value(
            qualifiers,
            key
        )

        if value:
            return clean_text(value)

    return None


def extract_raw_markers(record):
    """
    Extract unique molecular annotations while preserving order.
    """

    annotations = []

    for feature in get_gene_features(record):

        annotation = extract_feature_annotation(
            feature
        )

        if (
            annotation
            and
            annotation not in annotations
        ):

            annotations.append(
                annotation
            )

    if not annotations:

        description = clean_text(
            record.description
        )

        if description:
            annotations.append(
                description
            )

    return annotations


# ============================================================
# MOLECULAR ANNOTATION NORMALIZATION
# ============================================================

def matches_any(text, patterns):
    """
    Return True if text matches any supplied regular expression.
    """

    return any(
        re.search(pattern, text)
        for pattern in patterns
    )


def normalize_marker(annotation):
    """
    Interpret one molecular annotation.

    Returns
    -------
    normalized_name
    category
    status

    Statuses
    --------
    NORMALIZED
        Annotation supports a standardized marker designation.

    RECOGNIZED
        Annotation is biologically interpretable but is not forced
        into a traditional phylogenetic-marker abbreviation.

    AMBIGUOUS
        Annotation is meaningful but insufficiently specific for a
        defensible standardized assignment.

    UNRESOLVED
        Annotation cannot currently be interpreted reliably.
    """

    if annotation is None:

        return (
            None,
            "unknown",
            "UNRESOLVED"
        )

    raw = clean_text(annotation)
    text = raw.lower()


    # ========================================================
    # COI
    # ========================================================

    if matches_any(
        text,
        [
            r"\bcoi\b",
            r"\bco1\b",
            r"\bcox1\b",
            r"\bcox i\b",
            r"cytochrome c oxidase subunit i\b",
            r"cytochrome c oxidase subunit 1\b",
            r"cytochrome oxidase subunit i\b",
            r"cytochrome oxidase subunit 1\b"
        ]
    ):

        return (
            "COI",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # COII
    # ========================================================

    if matches_any(
        text,
        [
            r"\bcoii\b",
            r"\bco2\b",
            r"\bcox2\b",
            r"\bcox ii\b",
            r"cytochrome c oxidase subunit ii\b",
            r"cytochrome c oxidase subunit 2\b",
            r"cytochrome oxidase subunit ii\b",
            r"cytochrome oxidase subunit 2\b"
        ]
    ):

        return (
            "COII",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # COIII
    # ========================================================

    if matches_any(
        text,
        [
            r"\bcoiii\b",
            r"\bco3\b",
            r"\bcox3\b",
            r"\bcox iii\b",
            r"cytochrome c oxidase subunit iii\b",
            r"cytochrome c oxidase subunit 3\b",
            r"cytochrome oxidase subunit iii\b",
            r"cytochrome oxidase subunit 3\b"
        ]
    ):

        return (
            "COIII",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # CYTOCHROME B
    # ========================================================

    if matches_any(
        text,
        [
            r"\bcytb\b",
            r"\bcyt b\b",
	    r"\bcob\b",
            r"cytochrome b\b"
        ]
    ):

        return (
            "CytB",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # ATP6
    # ========================================================

    if matches_any(
        text,
        [
            r"\batp6\b",
            r"\batp 6\b",
            r"atp synthase subunit 6\b",
            r"atpase subunit 6\b"
        ]
    ):

        return (
            "ATP6",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # ATP8
    # ========================================================

    if matches_any(
        text,
        [
            r"\batp8\b",
            r"\batp 8\b",
            r"atp synthase subunit 8\b",
            r"atpase subunit 8\b"
        ]
    ):

        return (
            "ATP8",
            "mitochondrial_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # NADH DEHYDROGENASE GENES
    # ========================================================

    nd_patterns = {

        "ND4L": [
            r"\bnd4l\b",
            r"\bnad4l\b",
            r"nadh dehydrogenase subunit 4l\b"
        ],

        "ND1": [
            r"\bnd1\b",
            r"\bnad1\b",
            r"nadh dehydrogenase subunit 1\b"
        ],

        "ND2": [
            r"\bnd2\b",
            r"\bnad2\b",
            r"nadh dehydrogenase subunit 2\b"
        ],

        "ND3": [
            r"\bnd3\b",
            r"\bnad3\b",
            r"nadh dehydrogenase subunit 3\b"
        ],

        "ND4": [
            r"\bnd4\b",
            r"\bnad4\b",
            r"nadh dehydrogenase subunit 4\b"
        ],

        "ND5": [
            r"\bnd5\b",
            r"\bnad5\b",
            r"nadh dehydrogenase subunit 5\b"
        ],

        "ND6": [
            r"\bnd6\b",
            r"\bnad6\b",
            r"nadh dehydrogenase subunit 6\b"
        ]
    }

    for marker in [
        "ND4L",
        "ND1",
        "ND2",
        "ND3",
        "ND4",
        "ND5",
        "ND6"
    ]:

        if matches_any(
            text,
            nd_patterns[marker]
        ):

            return (
                marker,
                "mitochondrial_protein_coding",
                "NORMALIZED"
            )


    # ========================================================
    # MITOCHONDRIAL rRNA
    # ========================================================

    if matches_any(
        text,
        [
            r"\b12s\b",
            r"12s ribosomal rna\b",
            r"12s rrna\b"
        ]
    ):

        return (
            "12S",
            "mitochondrial_rRNA",
            "NORMALIZED"
        )


    if matches_any(
        text,
        [
            r"\b16s\b",
            r"16s ribosomal rna\b",
            r"16s rrna\b"
        ]
    ):

        return (
            "16S",
            "mitochondrial_rRNA",
            "NORMALIZED"
        )


    # ========================================================
    # NUCLEAR rRNA
    # ========================================================

    if matches_any(
        text,
        [
            r"\b18s\b",
            r"18s ribosomal rna\b",
            r"18s rrna\b"
        ]
    ):

        return (
            "18S",
            "nuclear_rRNA",
            "NORMALIZED"
        )


    if matches_any(
        text,
        [
            r"\b28s\b",
            r"28s ribosomal rna\b",
            r"28s rrna\b"
        ]
    ):

        return (
            "28S",
            "nuclear_rRNA",
            "NORMALIZED"
        )


    # ========================================================
    # HISTONE H3
    # ========================================================

    if matches_any(
        text,
        [
            r"\bh3\b",
            r"histone h3\b"
        ]
    ):

        return (
            "H3",
            "nuclear_protein_coding",
            "NORMALIZED"
        )


    # ========================================================
    # MITOCHONDRIAL tRNA
    # ========================================================

    # Complete mitochondrial genomes frequently contain gene
    # symbols such as trnA, trnN, trnG, etc. These are legitimate
    # molecular annotations and should not be labeled unresolved.

    if re.fullmatch(
        r"trn[a-z0-9]+(?:\([a-z]+\))?",
        text
    ):

        return (
            raw,
            "mitochondrial_tRNA",
            "RECOGNIZED"
        )


    if (
        "transfer rna"
        in text
        or
        "trna"
        in text
    ):

        return (
            raw,
            "transfer_RNA",
            "RECOGNIZED"
        )


    # ========================================================
    # AMBIGUOUS RIBOSOMAL ANNOTATIONS
    # ========================================================

    if (
        "large subunit ribosomal rna"
        in text
        or
        "large subunit rrna"
        in text
    ):

        return (
            None,
            "ribosomal_RNA",
            "AMBIGUOUS"
        )


    if (
        "small subunit ribosomal rna"
        in text
        or
        "small subunit rrna"
        in text
    ):

        return (
            None,
            "ribosomal_RNA",
            "AMBIGUOUS"
        )


    # ========================================================
    # OTHER RECOGNIZED GENES
    # ========================================================

    recognized_gene_terms = [

        "globin",
        "hemoglobin",
        "haemoglobin",
        "hemerythrin",
        "haemerythrin",
        "aquaporin",
        "glutathione s-transferase",
        "superoxide dismutase",
        "elongation factor",
        "actin",
        "tubulin",
        "heat shock protein",
        "hsp70",
        "hsp90"
    ]

    if any(
        term in text
        for term in recognized_gene_terms
    ):

        return (
            raw,
            "nuclear_or_nonmitochondrial_gene",
            "RECOGNIZED"
        )


    # ========================================================
    # OTHERWISE UNRESOLVED
    # ========================================================

    return (
        None,
        "unknown",
        "UNRESOLVED"
    )


# ============================================================
# RECORD-LEVEL MOLECULAR INTERPRETATION
# ============================================================

def interpret_record_markers(record):
    """
    Interpret all molecular annotations associated with a record.
    """

    raw_markers = extract_raw_markers(
        record
    )

    if not raw_markers:

        return {

            "marker_raw":
                None,

            "marker_normalized":
                None,

            "marker_category":
                "unknown",

            "normalization_status":
                "UNRESOLVED"
        }


    normalized_markers = []
    categories = []
    statuses = []


    for raw_marker in raw_markers:

        (
            normalized,
            category,
            status

        ) = normalize_marker(
            raw_marker
        )


        if (
            normalized
            and
            normalized not in normalized_markers
        ):

            normalized_markers.append(
                normalized
            )


        if (
            category
            and
            category not in categories
        ):

            categories.append(
                category
            )


        if (
            status
            and
            status not in statuses
        ):

            statuses.append(
                status
            )


    marker_raw = "; ".join(
        raw_markers
    )


    marker_normalized = (

        "; ".join(
            normalized_markers
        )

        if normalized_markers

        else None
    )


    marker_category = (

        "; ".join(
            categories
        )

        if categories

        else "unknown"
    )


    normalization_status = (

        "; ".join(
            statuses
        )

        if statuses

        else "UNRESOLVED"
    )


    return {

        "marker_raw":
            marker_raw,

        "marker_normalized":
            marker_normalized,

        "marker_category":
            marker_category,

        "normalization_status":
            normalization_status
    }


# ============================================================
# REFERENCE METADATA
# ============================================================

def extract_reference_metadata(record):
    """
    Extract basic metadata from the first GenBank reference.
    """

    references = record.annotations.get(
        "references"
    )

    if not references:

        return {

            "reference_title":
                None,

            "reference_journal":
                None,

            "reference_authors":
                None
        }


    reference = references[0]


    return {

        "reference_title":
            clean_text(
                getattr(
                    reference,
                    "title",
                    None
                )
            ),

        "reference_journal":
            clean_text(
                getattr(
                    reference,
                    "journal",
                    None
                )
            ),

        "reference_authors":
            clean_text(
                getattr(
                    reference,
                    "authors",
                    None
                )
            )
    }


# ============================================================
# GENBANK RECORD PARSING
# ============================================================

def parse_genbank_record(
    record,
    found_via_name=None,
    accepted_name=None
):
    """
    Convert one GenBank record into a standardized molecular
    evidence record.
    """

    source_feature = get_source_feature(
        record
    )

    qualifiers = (

        source_feature.qualifiers

        if source_feature

        else {}
    )


    # ========================================================
    # TAXONOMIC IDENTIFICATION
    # ========================================================

    deposited_name = first_value(
        qualifiers,
        "organism"
    )

    if deposited_name is None:

        deposited_name = (
            record.annotations.get(
                "organism"
            )
        )

    deposited_name = clean_text(
        deposited_name
    )


    # ========================================================
    # SOURCE METADATA
    # ========================================================

    geo_loc_name = first_value(
        qualifiers,
        "geo_loc_name"
    )

    country = first_value(
        qualifiers,
        "country"
    )

    locality = first_value(
        qualifiers,
        "locality"
    )

    lat_lon = first_value(
        qualifiers,
        "lat_lon"
    )

    collection_date = first_value(
        qualifiers,
        "collection_date"
    )

    collected_by = first_value(
        qualifiers,
        "collected_by"
    )

    identified_by = first_value(
        qualifiers,
        "identified_by"
    )

    specimen_voucher = first_value(
        qualifiers,
        "specimen_voucher"
    )

    isolate = first_value(
        qualifiers,
        "isolate"
    )

    strain = first_value(
        qualifiers,
        "strain"
    )

    clone = first_value(
        qualifiers,
        "clone"
    )

    isolation_source = first_value(
        qualifiers,
        "isolation_source"
    )

    source_note = join_values(
        qualifiers,
        "note"
    )

    db_xref = join_values(
        qualifiers,
        "db_xref"
    )


    # ========================================================
    # MOLECULAR EVIDENCE
    # ========================================================

    marker_info = interpret_record_markers(
        record
    )


    # ========================================================
    # REFERENCES
    # ========================================================

    reference_info = extract_reference_metadata(
        record
    )


    # ========================================================
    # SEQUENCE METADATA
    # ========================================================

    accession = record.id

    accession_base = (

        accession.split(".")[0]

        if accession

        else None
    )

    sequence_length = len(
        record.seq
    )

    molecule_type = (
        record.annotations.get(
            "molecule_type"
        )
    )

    topology = (
        record.annotations.get(
            "topology"
        )
    )


    # ========================================================
    # SEARCH PROVENANCE
    # ========================================================

    if isinstance(
        found_via_name,
        list
    ):

        found_via_name_text = "; ".join(
            found_via_name
        )

    else:

        found_via_name_text = (
            found_via_name
        )


    # ========================================================
    # VOUCHER METADATA
    # ========================================================

    voucher_present = bool(
        specimen_voucher
    )


    # ========================================================
    # RETURN
    # ========================================================

    return {

        "accession":
            accession,

        "accession_base":
            accession_base,

        "sequence_length":
            sequence_length,

        "molecule_type":
            molecule_type,

        "topology":
            topology,


        # Taxonomic provenance

        "deposited_name":
            deposited_name,

        "accepted_name":
            accepted_name,

        "found_via_name":
            found_via_name_text,


        # Molecular evidence

        "marker_raw":
            marker_info[
                "marker_raw"
            ],

        "marker_normalized":
            marker_info[
                "marker_normalized"
            ],

        "marker_category":
            marker_info[
                "marker_category"
            ],

        "normalization_status":
            marker_info[
                "normalization_status"
            ],


        # Geography and collection

        "geo_loc_name":
            clean_text(
                geo_loc_name
            ),

        "country":
            clean_text(
                country
            ),

        "locality":
            clean_text(
                locality
            ),

        "lat_lon":
            clean_text(
                lat_lon
            ),

        "collection_date":
            clean_text(
                collection_date
            ),

        "collected_by":
            clean_text(
                collected_by
            ),

        "identified_by":
            clean_text(
                identified_by
            ),


        # Specimen/sample provenance

        "voucher_raw":
            clean_text(
                specimen_voucher
            ),

        "voucher_present":
            voucher_present,

        "isolate":
            clean_text(
                isolate
            ),

        "strain":
            clean_text(
                strain
            ),

        "clone":
            clean_text(
                clone
            ),

        "isolation_source":
            clean_text(
                isolation_source
            ),

        "source_note":
            clean_text(
                source_note
            ),

        "db_xref":
            clean_text(
                db_xref
            ),


        # Publication provenance

        "reference_title":
            reference_info[
                "reference_title"
            ],

        "reference_journal":
            reference_info[
                "reference_journal"
            ],

        "reference_authors":
            reference_info[
                "reference_authors"
            ]
    }


# ============================================================
# ACCESSION / SEARCH-PROVENANCE MAPPING
# ============================================================

def build_accession_search_map(
    batch_ids,
    search_name_map
):
    """
    Associate accession versions with the taxonomic names through
    which their NCBI IDs were retrieved.
    """

    accession_search_map = {}

    if not batch_ids:
        return accession_search_map


    handle = Entrez.esummary(

        db="nucleotide",

        id=",".join(
            batch_ids
        ),

        retmode="xml"
    )

    summaries = Entrez.read(
        handle
    )

    handle.close()


    for index, summary in enumerate(
        summaries
    ):

        accession = summary.get(
            "AccessionVersion"
        )

        uid = (

            summary.get("Id")
            or
            summary.get("uid")
            or
            summary.get("Uid")
        )

        if uid is not None:

            uid = str(uid)


        if (
            uid is None
            and
            index < len(batch_ids)
        ):

            uid = str(
                batch_ids[index]
            )


        if accession and uid:

            accession_search_map[
                accession
            ] = search_name_map.get(
                uid,
                []
            )


    return accession_search_map


# ============================================================
# DOWNLOAD
# ============================================================

def fetch_genbank_records(
    ids,
    search_name_map=None,
    accepted_name=None,
    batch_size=100
):
    """
    Download and parse each unique GenBank record once.
    """

    if search_name_map is None:
        search_name_map = {}


    ids = [
        str(ncbi_id)
        for ncbi_id in ids
    ]

    total = len(ids)


    print(
        f"Downloading {total:,} unique "
        f"GenBank records..."
    )


    parsed_records = []


    for start in range(
        0,
        total,
        batch_size
    ):

        end = min(
            start + batch_size,
            total
        )

        batch_ids = ids[
            start:end
        ]


        print(
            f"  Downloading "
            f"{start + 1:,}-{end:,} "
            f"of {total:,}..."
        )


        accession_search_map = (
            build_accession_search_map(

                batch_ids,
                search_name_map
            )
        )


        handle = Entrez.efetch(

            db="nucleotide",

            id=",".join(
                batch_ids
            ),

            rettype="gb",

            retmode="text"
        )


        records = SeqIO.parse(
            handle,
            "genbank"
        )


        for record in records:

            found_via_name = (
                accession_search_map.get(
                    record.id,
                    []
                )
            )


            parsed_records.append(

                parse_genbank_record(

                    record=record,

                    found_via_name=(
                        found_via_name
                    ),

                    accepted_name=(
                        accepted_name
                    )
                )
            )


        handle.close()

        time.sleep(0.1)


    return pd.DataFrame(
        parsed_records
    )