"""
GBIF interface for TaxonWeave.

This module resolves scientific names against GBIF and retrieves
occurrence records.

The occurrence downloader is designed to tolerate temporary API
or network failures through automatic retries and exponential
backoff.
"""

import time

import pandas as pd
import requests


GBIF_BASE_URL = "https://api.gbif.org/v1"


# ============================================================
# GBIF TAXON RESOLUTION
# ============================================================

def resolve_gbif_taxon(scientific_name):
    """
    Resolve a scientific name against the GBIF species backbone.

    Parameters
    ----------
    scientific_name : str
        Scientific name to resolve.

    Returns
    -------
    dict
        GBIF species match.
    """

    url = f"{GBIF_BASE_URL}/species/match"

    params = {
        "name": scientific_name
    }

    response = requests.get(
        url,
        params=params,
        timeout=(10, 60)
    )

    response.raise_for_status()

    return response.json()


# ============================================================
# OCCURRENCE COUNT
# ============================================================

def get_occurrence_count(taxon_key):
    """
    Return the number of GBIF occurrence records associated
    with a GBIF taxon key.
    """

    url = f"{GBIF_BASE_URL}/occurrence/search"

    params = {
        "taxon_key": taxon_key,
        "limit": 0
    }

    response = requests.get(
        url,
        params=params,
        timeout=(10, 60)
    )

    response.raise_for_status()

    data = response.json()

    return data.get(
        "count",
        0
    )


# ============================================================
# SINGLE PAGE REQUEST WITH RETRY
# ============================================================

def request_occurrence_page(
    taxon_key,
    offset,
    limit=300,
    max_retries=5,
    initial_backoff=2
):
    """
    Retrieve one page of GBIF occurrence records.

    Temporary network failures and selected HTTP errors are
    retried automatically using exponential backoff.

    Parameters
    ----------
    taxon_key : int
        GBIF taxon key.

    offset : int
        Starting record offset.

    limit : int
        Number of records requested.

    max_retries : int
        Maximum number of attempts.

    initial_backoff : int or float
        Initial waiting period in seconds.

    Returns
    -------
    dict
        GBIF occurrence-search response.
    """

    url = f"{GBIF_BASE_URL}/occurrence/search"

    params = {
        "taxon_key": taxon_key,
        "limit": limit,
        "offset": offset
    }


    for attempt in range(
        1,
        max_retries + 1
    ):

        try:

            response = requests.get(
                url,
                params=params,

                # Separate connection and read timeouts.
                #
                # Connect timeout = 10 seconds
                # Read timeout = 90 seconds

                timeout=(10, 90)
            )


            # ------------------------------------------------
            # RETRY SELECTED SERVER/RATE-LIMIT RESPONSES
            # ------------------------------------------------

            if response.status_code in {
                429,
                500,
                502,
                503,
                504
            }:

                raise requests.exceptions.HTTPError(
                    f"Temporary GBIF HTTP "
                    f"{response.status_code}"
                )


            response.raise_for_status()

            return response.json()


        except (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
            requests.exceptions.HTTPError
        ) as error:

            # -----------------------------------------------
            # FINAL ATTEMPT FAILED
            # -----------------------------------------------

            if attempt == max_retries:

                raise RuntimeError(

                    "\nGBIF retrieval failed after "
                    f"{max_retries} attempts.\n"
                    f"Taxon key: {taxon_key}\n"
                    f"Offset: {offset:,}\n"
                    f"Original error: {error}"

                ) from error


            # -----------------------------------------------
            # EXPONENTIAL BACKOFF
            # -----------------------------------------------

            wait_seconds = (
                initial_backoff
                * (2 ** (attempt - 1))
            )


            print(
                f"GBIF request failed at "
                f"offset {offset:,} "
                f"(attempt {attempt}/{max_retries})."
            )


            print(
                f"Retrying in "
                f"{wait_seconds} seconds..."
            )


            time.sleep(
                wait_seconds
            )


# ============================================================
# OCCURRENCE PARSING
# ============================================================

def parse_occurrence(record):
    """
    Convert one GBIF occurrence record into the standardized
    TaxonWeave occurrence schema.
    """

    return {

        "gbif_id":
            record.get(
                "key"
            ),

        "scientific_name":
            record.get(
                "scientificName"
            ),

        "accepted_name":
            record.get(
                "acceptedScientificName"
            ),

        "basis_of_record":
            record.get(
                "basisOfRecord"
            ),

        "country":
            record.get(
                "country"
            ),

        "country_code":
            record.get(
                "countryCode"
            ),

        "state_province":
            record.get(
                "stateProvince"
            ),

        "locality":
            record.get(
                "locality"
            ),

        "latitude":
            record.get(
                "decimalLatitude"
            ),

        "longitude":
            record.get(
                "decimalLongitude"
            ),

        "event_date":
            record.get(
                "eventDate"
            ),

        "year":
            record.get(
                "year"
            ),

        "institution_code":
            record.get(
                "institutionCode"
            ),

        "collection_code":
            record.get(
                "collectionCode"
            ),

        "catalog_number":
            record.get(
                "catalogNumber"
            ),

        "identified_by":
            record.get(
                "identifiedBy"
            ),

        "recorded_by":
            record.get(
                "recordedBy"
            ),

        "dataset_name":
            record.get(
                "datasetTitle"
            )
    }


# ============================================================
# COMPLETE OCCURRENCE RETRIEVAL
# ============================================================

def get_occurrences(
    taxon_key,
    limit=300,
    max_retries=5,
    initial_backoff=2
):
    """
    Retrieve all GBIF occurrence records for a taxon.

    GBIF occurrence-search results are retrieved in pages.

    If a temporary network or server failure occurs, the failed
    page is retried automatically without discarding records
    already retrieved during the current run.

    Parameters
    ----------
    taxon_key : int
        GBIF taxon key.

    limit : int
        Number of records requested per page.

    max_retries : int
        Maximum attempts for an individual page.

    initial_backoff : int or float
        Initial retry delay.

    Returns
    -------
    pandas.DataFrame
        Standardized GBIF occurrence records.
    """

    records = []

    offset = 0


    while True:

        data = request_occurrence_page(

            taxon_key=taxon_key,

            offset=offset,

            limit=limit,

            max_retries=max_retries,

            initial_backoff=initial_backoff
        )


        results = data.get(
            "results",
            []
        )


        # ----------------------------------------------------
        # PARSE PAGE
        # ----------------------------------------------------

        for record in results:

            records.append(
                parse_occurrence(
                    record
                )
            )


        # ----------------------------------------------------
        # PROGRESS
        # ----------------------------------------------------

        print(
            f"Retrieved "
            f"{len(records):,} "
            f"GBIF records..."
        )


        # ----------------------------------------------------
        # DETERMINE WHETHER DOWNLOAD IS COMPLETE
        # ----------------------------------------------------

        end_of_records = data.get(
            "endOfRecords",
            False
        )


        if end_of_records:
            break


        if not results:
            break


        # Move to the next page using the actual number
        # returned rather than assuming GBIF returned `limit`.

        offset += len(
            results
        )


    return pd.DataFrame(
        records
    )