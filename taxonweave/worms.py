"""
WoRMS interface for TaxonWeave.

Functions in this module resolve scientific names against
the World Register of Marine Species (WoRMS), retrieve their
taxonomic classification, retrieve known synonyms, and resolve
unaccepted names to their currently accepted taxon.
"""

import requests


WORMS_BASE_URL = "https://www.marinespecies.org/rest"


def resolve_name(scientific_name):
    """
    Resolve a scientific name against WoRMS.

    Parameters
    ----------
    scientific_name : str
        Scientific name supplied by the user.

    Returns
    -------
    dict or None
        WoRMS taxonomic record, or None if no match is found.
    """

    url = f"{WORMS_BASE_URL}/AphiaRecordsByName/{scientific_name}"

    params = {
        "like": "false",
        "marine_only": "false",
    }

    response = requests.get(
        url,
        params=params,
        timeout=30
    )

    response.raise_for_status()

    records = response.json()

    if not records:
        return None

    return records[0]


def get_record_by_aphia_id(aphia_id):
    """
    Retrieve a WoRMS taxonomic record using its AphiaID.
    """

    url = f"{WORMS_BASE_URL}/AphiaRecordByAphiaID/{aphia_id}"

    response = requests.get(
        url,
        timeout=30
    )

    response.raise_for_status()

    return response.json()


def get_accepted_record(record):
    """
    Resolve a WoRMS record to its currently accepted taxon.

    If the supplied record is already accepted, the original
    record is returned.

    If it is unaccepted, the valid_AphiaID supplied by WoRMS
    is used to retrieve the accepted record.
    """

    if record is None:
        return None

    if record.get("status") == "accepted":
        return record

    valid_aphia_id = record.get("valid_AphiaID")

    if valid_aphia_id is None:
        return record

    return get_record_by_aphia_id(valid_aphia_id)


def get_classification(aphia_id):
    """
    Retrieve the complete taxonomic classification for a WoRMS AphiaID.
    """

    url = f"{WORMS_BASE_URL}/AphiaClassificationByAphiaID/{aphia_id}"

    response = requests.get(
        url,
        timeout=30
    )

    response.raise_for_status()

    return response.json()


def flatten_classification(classification):
    """
    Convert the nested WoRMS classification into a simple list.

    Returns
    -------
    list of dict
        Each dictionary contains rank, scientific name, and AphiaID.
    """

    flattened = []

    current = classification

    while current is not None:

        flattened.append(
            {
                "rank": current.get("rank"),
                "scientificname": current.get("scientificname"),
                "AphiaID": current.get("AphiaID"),
            }
        )

        current = current.get("child")

    return flattened


def get_synonyms(aphia_id):
    """
    Retrieve synonyms associated with a WoRMS AphiaID.
    """

    url = f"{WORMS_BASE_URL}/AphiaSynonymsByAphiaID/{aphia_id}"

    response = requests.get(
        url,
        timeout=30
    )

    response.raise_for_status()

    records = response.json()

    if not records:
        return []

    return records