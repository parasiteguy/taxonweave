"""
WoRMS interface for TaxonWeave.

Functions in this module resolve scientific names against
the World Register of Marine Species (WoRMS), retrieve their
taxonomic classification, retrieve known synonyms, and resolve
unaccepted names to their currently accepted taxon.

WoRMS requests are handled through a shared request function that
provides resilience to transient API failures, including rate
limiting, server errors, timeouts, empty responses, and non-JSON
responses.

HTTP 204 (No Content) responses are treated as valid responses
indicating that WoRMS has no records to return for the requested
endpoint.
"""

import time
from urllib.parse import quote

import requests


WORMS_BASE_URL = "https://www.marinespecies.org/rest"

# Request behavior
WORMS_TIMEOUT = 30
WORMS_MAX_ATTEMPTS = 3
WORMS_RETRY_STATUS_CODES = {429, 500, 502, 503, 504}
WORMS_RETRY_DELAY = 1.0


class WoRMSError(RuntimeError):
    """
    Raised when TaxonWeave cannot obtain a usable response from WoRMS.
    """


def _worms_get_json(url, params=None):
    """
    Send a GET request to the WoRMS REST API and return decoded JSON.

    The request is retried when WoRMS returns a transient HTTP error,
    an empty response, invalid JSON, or a temporary connection problem.

    HTTP 204 (No Content) is treated as a valid response and returns
    None. This allows endpoints such as synonym retrieval to indicate
    legitimately that no records are available.

    Parameters
    ----------
    url : str
        Complete WoRMS REST API endpoint.

    params : dict, optional
        Query parameters supplied to the endpoint.

    Returns
    -------
    object or None
        JSON-decoded WoRMS response, or None when WoRMS returns
        HTTP 204 (No Content).

    Raises
    ------
    WoRMSError
        If a usable WoRMS response cannot be obtained after the
        configured number of attempts.
    """

    last_error = None

    for attempt in range(1, WORMS_MAX_ATTEMPTS + 1):

        try:
            response = requests.get(
                url,
                params=params,
                timeout=WORMS_TIMEOUT,
            )

            # HTTP 204 means that WoRMS successfully processed the
            # request but has no content to return. This is a valid
            # response for endpoints such as synonym retrieval when
            # no records are available.
            if response.status_code == 204:
                return None

            # Retry HTTP responses that are commonly transient.
            if response.status_code in WORMS_RETRY_STATUS_CODES:
                last_error = (
                    f"WoRMS returned HTTP {response.status_code}"
                )

                if attempt < WORMS_MAX_ATTEMPTS:
                    time.sleep(WORMS_RETRY_DELAY * attempt)
                    continue

                raise WoRMSError(
                    "TaxonWeave could not retrieve data from WoRMS "
                    f"after {WORMS_MAX_ATTEMPTS} attempts. "
                    f"The WoRMS service returned HTTP "
                    f"{response.status_code}. Please try again later."
                )

            # Convert other HTTP errors into a clearer
            # TaxonWeave/WoRMS error.
            try:
                response.raise_for_status()

            except requests.exceptions.HTTPError as exc:
                raise WoRMSError(
                    "TaxonWeave could not retrieve data from WoRMS. "
                    f"The WoRMS service returned HTTP "
                    f"{response.status_code}."
                ) from exc

            # A successful HTTP status does not necessarily guarantee
            # that WoRMS returned a usable response body.
            #
            # Note that HTTP 204 has already been handled above.
            # Therefore, a 200 response with an empty body is treated
            # as an unexpected/transient response rather than as a
            # legitimate absence of records.
            if not response.text or not response.text.strip():
                last_error = "WoRMS returned an empty response"

                if attempt < WORMS_MAX_ATTEMPTS:
                    time.sleep(WORMS_RETRY_DELAY * attempt)
                    continue

                raise WoRMSError(
                    "TaxonWeave could not retrieve data from WoRMS "
                    f"after {WORMS_MAX_ATTEMPTS} attempts. "
                    "The WoRMS service returned an unexpected empty "
                    "response. Please try again later."
                )

            # WoRMS normally returns JSON. Occasionally an upstream
            # service or proxy may return HTML or another non-JSON
            # response even with a successful HTTP status.
            try:
                return response.json()

            except requests.exceptions.JSONDecodeError as exc:
                last_error = "WoRMS returned a non-JSON response"

                if attempt < WORMS_MAX_ATTEMPTS:
                    time.sleep(WORMS_RETRY_DELAY * attempt)
                    continue

                raise WoRMSError(
                    "TaxonWeave could not retrieve data from WoRMS "
                    f"after {WORMS_MAX_ATTEMPTS} attempts. "
                    "The WoRMS service returned an invalid or non-JSON "
                    "response. Please try again later."
                ) from exc

        except (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
        ) as exc:

            last_error = str(exc)

            if attempt < WORMS_MAX_ATTEMPTS:
                time.sleep(WORMS_RETRY_DELAY * attempt)
                continue

            raise WoRMSError(
                "TaxonWeave could not connect to WoRMS "
                f"after {WORMS_MAX_ATTEMPTS} attempts. "
                "Please check your internet connection or try again later."
            ) from exc

        except requests.exceptions.RequestException as exc:
            raise WoRMSError(
                "TaxonWeave encountered an error while communicating "
                f"with WoRMS: {exc}"
            ) from exc

    # Defensive fallback. Execution should normally never reach here.
    raise WoRMSError(
        "TaxonWeave could not retrieve data from WoRMS. "
        f"Last error: {last_error}"
    )


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

    encoded_name = quote(scientific_name, safe="")

    url = f"{WORMS_BASE_URL}/AphiaRecordsByName/{encoded_name}"

    params = {
        "like": "false",
        "marine_only": "false",
    }

    records = _worms_get_json(
        url,
        params=params,
    )

    if not records:
        return None

    return records[0]


def get_record_by_aphia_id(aphia_id):
    """
    Retrieve a WoRMS taxonomic record using its AphiaID.
    """

    url = f"{WORMS_BASE_URL}/AphiaRecordByAphiaID/{aphia_id}"

    return _worms_get_json(url)


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

    return _worms_get_json(url)


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

    WoRMS may return HTTP 204 (No Content) when no synonyms are
    associated with the requested AphiaID. In that case this
    function returns an empty list.
    """

    url = f"{WORMS_BASE_URL}/AphiaSynonymsByAphiaID/{aphia_id}"

    records = _worms_get_json(url)

    if not records:
        return []

    return records