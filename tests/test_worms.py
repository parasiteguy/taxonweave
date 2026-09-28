"""
Tests for the TaxonWeave WoRMS interface.

These tests mock HTTP responses so that the test suite does not depend
on the availability or behavior of the live WoRMS REST API.

The tests cover:

1. Normal HTTP 200 responses containing valid JSON.
2. HTTP 204 No Content responses.
3. WoRMS synonym retrieval when no synonyms are returned.
4. Retry behavior for HTTP 429 rate limiting.
5. Retry behavior for transient server errors.
6. Empty HTTP 200 responses.
7. Invalid/non-JSON HTTP 200 responses.
8. Network timeouts.
9. Connection errors.
"""

from unittest.mock import Mock, patch

import pytest
import requests

from taxonweave.worms import (
    WoRMSError,
    _worms_get_json,
    get_synonyms,
    resolve_name,
)


# ---------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------


def make_response(status_code=200, json_data=None, text=None):
    """
    Construct a mocked requests.Response-like object.
    """

    response = Mock()
    response.status_code = status_code

    if text is None:
        if status_code == 204:
            response.text = ""
        elif json_data is not None:
            response.text = "mock-json-content"
        else:
            response.text = ""
    else:
        response.text = text

    if status_code >= 400:
        response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            f"HTTP {status_code}"
        )
    else:
        response.raise_for_status.return_value = None

    if json_data is not None:
        response.json.return_value = json_data

    return response


# ---------------------------------------------------------------------
# Normal successful responses
# ---------------------------------------------------------------------


@patch("taxonweave.worms.requests.get")
def test_worms_get_json_valid_response(mock_get):
    """
    A normal HTTP 200 response containing valid JSON should be returned.
    """

    expected = {
        "AphiaID": 153847,
        "scientificname": "Polydora websteri",
        "status": "accepted",
    }

    mock_get.return_value = make_response(
        status_code=200,
        json_data=expected,
    )

    result = _worms_get_json("https://example.org/test")

    assert result == expected
    assert mock_get.call_count == 1


@patch("taxonweave.worms.requests.get")
def test_resolve_name_valid_response(mock_get):
    """
    resolve_name() should return the first matching WoRMS record.
    """

    record = {
        "AphiaID": 153847,
        "scientificname": "Polydora websteri",
        "status": "accepted",
        "valid_AphiaID": 153847,
        "valid_name": "Polydora websteri",
    }

    mock_get.return_value = make_response(
        status_code=200,
        json_data=[record],
    )

    result = resolve_name("Polydora websteri")

    assert result == record
    assert result["AphiaID"] == 153847
    assert result["scientificname"] == "Polydora websteri"


# ---------------------------------------------------------------------
# HTTP 204 No Content
# ---------------------------------------------------------------------


@patch("taxonweave.worms.requests.get")
def test_worms_get_json_handles_204(mock_get):
    """
    HTTP 204 is a legitimate WoRMS No Content response and should
    return None rather than raising an exception.
    """

    mock_get.return_value = make_response(
        status_code=204,
    )

    result = _worms_get_json("https://example.org/test")

    assert result is None
    assert mock_get.call_count == 1


@patch("taxonweave.worms.requests.get")
def test_get_synonyms_handles_204_no_content(mock_get):
    """
    Regression test for the Polydora websteri case.

    WoRMS returns HTTP 204 for the synonym endpoint when no synonym
    records are available. TaxonWeave should interpret this as an
    empty synonym list rather than as an API failure.
    """

    mock_get.return_value = make_response(
        status_code=204,
    )

    result = get_synonyms(153847)

    assert result == []
    assert mock_get.call_count == 1


# ---------------------------------------------------------------------
# Retry behavior
# ---------------------------------------------------------------------


@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_retries_after_429(mock_get, mock_sleep):
    """
    HTTP 429 Too Many Requests should trigger a retry.
    """

    rate_limited = make_response(
        status_code=429,
        text="Too Many Requests",
    )

    successful = make_response(
        status_code=200,
        json_data={"success": True},
    )

    mock_get.side_effect = [
        rate_limited,
        successful,
    ]

    result = _worms_get_json("https://example.org/test")

    assert result == {"success": True}
    assert mock_get.call_count == 2
    assert mock_sleep.call_count == 1


@pytest.mark.parametrize(
    "status_code",
    [500, 502, 503, 504],
)
@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_retries_after_server_error(
    mock_get,
    mock_sleep,
    status_code,
):
    """
    Common transient WoRMS/server errors should trigger a retry.
    """

    failed = make_response(
        status_code=status_code,
        text="Server error",
    )

    successful = make_response(
        status_code=200,
        json_data={"success": True},
    )

    mock_get.side_effect = [
        failed,
        successful,
    ]

    result = _worms_get_json("https://example.org/test")

    assert result == {"success": True}
    assert mock_get.call_count == 2
    assert mock_sleep.call_count == 1


# ---------------------------------------------------------------------
# Empty HTTP 200 response
# ---------------------------------------------------------------------


@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_empty_200_response_raises_error(
    mock_get,
    mock_sleep,
):
    """
    HTTP 200 with an empty body is different from HTTP 204.

    A 200 response is expected to contain JSON. If it is empty,
    TaxonWeave should retry and ultimately raise WoRMSError.
    """

    empty_response = make_response(
        status_code=200,
        text="",
    )

    mock_get.return_value = empty_response

    with pytest.raises(WoRMSError, match="empty"):
        _worms_get_json("https://example.org/test")

    assert mock_get.call_count == 3
    assert mock_sleep.call_count == 2


# ---------------------------------------------------------------------
# Invalid JSON
# ---------------------------------------------------------------------


@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_invalid_json_raises_error(
    mock_get,
    mock_sleep,
):
    """
    HTTP 200 containing non-JSON content should be retried and should
    ultimately produce a controlled WoRMSError.
    """

    response = make_response(
        status_code=200,
        text="<html>Service unavailable</html>",
    )

    response.json.side_effect = requests.exceptions.JSONDecodeError(
        "Expecting value",
        "<html>",
        0,
    )

    mock_get.return_value = response

    with pytest.raises(WoRMSError, match="invalid or non-JSON"):
        _worms_get_json("https://example.org/test")

    assert mock_get.call_count == 3
    assert mock_sleep.call_count == 2


# ---------------------------------------------------------------------
# Timeout
# ---------------------------------------------------------------------


@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_timeout_raises_controlled_error(
    mock_get,
    mock_sleep,
):
    """
    Repeated network timeouts should result in a controlled WoRMSError.
    """

    mock_get.side_effect = requests.exceptions.Timeout(
        "Request timed out"
    )

    with pytest.raises(WoRMSError, match="could not connect"):
        _worms_get_json("https://example.org/test")

    assert mock_get.call_count == 3
    assert mock_sleep.call_count == 2


# ---------------------------------------------------------------------
# Connection failure
# ---------------------------------------------------------------------


@patch("taxonweave.worms.time.sleep", return_value=None)
@patch("taxonweave.worms.requests.get")
def test_worms_connection_error_raises_controlled_error(
    mock_get,
    mock_sleep,
):
    """
    Repeated connection failures should result in a controlled
    WoRMSError rather than exposing the raw requests exception.
    """

    mock_get.side_effect = requests.exceptions.ConnectionError(
        "Connection failed"
    )

    with pytest.raises(WoRMSError, match="could not connect"):
        _worms_get_json("https://example.org/test")

    assert mock_get.call_count == 3
    assert mock_sleep.call_count == 2