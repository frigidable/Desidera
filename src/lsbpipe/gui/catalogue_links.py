"""Deterministic human-viewer links for verified SIMBAD/NED identities.

URLs use only official HTTPS domains and the exact *confirmed* catalogue
identifier. These are navigation helpers, not catalogue query adapters.
Unverified cone-search hits must never be promoted into hyperlinks.
"""
from __future__ import annotations

from html import escape
from urllib.parse import urlencode

SIMBAD_STRASBOURG = "simbad.cds.unistra.fr"
SIMBAD_HARVARD = "simbad.harvard.edu"
NED_HOST = "ned.ipac.caltech.edu"
_MAX_ID_LENGTH = 256


def confirmed_object_url(
    provider: str, main_id: str, *, endpoint_host: str | None = None
) -> str:
    """Construct official object-viewer URL for an already-confirmed identity.

    This function does not verify an identity itself. Only call it after
    immutable confirmation evidence has been validated by the read model.
    """
    if not isinstance(main_id, str):
        raise ValueError("Confirmed catalogue identifier is required")
    identifier = main_id.strip()
    if (not identifier or len(identifier) > _MAX_ID_LENGTH or
            any(ord(char) < 32 or ord(char) == 127 for char in identifier)):
        raise ValueError("Invalid confirmed catalogue identifier")
    if provider == "SIMBAD":
        # Use Harvard when that was the recorded *official* SIMBAD endpoint,
        # otherwise default to Strasbourg. Never trust arbitrary hosts from
        # external catalogue or provenance fields.
        host = (SIMBAD_HARVARD if endpoint_host == SIMBAD_HARVARD
                else SIMBAD_STRASBOURG)
        return f"https://{host}/simbad/sim-id?" + urlencode(
            {"Ident": identifier})
    if provider == "NED":
        return f"https://{NED_HOST}/byname?" + urlencode(
            {"objname": identifier})
    raise ValueError("Unsupported catalogue provider")


def confirmed_name_link(
    provider: str, main_id: str, *, endpoint_host: str | None = None
) -> str:
    """HTML-escaped anchor for Qt rich-text QLabel, with fixed trusted origin."""
    url = confirmed_object_url(
        provider, main_id, endpoint_host=endpoint_host)
    return (f'<a href="{escape(url, quote=True)}">'
            f'{escape(main_id)} ↗</a>')
