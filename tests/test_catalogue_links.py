"""P10f navigation links: no web calls or unverified scientific adoption."""
from html import unescape
from urllib.parse import parse_qs, urlsplit

import pytest

from lsbpipe.gui.catalogue_links import (
    confirmed_object_url, confirmed_name_link,
)


@pytest.mark.parametrize("provider,name,host,path,key", [
    ("SIMBAD", "NAME Leo P", "simbad.cds.unistra.fr",
     "/simbad/sim-id", "Ident"),
    ("SIMBAD", "[ABC2019] LP 26", "simbad.cds.unistra.fr",
     "/simbad/sim-id", "Ident"),
    ("NED", "PGC1 5065058 NED001", "ned.ipac.caltech.edu",
     "/byname", "objname"),
    ("NED", "Leo P", "ned.ipac.caltech.edu",
     "/byname", "objname"),
])
def test_official_https_name_links_roundtrip(provider, name, host, path, key):
    uri = confirmed_object_url(provider, name)
    parts = urlsplit(uri)
    assert parts.scheme == "https"
    assert parts.hostname == host
    assert parts.path == path
    assert parse_qs(parts.query) == {key: [name]}
    assert " " not in uri
    assert confirmed_name_link(provider, name).count("<a href=") == 1


def test_simbad_harvard_is_allowed_but_arbitrary_snapshot_host_is_not():
    harvard = confirmed_object_url(
        "SIMBAD", "NAME Leo P", endpoint_host="simbad.harvard.edu")
    assert urlsplit(harvard).hostname == "simbad.harvard.edu"
    for forged in ["evil.example", "simbad.harvard.edu.evil.example",
                   "http://simbad.harvard.edu", "simbad.cds.unistra.fr@evil"]:
        url = confirmed_object_url(
            "SIMBAD", "NAME Leo P", endpoint_host=forged)
        assert urlsplit(url).hostname == "simbad.cds.unistra.fr"


def test_identifier_metacharacters_are_escaped_as_text_and_url_query():
    malicious = '[TEST] 1"><img src=x onerror="alert(1)"'
    rich = confirmed_name_link("SIMBAD", malicious)
    assert "<img" not in rich
    assert "&lt;img" in rich
    assert "&quot;" in rich
    assert "%3C" in rich
    assert "onclick=" not in rich
    assert parse_qs(urlsplit(confirmed_object_url(
        "SIMBAD", malicious)).query)["Ident"] == [malicious]
    assert rich.startswith('<a href="https://simbad.cds.unistra.fr/')
    assert rich.endswith(' ↗</a>')


@pytest.mark.parametrize("provider,identifier", [
    ("SIMBAD", ""),
    ("NED", " "),
    ("NED", "NGC\n0001"),
    ("SIMBAD", "A" * 257),
    ("unknown", "Leo P"),
])
def test_invalid_identity_link_inputs_fail_closed(provider, identifier):
    with pytest.raises(ValueError):
        confirmed_object_url(provider, identifier)


def test_link_render_has_no_io(monkeypatch):
    import socket
    monkeypatch.setattr(socket, "create_connection", lambda *a, **k:
                        pytest.fail("Rendering a catalogue link opened network"))
    assert "Leo P" in confirmed_name_link("NED", "Leo P")
