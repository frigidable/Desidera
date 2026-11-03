"""P0.9a — deterministic network failover, with NO real network in CI."""
from __future__ import annotations

import sys
from types import ModuleType

import pytest

from lsbpipe.surveys.remote_failover import (
    SIMBAD_ENDPOINTS, AllRemoteEndpointsFailed,
    retryable_remote_failure, simbad_endpoint_order, with_simbad_failover,
)
from lsbpipe.surveys.galaxy_knowledge import query_simbad_candidates
from astropy.table import Table


def test_only_explicit_supported_simbad_mirrors():
    assert simbad_endpoint_order("auto") == SIMBAD_ENDPOINTS
    assert simbad_endpoint_order("prefer_harvard") == SIMBAD_ENDPOINTS[::-1]
    assert simbad_endpoint_order("harvard_only") == SIMBAD_ENDPOINTS[1:]
    assert simbad_endpoint_order("strasbourg_only") == SIMBAD_ENDPOINTS[:1]
    with pytest.raises(ValueError, match="unapproved mirror"):
        simbad_endpoint_order("https://attacker.example/simbad")


def test_timeout_fails_over_and_records_exact_origin():
    seen = []

    def action(host):
        seen.append(host)
        if host == SIMBAD_ENDPOINTS[0]:
            raise TimeoutError("offline")
        return {"results": [42]}

    output, receipt = with_simbad_failover(action)
    assert output["results"] == [42]
    assert seen == list(SIMBAD_ENDPOINTS)
    assert receipt["endpoint_host"] == SIMBAD_ENDPOINTS[1]
    assert receipt["attempted_hosts"] == list(SIMBAD_ENDPOINTS)
    assert receipt["failed_attempts"] == [{
        "host": SIMBAD_ENDPOINTS[0], "error_type": "TimeoutError"}]
    assert receipt["fallback_used"] is True


def test_prefer_harvard_success_skips_strasbourg():
    seen = []
    output, info = with_simbad_failover(
        lambda host: seen.append(host) or None, policy="prefer_harvard")
    assert output is None  # legitimate no counterpart, not a failure
    assert seen == [SIMBAD_ENDPOINTS[1]]
    assert info["fallback_used"] is False


def test_invalid_response_is_not_silently_hidden_by_mirror():
    seen = []

    def broken(host):
        seen.append(host)
        raise ValueError("inconsistent coordinates")

    with pytest.raises(ValueError, match="coordinates"):
        with_simbad_failover(broken)
    assert seen == [SIMBAD_ENDPOINTS[0]]


def test_all_offline_raises_without_publishing_science_data():
    with pytest.raises(AllRemoteEndpointsFailed) as err:
        with_simbad_failover(
            lambda host: (_ for _ in ()).throw(ConnectionError(host)))
    assert len(err.value.failures) == 2
    assert all(f.error_type == "ConnectionError" for f in err.value.failures)


def test_only_network_status_codes_retry():
    from urllib.error import HTTPError
    assert retryable_remote_failure(
        HTTPError("https://simbad.cds.unistra.fr", 503, "unavailable", {}, None))
    assert retryable_remote_failure(
        HTTPError("https://simbad.cds.unistra.fr", 429, "rate", {}, None))
    assert not retryable_remote_failure(
        HTTPError("https://simbad.cds.unistra.fr", 400, "bad query", {}, None))


def test_real_simbad_adapter_uses_second_server_after_primary_failure(
        monkeypatch):
    """Fake installed astroquery module to test the actual dispatch boundary."""
    seen = []

    class FakeSimbad:
        def __init__(self):
            self.server = None
            self.timeout = None
            self.ROW_LIMIT = None

        def query_object(self, name, wildcard=False):
            seen.append((self.server, "name", self.timeout))
            if self.server == SIMBAD_ENDPOINTS[0]:
                raise TimeoutError("host unreachable")
            return Table(
                rows=[("UGC 2162", 40.096208, 1.22925)],
                names=("main_id", "ra", "dec"))

        def query_region(self, centre, radius):
            seen.append((self.server, "cone", self.timeout))
            return Table(
                rows=[("UGC 2162", 40.096208, 1.22925)],
                names=("main_id", "ra", "dec"))

    fake_module = ModuleType("astroquery.simbad")
    fake_module.Simbad = FakeSimbad
    fake_root = ModuleType("astroquery")
    fake_root.__path__ = []
    monkeypatch.setitem(sys.modules, "astroquery", fake_root)
    monkeypatch.setitem(sys.modules, "astroquery.simbad", fake_module)
    result = query_simbad_candidates(
        name="UGC 2162", ra_deg=40.096208, dec_deg=1.22925)
    assert [(h, k) for h, k, timeout in seen] == [
        (SIMBAD_ENDPOINTS[0], "name"),
        (SIMBAD_ENDPOINTS[1], "name"),
        (SIMBAD_ENDPOINTS[1], "cone")]
    assert all(timeout == 15 for _, _, timeout in seen)
    assert result["endpoint_evidence"]["endpoint_host"] == SIMBAD_ENDPOINTS[1]
    assert result["endpoint_evidence"]["fallback_used"]
    assert len(result["candidates"]) == 1
    assert set(result["candidates"][0]["discovered_via"]) == {"cone", "name"}


def test_injected_client_invalid_policy_refuses_work():
    class Client:
        def query_object(self, *a, **k):
            pytest.fail("Invalid policy should prevent all queries")
    with pytest.raises(ValueError, match="unapproved mirror"):
        query_simbad_candidates(
            name="UGC 2162", ra_deg=None, dec_deg=None,
            client=Client(), mirror_policy="unexpected-domain")
