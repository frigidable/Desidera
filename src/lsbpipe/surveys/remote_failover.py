"""Strict provider endpoint routing with reproducible failover provenance.

Only explicitly reviewed endpoints may be attempted. Returning a successful
empty query is NOT grounds for trying another server. This module does no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, TypeVar
from urllib.error import HTTPError, URLError

T = TypeVar("T")

SIMBAD_ENDPOINTS = (
    "simbad.cds.unistra.fr", "simbad.harvard.edu",
)
SIMBAD_POLICIES = (
    "auto", "prefer_harvard", "strasbourg_only", "harvard_only",
)


def simbad_endpoint_order(policy: str) -> tuple[str, ...]:
    if policy == "auto":
        return SIMBAD_ENDPOINTS
    if policy == "prefer_harvard":
        return SIMBAD_ENDPOINTS[::-1]
    if policy == "strasbourg_only":
        return SIMBAD_ENDPOINTS[:1]
    if policy == "harvard_only":
        return SIMBAD_ENDPOINTS[1:]
    raise ValueError("Unknown SIMBAD endpoint policy; no unapproved mirror")


def retryable_remote_failure(exc: Exception) -> bool:
    """Switch on connection/service failure, not bad query or malformed data."""
    if isinstance(exc, HTTPError):
        return exc.code in (403, 404, 408, 425, 429) or exc.code >= 500
    if isinstance(exc, (URLError, TimeoutError, ConnectionError)):
        return True
    if isinstance(exc, OSError):
        return True
    # Astroquery uses requests and pyvo internally. Keep dependencies
    # optional for offline CI; avoid catching arbitrary programming errors.
    typ = type(exc)
    module, name = typ.__module__, typ.__name__
    if module.startswith(("requests.", "urllib3.")):
        status = getattr(getattr(exc, "response", None), "status_code", None)
        if status is not None:
            return status in (403, 404, 408, 425, 429) or status >= 500
        return name not in ("InvalidURL", "MissingSchema")
    if module.startswith("pyvo.") and name == "DALServiceError":
        return True
    if module.startswith("astroquery.") and name in (
            "RemoteServiceError", "TimeoutError"):
        return True
    return False


@dataclass(frozen=True)
class EndpointAttemptFailure:
    host: str
    error_type: str

    def as_dict(self) -> dict:
        return {"host": self.host, "error_type": self.error_type}


class AllRemoteEndpointsFailed(RuntimeError):
    def __init__(self, service: str, failures: list[EndpointAttemptFailure]):
        self.service = service
        self.failures = tuple(failures)
        super().__init__(
            f"{service} unreachable on all approved endpoints: " +
            ", ".join(f"{f.host} ({f.error_type})" for f in failures))


def with_simbad_failover(
        action: Callable[[str], T], *, policy: str = "auto"
) -> tuple[T, dict]:
    hosts = simbad_endpoint_order(policy)
    failed: list[EndpointAttemptFailure] = []
    for host in hosts:
        try:
            result = action(host)
        except Exception as exc:
            if not retryable_remote_failure(exc):
                raise
            failed.append(EndpointAttemptFailure(
                host=host, error_type=type(exc).__name__))
            continue
        return result, {
            "provider": "SIMBAD", "endpoint_host": host,
            "endpoint_policy": policy,
            "attempted_hosts": list(hosts[:len(failed) + 1]),
            "failed_attempts": [f.as_dict() for f in failed],
            "fallback_used": bool(failed),
        }
    raise AllRemoteEndpointsFailed("SIMBAD", failed)
