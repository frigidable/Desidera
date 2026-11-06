"""DESIDERA public facade over the preserved lsbpipe implementation.

Legacy imports and scientific provenance remain unchanged.
"""
from lsbpipe import *  # noqa: F401,F403
from lsbpipe import __all__ as _legacy_public_api
from lsbpipe import __version__

__all__ = [*_legacy_public_api, "__version__"]
