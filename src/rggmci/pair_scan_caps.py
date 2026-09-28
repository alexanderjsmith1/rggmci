"""pair_scan_caps.py — bound all-vs-all O(n^2) BGC-pair (and O(n^3) ANI) scans.

v9.7.409 (DEEP_AUDIT2_resource_dos #3). Several comparative scans pair every BGC with every other
BGC (``scan_efls``, ``comparative_pairs``, ``efls``, ``rggmci``, cross-strain threads) and the ANI
heatmap does a triple-nested intransitivity check over strains. Real actinomycete genomes carry
~10-50 BGCs, so these are cheap in practice — but the BGC/strain count is NOT capped at parse, so a
degenerate or crafted input declaring thousands of regions turns an O(n^2) (or O(n^3)) loop and its
in-memory candidate list into a resource-exhaustion vector. Truncating the *output* (several call
sites already slice ``[:500]``) does not bound the *work*.

This module bounds the WORK by capping the number of items fed into the pair loop, emitting a loud,
typed ``PAIR_SCAN_CAPPED`` warning to stderr when it trims. It is a robustness bound only: it changes
nothing on any real-sized input (well under the cap), and when it does trim it says so rather than
silently dropping pairs. Caps are env-overridable so a genuinely large cohort can raise them.

Claim-safety: this module carries no biological meaning. It bounds compute; it makes no claim.
"""
from __future__ import annotations

import os
import sys
from typing import Sequence, TypeVar

T = TypeVar("T")

#: Max items in an all-vs-all pair scan. 800 -> ~320k pairs, tractable; far above any real genome.
_DEFAULT_MAX_PAIR_SCAN_ITEMS = 800
#: Max strains in the O(n^3) ANI intransitivity check. 200 -> 8M triples worst case; still bounded.
_DEFAULT_MAX_ANI_STRAINS = 200


def _read_int(name: str, default: int) -> int:
    try:
        val = int(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default
    return val if val > 0 else default


def max_pair_scan_items() -> int:
    """Cap on items entering an all-vs-all pair scan (env: MAMEY_MAX_PAIR_SCAN_ITEMS)."""
    return _read_int("MAMEY_MAX_PAIR_SCAN_ITEMS", _DEFAULT_MAX_PAIR_SCAN_ITEMS)


def max_ani_strains() -> int:
    """Cap on strains entering the O(n^3) ANI intransitivity scan (env: MAMEY_MAX_ANI_STRAINS)."""
    return _read_int("MAMEY_MAX_ANI_STRAINS", _DEFAULT_MAX_ANI_STRAINS)


def cap_pair_scan_items(items: Sequence[T], *, label: str, cap: int | None = None) -> list[T]:
    """Return ``items`` bounded to ``cap`` (default :func:`max_pair_scan_items`).

    When ``items`` is longer than the cap, emit a typed ``PAIR_SCAN_CAPPED`` warning to stderr and
    return only the first ``cap`` entries, so the caller's O(n^2) loop is bounded. Deterministic:
    it trims a prefix of the caller's already-ordered sequence, never reorders.
    """
    limit = cap if cap is not None else max_pair_scan_items()
    seq = list(items)
    if len(seq) <= limit:
        return seq
    print(
        f"PAIR_SCAN_CAPPED: {label} received {len(seq)} items; the all-vs-all pair scan is capped "
        f"to the first {limit} to bound O(n^2) work (raise MAMEY_MAX_PAIR_SCAN_ITEMS to include "
        f"more). Pairs beyond the cap are not evaluated.",
        file=sys.stderr,
    )
    return seq[:limit]
