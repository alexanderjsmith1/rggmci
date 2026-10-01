"""RG-GMCI: reference-guided linkage candidates for fragmented antiSMASH assemblies.

A homology-guided candidate-linkage method, not a physical contig joiner. A linked pair is a
hypothesis to check at gene level, not evidence that two fragments are one pathway.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .core import compute_rggmci, parse_clusterblast_reference_map, run_rggmci  # noqa: F401


def _provenance() -> dict[str, Any]:
    import json
    from importlib import resources
    try:
        return json.loads(resources.files(__package__).joinpath("PROVENANCE.json").read_text())
    except (FileNotFoundError, ValueError, OSError):
        return {}


__version__ = _provenance().get("package_version", "unknown")
from .groups import both_at_contig_ends, candidate_groups
from .reader import read_regions


def _input_record(zip_path: str | Path) -> dict[str, Any]:
    """Which ZIP and which package build produced a result, so a result file can be traced on its own."""
    import hashlib
    prov = _provenance()
    h = hashlib.sha256()
    with open(zip_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return {"antismash_zip": Path(zip_path).name, "antismash_zip_sha256": h.hexdigest(), "rggmci_version": __version__,
            "rggmci_built_from_bundle": prov.get("built_from_bundle"), "rggmci_bundle_build": prov.get("build"),
            "rggmci_engine": prov.get("engine")}


def run(zip_path: str | Path, **options: Any) -> dict[str, Any]:
    """Read an antiSMASH ZIP and score every region pair.

    Returns the scorer's own result, plus `regions`, a `both_at_contig_ends` description on each pair and
    `candidate_groups`. The added keys never change a pair's score or confidence. `options` pass to the scorer: the
    reference-guided completion settings (`reference_completion`, `mibig_db`, `sensitivity`, `threads`, `label`,
    `pfam_hmm`) and
    the residue-tiling ones (`diamond_db`, `diamond`, `residue_scope`, `residue_work_dir`).
    """
    bgcs = read_regions(zip_path)
    options.setdefault("label", Path(zip_path).stem)
    result = run_rggmci(zip_path, bgcs, **options)
    result["input"] = _input_record(zip_path)
    result["regions"] = [{"bgc_id": b.bgc_id, "contig": b.contig, "region": f"region{b.region_number:03d}",
                          "products": b.products, "edge_status": b.edge_status} for b in bgcs]
    pairs = result.get("ranked_pairs") or []
    for p in pairs:
        p["both_at_contig_ends"] = both_at_contig_ends(p)
    result["candidate_groups"] = candidate_groups(pairs)
    return result
