"""Lifted by build_rggmci_package.py from Sapote-Mamey source. Do not edit here."""
from __future__ import annotations
import json, os, re, zipfile, io, warnings
from pathlib import Path
from .ziputil import regular_file_names


# --- from mamey/parsers.py: FASTA_EXTS, GBK_EXTS, GbkSizeGuardRefusal, QUALIFIER_MAX_CHARS, _GBK_MAX_COMPRESSION_RATIO, _GBK_MAX_UNCOMPRESSED_BYTES, _cap_qual, _contig_length_map_from_zip, _edge_status, _feature_products, _gbk_guard_limits, _gbk_size_guard, _guarded_read_bytes, _record_contig_id, _region_orig_bounds_from_zip, _replicon_intake_key, _require_seqio, _safe_read_text, is_macos_cruft, read_fasta_sequences_from_zip, read_genbank_records
FASTA_EXTS = (".fasta", ".fa", ".fna", ".ffn")

GBK_EXTS = (".gbk", ".gbff", ".gb")

_GBK_MAX_UNCOMPRESSED_BYTES = 100_000_000       # 100 MB: far above any real region/full-assembly GBK

_GBK_MAX_COMPRESSION_RATIO = 200                # uncompressed/compressed ratio above this = bomb-like

def _gbk_guard_limits() -> tuple[int, int]:
    """``(max_bytes, max_ratio)`` for the GBK size guard, read from the environment at call time."""
    try:
        max_bytes = int(os.environ.get("MAMEY_GBK_MAX_BYTES", str(_GBK_MAX_UNCOMPRESSED_BYTES)))
    except (TypeError, ValueError):
        max_bytes = _GBK_MAX_UNCOMPRESSED_BYTES
    try:
        max_ratio = int(os.environ.get("MAMEY_GBK_MAX_RATIO", str(_GBK_MAX_COMPRESSION_RATIO)))
    except (TypeError, ValueError):
        max_ratio = _GBK_MAX_COMPRESSION_RATIO
    return max_bytes, max_ratio

def _gbk_size_guard(info: "zipfile.ZipInfo") -> str | None:
    """Return a refusal reason if this zip member is too large / bomb-like to read, else None."""
    max_bytes, max_ratio = _gbk_guard_limits()
    size = getattr(info, "file_size", 0) or 0
    comp = getattr(info, "compress_size", 0) or 0
    if size > max_bytes:
        return f"uncompressed size {size} exceeds MAMEY_GBK_MAX_BYTES={max_bytes}"
    if comp > 0 and size / comp > max_ratio:
        return f"compression ratio {size // max(comp, 1)}x exceeds MAMEY_GBK_MAX_RATIO={max_ratio} (decompression-bomb guard)"
    return None

class GbkSizeGuardRefusal(ValueError):
    """Typed refusal: a zip member failed the GBK size / compression-ratio guard and was NOT loaded.

    v9.7.410 (CLAUDE_410_gbk_size_guard_all_reads): the .409 guard was wired to only the SeqIO
    read in read_genbank_records(); every other member read (_safe_read_text → FASTA / region
    bounds / version JSON; antismash_input._first_region_header; detect_strictness saccharide
    fallback) still pulled the full uncompressed member into RAM. All of them now go through
    _guarded_read_bytes() and surface this one typed error (message prefix
    ``GBK_SIZE_GUARD_REFUSED:``) instead of a silent full load.
    """

    def __init__(self, member: str, reason: str):
        super().__init__(f"GBK_SIZE_GUARD_REFUSED: {member}: {reason}")
        self.member = member
        self.reason = reason

def _guarded_read_bytes(zf: zipfile.ZipFile, name: str, limit: int | None = None) -> bytes:
    """Read a zip member with the GBK size guard applied BEFORE and DURING the read.

    1. ``getinfo()`` preflight through :func:`_gbk_size_guard` (declared size / ratio).
    2. Bounded streaming read: at most ``limit`` bytes when given (header peeks), otherwise
       ``cap + 1`` bytes so a member whose central-directory size lies about its real size is
       still refused rather than decompressed in full.
    Raises :class:`GbkSizeGuardRefusal`; never returns more than the cap.
    """
    info = zf.getinfo(name)
    reason = _gbk_size_guard(info)
    if reason:
        raise GbkSizeGuardRefusal(name, reason)
    max_bytes, _ = _gbk_guard_limits()
    with zf.open(name) as fh:
        if limit is not None:
            return fh.read(max(0, min(limit, max_bytes)))
        data = fh.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise GbkSizeGuardRefusal(
            name,
            f"streamed size exceeds MAMEY_GBK_MAX_BYTES={max_bytes} "
            f"(central directory declared {info.file_size} bytes)")
    return data

def _safe_read_text(zf: zipfile.ZipFile, name: str) -> str:
    """Decode a zip member as UTF-8 text. Raises GbkSizeGuardRefusal for an over-cap / bomb-like
    member instead of loading it (v9.7.410; previously an unbounded ``zf.read``)."""
    return _guarded_read_bytes(zf, name).decode("utf-8", errors="replace")

def is_macos_cruft(name: str) -> bool:
    """True for macOS archive artifacts that are not real package files.

    Covers ``__MACOSX/`` resource-fork trees, AppleDouble ``._`` sidecars, and
    ``.DS_Store``. These satisfy naive extension/substring filters (e.g. a
    ``._NODE_x.region001.gbk`` shadow ends in ``.gbk`` and contains 'region'),
    so any consumer that counts or parses ZIP entries by a bare
    ``endswith()`` + substring predicate would otherwise double-count them or
    pollute empty/errored-record diagnostics. Strip at the point names first
    enter the system so no downstream counter has to know about AppleDouble.
    (v9.7.152 / engine 1.9.101 — see AUDIT_AS-XXX_macosx_appledouble.)
    """
    base = name.rsplit("/", 1)[-1]
    return (
        name.startswith("__MACOSX/")
        or "/__MACOSX/" in name
        or base.startswith("._")
        or base == ".DS_Store"
    )

def _record_contig_id(rec) -> str:
    """Return the physical contig identifier from a GenBank record.

    Biopython uses the VERSION line as ``rec.id``. For SPAdes-style antiSMASH
    region GBKs, a contig such as ``NODE_1_length_457136_cov_55.054054`` can
    become ``NODE_1_length_457136_cov_55.54054`` because ACCESSION/VERSION
    treats the coverage decimal as a sequence version. The LOCUS token survives
    as ``rec.name``. Prefer it for NODE_*_length_*_cov_* records so FASTA/GBK
    joins and true contig-length lookup remain exact. For normal accession
    records, keep Biopython's ``rec.id`` behavior.
    """
    rec_id = str(getattr(rec, "id", "") or "").strip()
    rec_name = str(getattr(rec, "name", "") or "").strip()
    if rec_name and rec_name not in {".", "<unknown name>", "<unknown>"}:
        if re.match(r"^NODE_\d+_length_\d+_cov_", rec_name, flags=re.I):
            return rec_name
    return rec_id or rec_name

def _replicon_intake_key(item) -> tuple[int, str]:
    """Order chromosome records before plasmids without rereading the archive.

    antiSMASH exports do not provide one universal replicon field, so this uses
    only explicit text already present on the parsed record and source member.
    Unknown or contradictory labels retain deterministic lexical order after
    the recognized chromosome and plasmid groups.
    """
    source_name, rec = item
    annotations = getattr(rec, "annotations", {}) or {}
    text = " ".join(
        str(value or "")
        for value in (
            source_name,
            getattr(rec, "id", ""),
            getattr(rec, "name", ""),
            getattr(rec, "description", ""),
            annotations.get("replicon", ""),
            annotations.get("chromosome", ""),
        )
    ).lower()
    has_chromosome = bool(re.search(r"\bchromosome\b", text))
    has_plasmid = bool(re.search(r"\bplasmid\b", text))
    priority = 0 if has_chromosome and not has_plasmid else 1 if has_plasmid and not has_chromosome else 2
    return priority, str(source_name).casefold()

def read_fasta_sequences_from_zip(zip_path: str | Path) -> dict[str, str]:
    """Return the complete admitted FASTA sequence set from an archive.

    A member-level size/refusal guard cannot safely degrade to a partial return:
    every production caller treats a non-empty mapping as the complete assembly
    and suppresses its GenBank fallback.  Propagate the typed refusal so callers
    cannot publish subset-derived assembly metrics, boundary lengths, or scans.
    """
    seqs = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in regular_file_names(zf):
            if name.lower().endswith(FASTA_EXTS) and not is_macos_cruft(name):
                # v9.7.410 correction: _safe_read_text propagates
                # GbkSizeGuardRefusal. Returning the other members would make a
                # partial assembly indistinguishable from a complete one.
                text = _safe_read_text(zf, name)
                cur, parts = None, []
                for line in text.splitlines():
                    if line.startswith(">"):
                        if cur:
                            seqs[cur] = "".join(parts)
                        # PARSE-02: tolerate a bare ">" header (no id token) rather than
                        # raising IndexError on split()[0]; synthesize a stable placeholder id.
                        _tok = line[1:].split()
                        cur = _tok[0] if _tok else f"unnamed_{len(seqs)+1}"
                        parts = []
                    elif cur:
                        parts.append(line.strip())
                if cur:
                    seqs[cur] = "".join(parts)
    return seqs

def _require_seqio():
    """Import Biopython lazily.

    This keeps `python -m mamey validate ...` and JSON/TXT-only utilities usable
    in restricted sandboxes that do not have Biopython installed. GenBank-backed
    run paths still require Biopython and fail with a clear message.
    """
    try:
        from Bio import SeqIO  # type: ignore
        return SeqIO
    except ImportError:
        return None

def read_genbank_records(zip_path: str | Path, region_only: bool = False,
                         exclude_regions: bool = False):
    SeqIO = _require_seqio()
    records = []
    errored: list[str] = []        # GBKs that raised during parse
    empty: list[str] = []          # GBKs that parsed without error but yielded 0 records
    # v9.7.173: fail clearly on a directory input. A sealed Mamey package dir has no GBKs and
    # stores parsed JSON (not translations), so `mamey compare <package_dir>` used to crash with
    # a raw IsADirectoryError from zipfile.ZipFile below. Give an actionable message instead.
    if Path(zip_path).is_dir():
        raise ValueError(
            f"read_genbank_records expects an antiSMASH output ZIP, got a directory: {zip_path}. "
            f"A sealed Mamey package stores parsed JSON (no GenBank records / protein "
            f"translations); pass the strain's antiSMASH output ZIP instead.")
    # v9.7.184 P1: accept a single .gbk/.gb/.gbff file directly. The blastp-online --package help
    # advertises "antiSMASH region GBK/ZIP", but a bare GBK fell straight into ZipFile below and
    # raised BadZipFile. Parse the file in place instead, mirroring the zip path's record shape.
    _p = Path(zip_path)
    if _p.is_file() and _p.suffix.lower() in GBK_EXTS:
        name = _p.name
        if region_only and "region" not in name.lower():
            return []
        try:
            if SeqIO is not None:
                with open(_p, encoding="utf-8", errors="replace") as handle:
                    for rec in SeqIO.parse(handle, "genbank"):
                        records.append((name, rec))
            else:
                from ._gbk_shim import parse_genbank_text
                for rec in parse_genbank_text(_p.read_text(encoding="utf-8", errors="replace")):
                    records.append((name, rec))
        except Exception as exc:
            raise ValueError(f"could not parse GenBank file {name}: {exc}") from exc
        return records
    with zipfile.ZipFile(zip_path) as zf:
        gbks = [n for n in regular_file_names(zf) if n.lower().endswith(GBK_EXTS) and not is_macos_cruft(n)]
        if region_only:
            gbks = [n for n in gbks if "region" in Path(n).name.lower()]
        elif exclude_regions:
            # v9.7.105 fix: the genome-wide CDS inventory must come from the full-assembly
            # GBK(s) only. Per-region GBKs carry REGION-LOCAL coordinates; mixing them in
            # duplicates every CDS at a shifted frame, which piles low-coordinate copies into
            # the first region (closed-genome scoping bug). Every region CDS is already present
            # in the full assembly, so excluding region files loses nothing.
            _full = [n for n in gbks if "region" not in Path(n).name.lower()]
            if _full:                      # only exclude when a full-assembly GBK exists
                gbks = _full
        total_gbks = len(gbks)
        for name in sorted(gbks):
            n_before = len(records)
            # v9.7.409 (DEEP_AUDIT2_resource_dos #3): refuse an over-large / bomb-like GBK BEFORE
            # loading it into RAM. Treated as a skipped (errored) member so the existing skip
            # diagnostic surfaces it — one refused file must not abort the rest of the package.
            try:
                _too_big = _gbk_size_guard(zf.getinfo(name))
            except Exception:
                _too_big = None
            if _too_big:
                errored.append(name)
                continue
            try:
                if SeqIO is not None:
                    with zf.open(name) as raw:
                        handle = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                        for rec in SeqIO.parse(handle, "genbank"):
                            records.append((name, rec))
                else:
                    from ._gbk_shim import parse_genbank_text
                    text = zf.read(name).decode("utf-8", errors="replace")
                    for rec in parse_genbank_text(text):
                        records.append((name, rec))
            except Exception:
                # Resilient by design: one malformed GBK must not abort the others.
                errored.append(name)
                continue
            # A GBK that parses cleanly but produces zero records is *also* a silent
            # drop — and the more insidious one, because no exception is raised. This
            # is exactly what produces downstream BGC-count mismatches the receipt
            # audit then has to diagnose from a distance, so surface it too.
            if len(records) == n_before:
                empty.append(name)
    skipped = errored + empty
    if skipped:
        def _names(lst):
            return ", ".join(Path(f).name for f in lst[:8]) + (" …" if len(lst) > 8 else "")
        detail = []
        if errored:
            detail.append(f"{len(errored)} errored [{_names(errored)}]")
        if empty:
            detail.append(f"{len(empty)} parsed-empty [{_names(empty)}]")
        warnings.warn(
            f"read_genbank_records: {len(skipped)}/{total_gbks} GBK file(s) yielded no "
            f"records in {zip_path}: " + "; ".join(detail),
            stacklevel=2,
        )
    return records

def _feature_products(feature) -> list[str]:
    vals = []
    # Region/protocluster/candidate-cluster category is an antiSMASH umbrella,
    # not another product class. Preserve the legacy CDS-only fallback below.
    keys = (("product", "products", "category", "aSDomain", "domain")
            if getattr(feature, "type", None) == "CDS" else ("product", "products"))
    for key in keys:
        if key in feature.qualifiers:
            q = feature.qualifiers[key]
            vals.extend(q if isinstance(q, list) else [str(q)])
    out = []
    for v in vals:
        for part in re.split(r"[,;/]+", _cap_qual(str(v))):
            part = _cap_qual(part.strip())   # H12: a product token is a class label; never a 300 KB figure label
            if part and part not in out:
                out.append(part)
    return out

def _edge_status(start: int, end: int, contig_len: int, flank_bp: int = 5000,
                 is_circular: bool = False) -> str:
    if contig_len <= 0:
        return "Unknown"
    length = end - start + 1
    if length >= 0.95 * contig_len:
        return "Full-contig"
    # v9.7.87 P0-a: on a closed/circular replicon the origin is not a truncation point — a BGC
    # adjacent to the origin wraps, it is not Edge. Only call Edge near a contig boundary when the
    # topology is linear (or unknown). This fixes closed genomes getting a false Edge at the origin
    # (parsing-control FAIL, corrected < raw); e.g. N. nova NZ_CP006850, Solwaraspora_WMMA2065.
    if is_circular:
        return "Interior"
    if start <= flank_bp or (contig_len - end) <= flank_bp:
        return "Edge"
    return "Interior"

def _region_orig_bounds_from_zip(zip_path: str | Path) -> dict[str, tuple[int, int]]:
    """Return absolute antiSMASH region bounds keyed by GBK filename.

    antiSMASH region GBKs are clipped records whose feature coordinates restart
    at 1.  The absolute source-record coordinates are stored in the comment as
    ``Orig. start`` / ``Orig. end``.  Boundary classification must use these
    absolute coordinates plus the true contig length; otherwise a closed
    chromosome with many clipped region GBKs is falsely labeled Edge/Full-contig.
    """
    out: dict[str, tuple[int, int]] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in regular_file_names(zf):
            if is_macos_cruft(name):
                continue
            if not name.lower().endswith(GBK_EXTS) or "region" not in Path(name).name.lower():
                continue
            try:
                text = _safe_read_text(zf, name)
            except GbkSizeGuardRefusal as exc:
                # v9.7.410: refused before any load (the same member is refused again, and
                # counted as errored, by read_genbank_records). Warn so the skip is not silent.
                warnings.warn(f"_region_orig_bounds_from_zip: {exc}; member skipped", stacklevel=2)
                continue
            except Exception:
                continue
            # v9.7.438: GenBank marks a partial boundary with `<` or `>` and antiSMASH passes it
            # through, e.g. `Orig. end :: >247156`. `(\d+)` will not match across the marker, so the
            # pair was dropped, the caller fell back to the clipped record's LOCAL coordinates
            # (which restart at 1), and `_edge_status` then saw `start <= flank_bp` and returned
            # Edge. Observed on a region sitting 180 kb clear of both ends of a 428 kb contig
            # (N. macrotermitis NZ_WEGK01000006.1 region003). The number after the marker is the
            # real coordinate -- `>247156` means truncated AT 247156, not unknown -- so accept and
            # discard the marker. Circular records were unaffected: _edge_status returns Interior
            # for them before it reads the start coordinate.
            sm = re.search(r"Orig\.\s*start\s*::\s*[<>]?(\d+)", text, flags=re.I)
            em = re.search(r"Orig\.\s*end\s*::\s*[<>]?(\d+)", text, flags=re.I)
            if sm and em:
                out[name] = (int(sm.group(1)), int(em.group(1)))
    return out

def _contig_length_map_from_zip(zip_path: str | Path) -> dict[str, int]:
    """True source-record contig lengths from FASTA/full GenBank records."""
    seqs = read_fasta_sequences_from_zip(zip_path)
    if seqs:
        return {k: len(v) for k, v in seqs.items()}
    lengths: dict[str, int] = {}
    for name, rec in read_genbank_records(zip_path, region_only=False):
        # Skip clipped antiSMASH region GBKs; they are not contig-length records.
        if "region" in Path(name).name.lower():
            continue
        _cid = _record_contig_id(rec)
        lengths[_cid] = max(lengths.get(_cid, 0), len(rec.seq))
    return lengths

QUALIFIER_MAX_CHARS = 4000

def _cap_qual(value):
    """Truncate an over-long qualifier string; leaves non-strings and normal strings untouched."""
    if isinstance(value, str) and len(value) > QUALIFIER_MAX_CHARS:
        return value[:QUALIFIER_MAX_CHARS] + f" …[truncated {len(value) - QUALIFIER_MAX_CHARS} chars]"
    return value
