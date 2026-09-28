"""Lifted by build_rggmci_package.py from Sapote-Mamey source. Do not edit here."""
from __future__ import annotations
import re
from collections import defaultdict
from typing import Iterable, Any
import csv
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


# --- from mamey/bgc_blastp_panel.py: _residues, _safe_token, assign_rounds, fasta_header, wrap_fasta
def _safe_token(s: Any, fallback: str = "NA") -> str:
    t = re.sub(r"[^A-Za-z0-9_.:+-]+", "_", str(s or "")).strip("_")
    return t or fallback

def fasta_header(strain: str, row: dict) -> str:
    fields = [
        _safe_token(strain),
        _safe_token(row.get("bgc_id")),
        f"slot={row.get('slot')}",
        f"role={_safe_token(row.get('selection_role'))}",
        f"gene={_safe_token(row.get('locus_tag') or row.get('protein_id'))}",
        f"node={_safe_token(row.get('node_id') or row.get('contig'))}",
        f"region={_safe_token(row.get('antismash_region'))}",
        f"coords={row.get('start')}-{row.get('end')}",
        f"aa={row.get('aa_len')}",
        f"reason={_safe_token(row.get('selection_reason'))}",
    ]
    return ">" + "|".join(fields)

def wrap_fasta(header: str, seq: str, width: int = 60) -> str:
    lines = [header]
    lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + "\n"

def _residues(row: dict) -> int:
    return len(row.get("sequence", ""))

def assign_rounds(rows: list[dict], strain: str, proteins_per_file: int = 20,
                  max_residues: int = 85000, giant_aa_threshold: int = 2500,
                  preserve_bgc_boundaries: bool = True,
                  isolate_giants: bool = False) -> list[list[dict]]:
    """Pack rows into NCBI-safe BLASTP rounds.

    NCBI web BLASTP rejects queries by total amino-acid residues, so the load-
    bearing cap is ``max_residues``. The record-count cap keeps ChatGPT/manual
    review manageable. Giant proteins are flagged in the manifest by default,
    but not isolated unless ``isolate_giants`` is requested. This keeps the
    default output compact while still supporting reduced troubleshooting batches
    if NCBI reports CPU-limit/timeouts.
    """
    rows = list(rows)
    if not rows:
        return []

    groups: list[list[dict]]
    if preserve_bgc_boundaries:
        grouped: dict[str, list[dict]] = defaultdict(list)
        order: list[str] = []
        for row in rows:
            bid = str(row.get("bgc_id", ""))
            if bid not in grouped:
                order.append(bid)
            grouped[bid].append(row)
        groups = [grouped[bid] for bid in order]
    else:
        groups = [[r] for r in rows]

    rounds: list[list[dict]] = []
    cur: list[dict] = []
    cur_residues = 0

    def flush():
        nonlocal cur, cur_residues
        if cur:
            rounds.append(cur)
            cur = []
            cur_residues = 0

    for group in groups:
        group_residues = sum(_residues(r) for r in group)
        # Optional troubleshooting mode: isolate giant multidomain proteins when
        # NCBI reports CPU-limit or timeout problems. The default is to flag,
        # not isolate, because the residue cap is the load-bearing web limit.
        #
        # v9.7.240 (P4a): this previously emitted EVERY protein of a giant-bearing BGC
        # as its own one-protein round, not just the giant. On AS-XXX with
        # --genes-per-bgc 8 that turned 13 rounds into 48, most of them singletons of a
        # ~200 aa regulator -- unusable as a manual BLASTp handoff. Isolate the giants;
        # let the non-giants pack normally.
        if isolate_giants and any(_residues(r) >= giant_aa_threshold for r in group):
            giants = [r for r in group if _residues(r) >= giant_aa_threshold]
            rest = [r for r in group if _residues(r) < giant_aa_threshold]
            flush()
            for row in giants:
                rounds.append([row])
            if not rest:
                continue
            group = rest
            group_residues = sum(_residues(r) for r in group)
        if cur and (len(cur) + len(group) > proteins_per_file or cur_residues + group_residues > max_residues):
            flush()
        if len(group) > proteins_per_file or group_residues > max_residues:
            # A BGC group itself is too large; split within group and make the
            # break visible via manifest warnings already carried by rows.
            for row in group:
                if cur and (len(cur) >= proteins_per_file or cur_residues + _residues(row) > max_residues):
                    flush()
                cur.append(row)
                cur_residues += _residues(row)
            continue
        cur.extend(group)
        cur_residues += group_residues
    flush()
    return rounds


# --- from mamey/blastp_followup.py: CLAIM_SAFETY, HitRecord, _accession_key, _cov_str, _find_int, _find_text, _looks_like_hit_table_header, _safe_float, _safe_float_or_none, _safe_int, _text, best_hits_by_query, merge_hit_xml, parse_hit_table_csv, parse_query_id, parse_xml2, query_len_from_id
CLAIM_SAFETY = (
    "BLASTP evidence is sequence similarity only; it does not prove compound "
    "identity, pathway completeness, expression, or bioactivity."
)

@dataclass(frozen=True)
class HitRecord:
    query_id: str
    subject_id: str
    pct_identity: float
    align_len: int
    mismatches: int
    gap_opens: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float
    pct_positive: float | None = None
    query_len: int | None = None
    subject_title: str = ""
    subject_accession: str = ""
    subject_sciname: str = ""
    subject_taxid: str = ""
    hit_len: int | None = None

    @property
    def query_coverage(self) -> float | None:
        # BLP-04 (v9.7.338): when the query length is genuinely unknown, DO NOT fall back to the
        # alignment extent (max(qstart, qend, align_len)) as the denominator — align_len / align_len
        # forces coverage to ~1.0 and fabricates near-full coverage on exactly the hits where we
        # know the least. Report unknown coverage as None (emitted as "NA"), never a false 1.0.
        qlen = self.query_len or query_len_from_id(self.query_id)
        if not qlen:
            return None
        return min(1.0, self.align_len / qlen)

def _cov_str(cov: float | None, places: int = 3) -> str:
    """Format a query-coverage value for emission; unknown coverage (None) -> "NA" (BLP-04)."""
    return "NA" if cov is None else f"{cov:.{places}f}"

def _safe_float(x: Any, default: float = 0.0) -> float:
    try:
        return float(str(x).strip())
    except Exception:
        return default

def _safe_float_or_none(x: Any) -> float | None:
    """Like _safe_float but returns None (not a poison sentinel) on parse failure, so a
    malformed column never becomes a fabricated numeric value downstream. v9.7.371 fix: the
    pct_positive call site used _safe_float(x, default=-1.0) -- a *value*, not a None-checkable
    sentinel -- so a column misalignment (this file's own documented NCBI-CSV quoting edge case)
    silently became a fabricated -1.0 percent, rendered verbatim as '-1.00%' in reader-facing
    output instead of the 'NA' that _cov_str() already gives the equivalent coverage field."""
    try:
        return float(str(x).strip())
    except Exception:
        return None

def _safe_int(x: Any, default: int = 0) -> int:
    try:
        return int(float(str(x).strip()))
    except Exception:
        return default

def _text(x: Any) -> str:
    return re.sub(r"\s+", " ", str(x or "")).strip()

def parse_query_id(query_id: str) -> dict[str, str]:
    """Parse Sapote/Mamey FASTA header metadata from a BLASTP query id/title.

    Two header conventions are supported:

    1. Pipe-delimited (this bundle's own ``bgc_blastp_panel`` emission), e.g.
       ``strain|BGC002|slot=1|gene=ctg11_9|node=NODE_11|aa=317|...``.
    2. Space-delimited ``key=value`` (the convention used by manually-built
       NCBI BLASTP rounds, e.g. ``BGC008_ctg162_3 contig=NODE_162 node=NODE_162
       start=2859 end=3950 strand=+ kind=biosynthetic sec_met_domain=[...]``).
       Here the first whitespace token carries a combined ``BGC<n>_<gene>``
       identifier (PATCH-001).

    The two are disambiguated by the presence of a ``|``: a header with no
    pipe is parsed via the whitespace path so ``bgc_id``/``gene`` are still
    populated (without which ``merge_hit_xml``'s fallback match degenerates to
    ``None == None`` and silently mis-binds XML metadata — PATCH-001).
    """
    q = _text(query_id).lstrip(">")
    out: dict[str, str] = {"query_id": q}

    if "|" in q:
        parts = q.split("|")
        if parts:
            out["strain"] = parts[0]
        if len(parts) > 1 and parts[1].startswith("BGC"):
            out["bgc_id"] = parts[1]
        for part in parts[2:]:
            if "=" in part:
                k, v = part.split("=", 1)
                out[k.strip()] = v.strip()
            elif part.startswith("BGC") and "bgc_id" not in out:
                out["bgc_id"] = part
            elif part.startswith("NODE") and "node" not in out:
                out["node"] = part
            elif part.startswith("region") and "region" not in out:
                out["region"] = part
            elif "gene" not in out and ("_" in part or part.startswith("ctg")):
                out["gene"] = part
    else:
        # Space-delimited key=value convention (PATCH-001). The leading token
        # is a combined BGC+gene identifier; the remainder are key=value pairs
        # drawn from the same vocabulary as the pipe path (node=, region=,
        # start=, end=, strand=, kind=, aa=, etc.).
        tokens = q.split()
        if tokens:
            head = tokens[0]
            m = re.match(r"^(BGC\d+)_(\S+)$", head)
            if m:
                out["bgc_id"] = m.group(1)
                out["gene"] = m.group(2)
            elif head.startswith("BGC"):
                out["bgc_id"] = head
            # else: a first token that is neither BGC###_<gene> nor BGC###
            # is not a strain name — leave strain blank rather than guess
            # (honest-blank: the same standard applied to the combined-token
            # case; real manual headers always lead with BGC###_).
        for part in tokens[1:]:
            if "=" in part:
                k, v = part.split("=", 1)
                k = k.strip()
                v = v.strip()
                # First wins for keys that can legitimately recur (e.g. a bare
                # node= then a contig=); don't clobber an already-set value.
                if k not in out:
                    out[k] = v
                if k == "contig" and "node" not in out:
                    out["node"] = v
            elif part.startswith("BGC") and "bgc_id" not in out:
                out["bgc_id"] = part
            elif part.startswith("NODE") and "node" not in out:
                out["node"] = part

    if "aa" not in out:
        m = re.search(r"(?:^|[|\s])aa=(\d+)(?:[|\s]|$)", q)
        if m:
            out["aa"] = m.group(1)
    return out

def query_len_from_id(query_id: str) -> int | None:
    meta = parse_query_id(query_id)
    if meta.get("aa"):
        try:
            return int(meta["aa"])
        except Exception:
            return None
    return None

def _looks_like_hit_table_header(row: list[str]) -> bool:
    low = [c.strip().lower() for c in row]
    return "qseqid" in low or "query id" in " ".join(low) or "query_id" in low

def parse_hit_table_csv(path: str | Path) -> list[HitRecord]:
    """Parse NCBI BLASTP Hit Table CSV.

    NCBI often emits this file without a header. This parser preserves the first
    row if it looks like a real query row, which prevents the common off-by-one
    loss of the top hit.
    """
    p = Path(path)
    # BLAST-P03: mirror blastp_ingest.parse_hit_table — an HTML/QBlastInfo stub is not a Hit Table;
    # bail before csv-parsing it into junk rows (was silently yielding empty/garbage).
    _head = p.read_text(encoding="utf-8", errors="replace")[:400].lstrip()
    if _head.startswith("<") or "QBlastInfo" in _head:
        return []
    records: list[HitRecord] = []
    with p.open(newline="", encoding="utf-8-sig") as handle:
        rows = [r for r in csv.reader(handle) if r and any(c.strip() for c in r)]
    if not rows:
        return []
    header = None
    data_rows = rows
    if _looks_like_hit_table_header(rows[0]):
        header = [c.strip().lower().replace(" ", "_") for c in rows[0]]
        data_rows = rows[1:]
    for row in data_rows:
        if len(row) < 12:
            continue
        if header:
            d = {header[i]: row[i] for i in range(min(len(header), len(row)))}
            vals = [
                d.get("qseqid") or d.get("query_id") or d.get("query") or row[0],
                d.get("sseqid") or d.get("subject_id") or d.get("subject") or row[1],
                d.get("pident") or d.get("percent_identity") or d.get("pct_identity") or row[2],
                d.get("length") or d.get("align_len") or row[3],
                d.get("mismatch") or d.get("mismatches") or row[4],
                d.get("gapopen") or d.get("gap_opens") or row[5],
                d.get("qstart") or row[6],
                d.get("qend") or row[7],
                d.get("sstart") or row[8],
                d.get("send") or row[9],
                d.get("evalue") or row[10],
                d.get("bitscore") or row[11],
                d.get("ppos") or d.get("pct_positive") or (row[12] if len(row) > 12 else ""),
            ]
        else:
            # NCBI web BLASTP Hit Table CSV is often headerless. It also does
            # not reliably quote commas embedded in the query title. Sapote/Mamey
            # FASTA headers may include KCB labels such as
            # "..., complete_genome | Type: T1PKS"; a normal csv.reader then
            # splits the query id into extra fields. The BLAST columns at the
            # end are stable, so recover by treating the final 12 fields as
            # subject_id + metrics and joining everything before them back into
            # query_id. This preserves the first row and prevents false
            # zero-identity/weak-hit calls.
            if len(row) > 13:
                vals = [",".join(row[:len(row) - 12])] + row[len(row) - 12:]
            else:
                vals = row[:13] if len(row) >= 13 else row[:12] + [""]
        records.append(HitRecord(
            query_id=_text(vals[0]), subject_id=_text(vals[1]), pct_identity=_safe_float(vals[2]),
            align_len=_safe_int(vals[3]), mismatches=_safe_int(vals[4]), gap_opens=_safe_int(vals[5]),
            qstart=_safe_int(vals[6]), qend=_safe_int(vals[7]), sstart=_safe_int(vals[8]), send=_safe_int(vals[9]),
            evalue=_text(vals[10]), bitscore=_safe_float(vals[11]),
            pct_positive=_safe_float_or_none(vals[12]) if vals[12] != "" else None,
            query_len=query_len_from_id(_text(vals[0])),
        ))
    return records

def _find_text(elem: ET.Element, name: str) -> str:
    child = elem.find(f".//{{*}}{name}")
    return _text(child.text if child is not None else "")

def _find_int(elem: ET.Element, name: str) -> int | None:
    txt = _find_text(elem, name)
    if not txt:
        return None
    try:
        return int(float(txt))
    except Exception:
        return None

def parse_xml2(path: str | Path) -> dict[str, dict[str, Any]]:
    """Return query-level and hit-title metadata from NCBI XML2."""
    p = Path(path)
    # BLAST-P02: degrade to empty on an NCBI HTML/QBlastInfo error stub or a truncated/malformed
    # file instead of raising ParseError into ingest — the stub-guard the rest of the family has.
    try:
        if "QBlastInfo" in p.read_text(encoding="utf-8", errors="replace")[:400]:
            return {}
        root = ET.parse(p).getroot()
    except (ET.ParseError, OSError):
        return {}
    queries: dict[str, dict[str, Any]] = {}
    for search in root.findall(".//{*}Search"):
        qid = _find_text(search, "query-title") or _find_text(search, "query-id")
        if not qid:
            continue
        qlen = _find_int(search, "query-len")
        hits = []
        for hit in search.findall(".//{*}Hit"):
            descr = hit.find(".//{*}HitDescr")
            hsp = hit.find(".//{*}Hsp")
            if descr is None:
                continue
            hits.append({
                "subject_id": _find_text(descr, "id"),
                "subject_accession": _find_text(descr, "accession"),
                "subject_title": _find_text(descr, "title"),
                "subject_taxid": _find_text(descr, "taxid"),
                "subject_sciname": _find_text(descr, "sciname"),
                "hit_len": _find_int(hit, "len"),
                "bitscore": _safe_float(_find_text(hsp, "bit-score")) if hsp is not None else 0.0,
                "evalue": _find_text(hsp, "evalue") if hsp is not None else "",
                "identity": _find_int(hsp, "identity") if hsp is not None else None,
                "positive": _find_int(hsp, "positive") if hsp is not None else None,
                "align_len": _find_int(hsp, "align-len") if hsp is not None else None,
            })
        queries[qid] = {"query_len": qlen, "hits": hits}
    return queries

def _accession_key(subject_id: str) -> str:
    s = subject_id.strip()
    if "|" in s:
        parts = [p for p in s.split("|") if p]
        for p in parts:
            if re.match(r"^[A-Z]{1,4}_?\d", p) or re.match(r"^[A-Z]{2,}\d", p):
                return p.split(".")[0]
    return s.split(".")[0]

def merge_hit_xml(hit_records: list[HitRecord], xml_meta: dict[str, dict[str, Any]]) -> list[HitRecord]:
    if not xml_meta:
        return hit_records
    # Build exact query/title maps plus accession maps per query.
    merged: list[HitRecord] = []
    for r in hit_records:
        qmeta = xml_meta.get(r.query_id)
        if not qmeta:
            # Some NCBI hit tables may truncate query ids; try metadata-insensitive fallback.
            rmeta = parse_query_id(r.query_id)
            r_bgc = rmeta.get("bgc_id")
            r_gene = rmeta.get("gene")
            for xqid, xm in xml_meta.items():
                xmeta = parse_query_id(xqid)
                # PATCH-001: require the matched keys to be TRUTHY before
                # treating equality as a match. Without this guard, a header
                # that parse_query_id can't decompose leaves bgc_id/gene unset
                # on BOTH sides, so `None == None` is True for *every* XML
                # query and the loop silently binds the first-iterated query's
                # hits to every row.
                if not (r_bgc and r_gene):
                    break  # nothing to match on; leave qmeta None rather than mis-bind
                if xmeta.get("bgc_id") == r_bgc and xmeta.get("gene") == r_gene:
                    qmeta = xm
                    break
        qlen = qmeta.get("query_len") if qmeta else r.query_len
        hit_extra: dict[str, Any] = {}
        if qmeta:
            acc = _accession_key(r.subject_id)
            for h in qmeta.get("hits", []):
                hkeys = {_accession_key(str(h.get("subject_id", ""))), _accession_key(str(h.get("subject_accession", "")))}
                if acc in hkeys:
                    hit_extra = h
                    break
        merged.append(HitRecord(
            query_id=r.query_id, subject_id=r.subject_id, pct_identity=r.pct_identity,
            align_len=r.align_len, mismatches=r.mismatches, gap_opens=r.gap_opens,
            qstart=r.qstart, qend=r.qend, sstart=r.sstart, send=r.send, evalue=r.evalue,
            bitscore=r.bitscore, pct_positive=r.pct_positive, query_len=qlen or r.query_len,
            subject_title=hit_extra.get("subject_title", ""),
            subject_accession=hit_extra.get("subject_accession", ""),
            subject_sciname=hit_extra.get("subject_sciname", ""),
            subject_taxid=hit_extra.get("subject_taxid", ""),
            hit_len=hit_extra.get("hit_len"),
        ))
    return merged

def best_hits_by_query(records: Iterable[HitRecord]) -> dict[str, HitRecord]:
    best: dict[str, HitRecord] = {}
    for r in records:
        cur = best.get(r.query_id)
        if cur is None or (r.bitscore, r.pct_identity, r.align_len) > (cur.bitscore, cur.pct_identity, cur.align_len):
            best[r.query_id] = r
    return best
