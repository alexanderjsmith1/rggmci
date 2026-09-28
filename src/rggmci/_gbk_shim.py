"""Minimal GenBank parser shim — replaces Biopython SeqIO for antiSMASH region GBKs."""
from __future__ import annotations
import re
from dataclasses import dataclass, field


@dataclass
class _Location:
    start: int = 0
    end: int = 0
    strand: int = 1  # 1 = forward, -1 = complement

    def extract(self, seq):
        """Extract subsequence (mimics Biopython)."""
        subseq = seq[self.start:self.end]
        if self.strand == -1:
            comp = str.maketrans('ATCGatcg', 'TAGCtagc')
            subseq = subseq.translate(comp)[::-1]
        return subseq

@dataclass
class _Feature:
    type: str = ""
    qualifiers: dict[str, list[str]] = field(default_factory=dict)
    location: _Location = field(default_factory=_Location)

    def extract(self, seq):
        """Extract nucleotide subsequence (mimics Biopython feature.extract)."""
        return self.location.extract(seq)

@dataclass
class _Record:
    id: str = ""
    seq: str = ""
    features: list[_Feature] = field(default_factory=list)
    annotations: dict = field(default_factory=dict)
    # Biopython's record.name (the LOCUS name) and record.description (the DEFINITION line). Callers read them:
    # parsers._replicon_intake_key orders chromosome before plasmid records from the description.
    name: str = ""
    description: str = ""


def _parse_qualifiers(lines: list[str]) -> dict[str, list[str]]:
    """Parse /key="value" qualifiers, handling multi-line values."""
    qualifiers: dict[str, list[str]] = {}
    # Join continuation lines: qualifier lines start with /
    joined: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('/'):
            joined.append(stripped)
        elif joined:
            joined[-1] += ' ' + stripped
    
    for entry in joined:
        m = re.match(r'/(\w+)="(.*)"$', entry, re.DOTALL)
        if m:
            key, val = m.group(1), m.group(2)
            qualifiers.setdefault(key, []).append(val)
        else:
            # Flag qualifier without value
            m2 = re.match(r'/(\w+)$', entry)
            if m2:
                qualifiers.setdefault(m2.group(1), []).append("")
    # PARSE-01: continuation lines are joined with a literal space above, which corrupts
    # sequence-valued qualifiers whose multi-line values must be concatenated without gaps
    # (notably /translation). Collapse all interior whitespace on the translation qualifier
    # only; descriptive qualifiers (/product, /note, ...) legitimately contain spaces and are
    # left untouched.
    if "translation" in qualifiers:
        qualifiers["translation"] = ["".join(v.split()) for v in qualifiers["translation"]]
    return qualifiers


def _parse_features(text: str) -> list[_Feature]:
    """Parse FEATURES table from GenBank text."""
    features = []
    # Split into feature blocks: each starts with "     TYPE" (5 spaces + key)
    blocks = re.split(r'\n(?=     [A-Za-z_])', '\n' + text)
    
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        header = lines[0]
        
        # Parse type and location from first line
        hm = re.match(r'(\S+)\s+(.*)', header)
        if not hm:
            continue
        ftype = hm.group(1)
        loc_str = hm.group(2).strip()
        
        # Parse coordinates and strand
        is_complement = 'complement' in loc_str.lower()
        strand = -1 if is_complement else 1
        nums = re.findall(r'(\d+)\.\.(\d+)', loc_str)
        if nums:
            start = int(nums[0][0]) - 1  # 0-based like Biopython
            end = int(nums[-1][1])
        else:
            singles = re.findall(r'(\d+)', loc_str)
            if singles:
                start, end = int(singles[0]) - 1, int(singles[-1])
            else:
                start, end = 0, 0
        # v9.7.410 hostile audit: the regex above takes the outermost numbers of whatever the
        # location string holds, so `9..1`, `complement(join(9..1,abc))` or a negative start
        # produced an INVERTED span (start > end) that downstream `seq[start:end]` slices turn
        # into silent empties and length maths turn negative. Normalise to a well-formed,
        # possibly zero-length span; a well-formed antiSMASH GBK is untouched by this.
        start = max(0, start)
        if end < start:
            end = start

        qualifiers = _parse_qualifiers(lines[1:])
        features.append(_Feature(type=ftype, qualifiers=qualifiers,
                                 location=_Location(start=start, end=end, strand=strand)))
    return features


def parse_genbank_text(text: str) -> list[_Record]:
    """Parse GenBank records from text. Returns list of _Record."""
    records = []
    entries = re.split(r'\n//\s*', text)
    
    for entry in entries:
        entry = entry.strip()
        if not entry or 'LOCUS' not in entry[:500]:
            continue
        
        # ID
        ver_m = re.search(r'^VERSION\s+(\S+)', entry, re.MULTILINE)
        loc_m = re.search(r'^LOCUS\s+(\S+)', entry, re.MULTILINE)
        rec_id = ver_m.group(1) if ver_m else (loc_m.group(1) if loc_m else 'unknown')
        
        # Sequence
        origin_m = re.search(r'^ORIGIN\s*\n(.*)', entry, re.MULTILINE | re.DOTALL)
        seq = re.sub(r'[^a-zA-Z]', '', origin_m.group(1)) if origin_m else ''
        if not seq:
            len_m = re.search(r'^LOCUS\s+\S+\s+(\d+)\s+bp', entry, re.MULTILINE)
            if len_m:
                seq = 'N' * int(len_m.group(1))
        
        # Features
        feat_m = re.search(r'^FEATURES\s+Location/Qualifiers\n(.*?)(?=^ORIGIN|^CONTIG|^BASE COUNT|\Z)',
                           entry, re.MULTILINE | re.DOTALL)
        features = _parse_features(feat_m.group(1)) if feat_m else []
        
        # Annotations (topology from LOCUS line)
        annot: dict = {}
        topo_m = re.search(r'\b(linear|circular)\b', entry.split('\n', 1)[0], re.IGNORECASE)
        if topo_m:
            annot['topology'] = topo_m.group(1).lower()

        # DEFINITION, as Biopython reads it: continuation lines joined by one space, one final period dropped.
        def_m = re.search(r'^DEFINITION  (.*(?:\n {12}.*)*)', entry, re.MULTILINE)
        description = " ".join(def_m.group(1).split()) if def_m else ""
        if description.endswith("."):
            description = description[:-1]

        records.append(_Record(id=rec_id, seq=seq, features=features, annotations=annot,
                               name=loc_m.group(1) if loc_m else rec_id, description=description))
    return records


# --- SeqIO-compatible adapter (Bio-free fallback for tools) -------------------
# Lets a tool do:  try: from Bio import SeqIO
#                  except ImportError: from mamey._gbk_shim import SeqIO
# and keep calling SeqIO.parse(handle, "genbank") when biopython is absent.
class _ShimSeqIO:
    @staticmethod
    def parse(handle, fmt="genbank"):
        if hasattr(handle, "read"):
            text = handle.read()
        else:
            with open(handle, encoding="utf-8") as _fh:
                text = _fh.read()
        return iter(parse_genbank_text(text))


SeqIO = _ShimSeqIO()
