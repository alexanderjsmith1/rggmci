"""CSV/TSV formula-injection neutraliser (v9.7.409, DEEP_AUDIT2 export-injection;
v9.7.410 CLAUDE_410_csv_writer_coverage: full writer coverage + csv.writer analogue).

BC_409 hostile-gates H8 neutralised spreadsheet formulas for the .xlsx workbooks only
(via openpyxl cell typing in mamey/xlsx_determinism.py). The plain csv.DictWriter / TSV
writers are a separate, unguarded output surface: a gene/product/def cell copied verbatim
from antiSMASH that begins with '=', '+', '-', '@' (or a tab/CR) is executed as a formula
when the CSV/TSV is opened in Excel or LibreOffice (DDE / HYPERLINK exfiltration).

This module is the CSV analogue of H8: it prefixes any such leading character with a single
quote — the standard CSV-injection neutralisation (OWASP CSV Injection). ONLY cells whose
first character is dangerous change; every other cell is byte-identical, so machine readers
and legitimate content (e.g. "T1PKS; NRPS") are unaffected.

THE RULE (single source of truth; the HTML-widget ``exportCsv`` JavaScript in
mamey/widget_deliverable.py is generated from the constants below so both surfaces apply the
identical rule):

* A **string** cell whose first character is one of ``= + - @ TAB CR`` is prefixed with ``'``.
* Exemption for the sign leaders ``+`` / ``-`` only: a lone sign (the strand column ``+``/``-``)
  or a plain signed number (``-1.5``, ``+3``, ``-2e-4``; regex ``CSV_NUMERIC_TOKEN``) passes
  through unchanged. A spreadsheet reads those as a number or as inert text, never as a
  formula, and prefixing them would corrupt numeric round-trips (``float("'-1.5")`` fails).
* Non-string cells (int, float, None, bool) are never touched.

Coverage (v9.7.410): every ``csv.DictWriter`` / ``csv.writer`` site under ``mamey/`` and
``deliverable_tools/`` that emits user- or tool-derived text routes through ``SafeDictWriter`` /
``SafeWriter``. Deliberately left on the plain writers: raw BLAST result captures that must stay
byte-faithful to the external tool (``blastp_online.py``, ``blastp_ebi.py``,
``ebi_xml_to_outfmt10.py``) and the numeric-only phase-timing table (``timing.py``).
"""
from __future__ import annotations

import csv
import re

# The characters a spreadsheet treats as the start of a formula / command. Tab and CR are
# included because a leading whitespace control can smuggle a formula past a naive eye.
_CSV_INJECTION_LEADERS = ("=", "+", "-", "@", "\t", "\r")

# Plain signed number: the only '+'/'-'-led strings a spreadsheet will NOT evaluate as a
# formula. Kept to a syntax that is valid in both Python ``re`` and JavaScript ``RegExp`` so
# the widget JS can embed the very same pattern string.
CSV_NUMERIC_TOKEN = r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$"
# re.ASCII: JavaScript ``\d`` is ASCII-only, so Python must be too or ``-١`` (Arabic-Indic digit)
# would be exempted here and prefixed in the widget — the two surfaces must agree cell-for-cell.
_NUMERIC_TOKEN_RE = re.compile(CSV_NUMERIC_TOKEN, re.ASCII)

# JavaScript character-class body derived from the Python policy (used by widget_deliverable
# to build the client-side ``exportCsv`` guard from this module rather than from a hand-copied
# literal). Evaluates to ``=+\-@\t\r``.
CSV_LEADER_CLASS_JS = "".join(
    "\\t" if leader == "\t" else "\\r" if leader == "\r"
    else "\\-" if leader == "-" else leader
    for leader in _CSV_INJECTION_LEADERS
)


def csv_safe_cell(value):
    """Neutralise one cell. A string that would be read as a formula gets a leading "'";
    everything else (numbers, None, safe strings, lone strand signs, plain signed numbers)
    passes through unchanged."""
    if not (isinstance(value, str) and value and value[0] in _CSV_INJECTION_LEADERS):
        return value
    if value[0] in "+-" and (len(value) == 1 or _NUMERIC_TOKEN_RE.match(value)):
        return value
    return "'" + value


class SafeDictWriter(csv.DictWriter):
    """Drop-in replacement for csv.DictWriter that neutralises formula-leading string cells
    at write time. csv.DictWriter.writerows() does NOT delegate to writerow(), so both are
    overridden. Header/field handling is inherited unchanged."""

    def writerow(self, rowdict):
        return super().writerow({k: csv_safe_cell(v) for k, v in rowdict.items()})

    def writerows(self, rowdicts):
        # A LIST, not a generator expression. csv.DictWriter.writerows accepts any iterable, but
        # callers and tests that wrap or monkeypatch this method index into the argument — the
        # shipped tests/test_legacy_feature_gate_atomic_write_bc2_401.py does `rows[0]` to simulate
        # a mid-write interruption, and a generator turns that crash-safety test into a TypeError.
        # Materialising costs one list of rows the caller already holds in memory.
        return super().writerows(
            [{k: csv_safe_cell(v) for k, v in rd.items()} for rd in rowdicts]
        )


class SafeWriter:
    """Drop-in replacement for ``csv.writer(...)`` (positional-row writer) applying the same
    per-cell rule. ``csv.writer`` returns a C-level object that cannot be subclassed, so this
    is a thin wrapper; ``dialect`` and every ``**fmtparams`` are forwarded unchanged."""

    def __init__(self, csvfile, dialect="excel", **fmtparams):
        self._w = csv.writer(csvfile, dialect, **fmtparams)

    @property
    def dialect(self):
        return self._w.dialect

    def writerow(self, row):
        return self._w.writerow([csv_safe_cell(v) for v in row])

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)
