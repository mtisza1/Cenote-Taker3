"""Utility IO helpers for Cenote-Taker3."""
from __future__ import annotations

from pathlib import Path
from typing import Dict
import csv

__all__ = [
    "read_circular_tsv",
]


_BOOLEAN_TRUE = {"1", "true", "yes", "y", "t"}
_BOOLEAN_FALSE = {"0", "false", "no", "n", "f"}


def _str_to_bool(token: str) -> bool:
    """Convert a string token to boolean.

    Accepts common representations (case-insensitive):
        1/0, true/false, yes/no, y/n, t/f
    Raises ValueError on unrecognised tokens.
    """
    token_lc = token.strip().lower()
    if token_lc in _BOOLEAN_TRUE:
        return True
    if token_lc in _BOOLEAN_FALSE:
        return False
    raise ValueError(f"Invalid boolean value '{token}'. Expected one of "
                     f"{sorted(_BOOLEAN_TRUE | _BOOLEAN_FALSE)}")


def read_circular_tsv(path: str | Path) -> Dict[str, bool]:
    """Parse a TSV declaring circular contig records.

    Parameters
    ----------
    path:
        File system path to the TSV. The file must contain at least two columns:
        record_id and is_circular. A header row is optional. Additional columns
        are ignored.

    Returns
    -------
    dict[str, bool]
        Mapping of record identifier to circularity boolean.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If required columns are missing or duplicate record IDs are found.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    mapping: Dict[str, bool] = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        first_row_peeked = next(reader)
        # Detect header: if first cell is non-boolean and second is non-boolean
        # header assumed.
        header_present = False
        if len(first_row_peeked) >= 2:
            try:
                _str_to_bool(first_row_peeked[1])
            except ValueError:
                header_present = True
        # Reset reader to start
        handle.seek(0)
        if header_present:
            next(reader)  # skip header
        for row_num, row in enumerate(reader, start=2 if header_present else 1):
            if len(row) < 2:
                raise ValueError(f"Row {row_num} does not have 2 columns: {row}")
            record_id, circ_raw = row[0], row[1]
            try:
                circ_flag = _str_to_bool(circ_raw)
            except ValueError as exc:
                raise ValueError(f"Invalid boolean on row {row_num}: {exc}") from exc
            if record_id in mapping:
                raise ValueError(f"Duplicate record_id '{record_id}' (row {row_num})")
            mapping[record_id] = circ_flag
    return mapping
