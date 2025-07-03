"""Genome-related helper functions refactored from legacy shell/Python scripts.

This module centralises utilities that operate on contig sequences such as
terminal repeat detection and sequence trimming. Functions were migrated from
`python_modules/terminal_repeats.py` to allow clean, importable usage throughout
Cenote-Taker3's new Python pipeline.
"""
from __future__ import annotations

import os
import re
import string
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

__all__ = [
    "fetch_dtr",
    "reverse_complement",
    "fetch_itr",
    "assess_terminal_seqs",
    "get_min_length_contigs",
    "get_min_hallmark_contigs",
]


# ---------------------------------------------------------------------------
# Utilities copied (with min. tweaks) from CheckV via original terminal_repeats
# ---------------------------------------------------------------------------

def fetch_dtr(fullseq: str, min_length: int = 20) -> str:
    """Return direct terminal repeat (DTR) sequence or empty string.

    Heuristic from CheckV: search for the first *min_length* bases at the start
    of *fullseq* re-occurring in the second half of the sequence and extending
    to the contig end.
    """
    startseq = fullseq[0:min_length]
    # Find all occurrences of *startseq*; keep those in 2nd half of sequence.
    matches = [m.start() for m in re.finditer(f"(?={re.escape(startseq)})", fullseq)]
    matches = [_ for _ in matches if _ >= len(fullseq) / 2]
    for matchpos in matches:
        endseq = fullseq[matchpos:]
        if fullseq[: len(endseq)] == endseq:
            return endseq  # confirmed DTR
    return ""


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence (A/C/T/G only)."""
    trans = str.maketrans("ACTGactg", "TGACtgac")
    return seq.translate(trans)[::-1]


def fetch_itr(seq: str, min_len: int = 20, max_len: int = 1000) -> str:
    """Return inverted terminal repeat (ITR) sequence or empty string."""
    rev = reverse_complement(seq)
    if seq[:min_len] == rev[:min_len]:
        i = min_len + 1
        while i <= max_len and seq[:i] == rev[:i]:
            i += 1
        return seq[: i - 1]
    return ""


# ---------------------------------------------------------------------------
# High-level workflow function
# ---------------------------------------------------------------------------

def assess_terminal_seqs(
    fasta_path: str | Path,
    hallmark_table: str | Path,
    length_circ: int,
    length_lin: int,
    hall_circ: int,
    hall_lin: int,
    out_dir: str | Path,
    wrap: bool = True,
    max_length: int | float = 1_000_000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Analyse DTR/ITR for each contig and write summary tables.

    Parameters
    ----------
    fasta_path:
        FASTA file of contigs.
    hallmark_table:
        Path to `hallmarks_per_orig_contigs.tsv` produced earlier in the pipeline.
    length_circ / length_lin:
        Minimum contig length thresholds (nts) for circular / linear contigs.
    hall_circ / hall_lin:
        Minimum total hallmark counts required.
    out_dir:
        Directory to write output tables and trimmed FASTA.
    wrap:
        If True, trim detected DTR from circular contigs.
    max_length:
        Skip DTR trimming for contigs longer than this.

    Returns
    -------
    terminal_df, keep_df : tuple of pandas.DataFrame
        `terminal_df` - summary for every contig.
        `keep_df` - subset passing thresholds.
    """
    fasta_path = Path(fasta_path)
    hallmark_table = Path(hallmark_table)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = out_dir / "trimmed_TRs_hallmark_contigs.fasta"
    if output_fasta.exists():
        output_fasta.unlink()

    # Build summary list
    terminal_rows: List[Tuple[str, int, int, str, str]] = []
    with fasta_path.open() as handle, output_fasta.open("a") as fasta_out:
        for record in SeqIO.parse(handle, "fasta"):
            dtr_seq = fetch_dtr(str(record.seq))
            if not dtr_seq or len(record) > max_length:
                dtr_seq = "NA"

            itr_seq = fetch_itr(str(record.seq)) or "NA"

            if dtr_seq != "NA" and wrap:
                trimmed_seq = record.seq[:-len(dtr_seq)]
            else:
                trimmed_seq = record.seq

            # Write to output FASTA
            print(f">{record.id}", file=fasta_out)
            print(trimmed_seq, file=fasta_out)

            terminal_rows.append(
                (
                    record.id,
                    len(record.seq),
                    len(trimmed_seq),
                    dtr_seq,
                    itr_seq,
                )
            )

    terminal_df = pd.DataFrame(
        terminal_rows,
        columns=[
            "contig",
            "in_length_contig",
            "out_length_contig",
            "dtr_seq",
            "itr_seq",
        ],
    )
    terminal_df.to_csv(out_dir / "hallmark_contigs_terminal_repeat_summary.tsv", sep="\t", index=False)

    hall_df = pd.read_csv(hallmark_table, sep="\t")
    full_df = hall_df.merge(terminal_df, on="contig", how="left")

    keep_df = full_df[
        (
            full_df["dtr_seq"].notna()
            & (full_df["total_hit_count"] >= hall_circ)
            & (full_df["in_length_contig"] >= length_circ)
        )
        | (
            (full_df["total_hit_count"] >= hall_lin)
            & (full_df["in_length_contig"] >= length_lin)
        )
    ][[
        "contig",
        "in_length_contig",
        "out_length_contig",
        "dtr_seq",
        "itr_seq",
    ]]

    keep_df.to_csv(out_dir / "threshold_contigs_terminal_repeat_summary.tsv", sep="\t", index=False)
    keep_df["contig"].to_csv(out_dir / "contigs_over_threshold.txt", sep="\t", index=False, header=False)

    return terminal_df, keep_df

# ---------------------------------------------------------------------------
# Length filtering and header renaming (replacement for seqkit pipeline)
# ---------------------------------------------------------------------------

def get_min_length_contigs(
    input_fasta: str | Path,
    output_fasta: str | Path,
    min_length: int,
    run_title: str,
) -> Tuple[Dict[str, str], Path]:
    """Filter contigs by *min_length* and rename headers.

    Each retained record receives a new ID formatted as
    `{run_title}_{n}` where `n` is an incrementing integer (1-based order of
    appearance). The original ID is preserved after a space in the description
    line, matching legacy behaviour.

    Parameters
    ----------
    input_fasta:
        Path to original contigs FASTA.
    output_fasta:
        Destination FASTA path for filtered, renamed sequences.
    min_length:
        Minimum length in nucleotides.
    run_title:
        Prefix for new contig identifiers.

    Returns
    -------
    tuple
        (mapping, output_path)
        * mapping : dict[str, str] – new_id ➜ original_id mapping.
        * output_path : pathlib.Path – path to the written FASTA.
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    mapping: Dict[str, str] = {}

    with input_fasta.open() as in_handle, output_fasta.open("w") as out_handle:
        for idx, record in enumerate(SeqIO.parse(in_handle, "fasta"), start=1):
            if len(record.seq) < min_length:
                continue
            new_id = f"{run_title}_{idx}"
            original_id = record.id
            record.id = new_id
            record.description = f"{original_id}"
            SeqIO.write(record, out_handle, "fasta")
            mapping[new_id] = original_id
    return mapping, output_fasta


# ---------------------------------------------------------------------------
# Extract subset FASTA for contigs meeting hallmark criteria
# ---------------------------------------------------------------------------

def get_min_hallmark_contigs(
    input_fasta: str | Path,
    ids: List[str] | Path,
    output_fasta: str | Path,
) -> Path:
    """Write a FASTA containing only records whose IDs are in *ids*.

    Parameters
    ----------
    input_fasta:
        Source FASTA file.
    ids:
        Either an iterable/list of record IDs to keep **or** a path to a
        newline-delimited text file containing one ID per line (as produced by
        earlier steps like `contigs_to_keep.txt`).
    output_fasta:
        Destination FASTA path.

    Returns
    -------
    pathlib.Path
        Path to the written subset FASTA.
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    # Load ID list
    if isinstance(ids, (str, Path)):
        with Path(ids).open() as handle:
            keep_ids = {line.strip().split()[0] for line in handle if line.strip()}
    else:
        keep_ids = set(ids)

    with input_fasta.open() as in_handle, output_fasta.open("w") as out_handle:
        written = 0
        for record in SeqIO.parse(in_handle, "fasta"):
            if record.id in keep_ids:
                SeqIO.write(record, out_handle, "fasta")
                written += 1
        if written == 0:
            raise ValueError("No contigs matched provided ID list; output FASTA not created.")
    return output_fasta

