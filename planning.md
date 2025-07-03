# Cenote-Taker3 Development Plan

> Last updated: 2025-07-03

This document sketches high-level engineering plans for upcoming improvements to the **Cenote-Taker3** command-line tool.

---

## 1️⃣  Add support for circularity declaration via TSV input

Goal: allow users to pass a tab-separated values file that declares which input records are circular so that downstream modules can treat them accordingly.

### Proposed user interface
* **New CLI option**: `--circular-tsv <path>`
  * Mutually exclusive with existing circularity options (if any).
  * Optional; if omitted, behaviour is unchanged.

### Expected TSV format
| column | description                     |
|--------|---------------------------------|
| `record_id` | Identifier matching the record IDs in the primary input (FASTA/FASTQ). |
| `is_circular` | Boolean (`1/0`, `true/false`, or `yes/no`) indicating circularity. |

*Header row is optional.*  Additional columns are ignored.

### Implementation steps
1. **Argument parsing** – extend `argparse` configuration in `src/cenote/cenotetaker3.py` to accept the new option.
2. **TSV reader** – add helper in e.g. `utils/io.py` to read and validate the file, returning a `dict[str,bool]`.
3. **Workflow integration** – where the pipeline currently infers or assumes circularity, replace with lookup from the TSV (falling back to existing logic when unspecified).
4. **Validation & error handling** – raise descriptive errors for missing/duplicate IDs or malformed booleans.
5. **Unit & integration tests** – create sample TSV and FASTA datasets in `test_data/`, update CI workflow.
6. **Documentation** – update `README.md` and `--help` message.

---

## 2️⃣  Replace bedtools dependencies with pure-Python implementations

Goal: eliminate the external `bedtools` dependency by re-implementing the sub-commands currently used inside Cenote-Taker3.

### Commands to replace
1. **intersect** (`bedtools intersect -c`-style)
   * Produce overlap counts between two BED interval sets.
2. **getfasta** (`bedtools getfasta`)
   * Extract sequence segments from a reference FASTA using BED coordinates; respect strand when requested.

### Proposed Python module
Create `src/cenote/pybed.py` exposing:

```python
from typing import Union, Dict
import pandas as pd

def intersect(a_bed: Union[str, pd.DataFrame],
              b_bed: Union[str, pd.DataFrame],
              count: bool = True) -> pd.DataFrame:
    """Return intersected intervals (optionally with counts)."""


def get_fasta(bed: Union[str, pd.DataFrame],
              fasta_path: str,
              stranded: bool = False) -> Dict[str, str]:
    """Return interval sequences keyed by identifier."""
```

Leverage `pyranges` for interval math and `pyfaidx`/`biopython` for sequence extraction, but keep fallbacks to pure pandas/SeqIO for minimal installs.

### Implementation steps
1. Audit all `bedtools` invocations (primarily in `cenote_main.sh`).
2. Capture required flags/behaviours (e.g., `-c`, strandedness, 0-/1-based coords).
3. Implement and unit-test new functions against `bedtools` outputs on representative datasets.
4. Replace shell calls with the Python equivalents inside the workflow.
5. Benchmark on large datasets; optimise with `pyranges` or `numba` if needed.
6. Remove `bedtools` from installation docs and CI requirements.

---

## 3️⃣  Pythonic pipeline refactor (eliminate `cenote_main.sh`)

Goal: convert the current mix of shell and Python scripts into a cohesive, importable Python workflow.

### Strategy
1. **Workflow audit** – list every step executed in `cenote_main.sh`, noting inputs/outputs.
2. **Module extraction** – move scripts in `src/cenote/python_modules/` into proper functions within `src/cenote/`.
3. **Shell-tool replacement** – substitute `seqkit`, `awk`, etc. with Python (pandas, BioPython, regex, etc.).
4. **Pipeline orchestration** – create a `pipeline.py` (or extend `cenotetaker3.py`) that sequentially calls these functions.
5. **Argument handling** – keep CLI thin: parse args, configure logging, invoke pipeline.
6. **Parallelism** – where `xargs`/GNU parallel were used, evaluate `concurrent.futures` or `pandas` vectorisation.
7. **Testing & validation** – compare outputs against legacy workflow on regression datasets.
8. **Documentation** – update README and example commands; remove references to `cenote_main.sh`.

Constraints: avoid creating new custom classes unless essential; favour built-ins, pandas DataFrames, and BioPython objects.

---

## 4️⃣  Additional planned features
* Automatic database update utility.
* Performance profiling & optimisation pass.
* Switch to Pydantic-based config model.

(Feel free to add more ideas!)
