# LipidScore_DDLPS — analysis code

Reproducibility code accompanying the manuscript:

> **A six-gene lipid metabolism score tracks copy-number clusters independently
> of grade in dedifferentiated liposarcoma**

The repository holds the analysis code and the governance locks behind the
reported statistics. It is intentionally **not** a copy of the manuscript.

---

## Scope: what is here, and what is deliberately not

| | Item | Note |
|---|---|---|
| ✅ | `scripts/` | analysis code — two independent implementations + QC assertions |
| ✅ | `governance/` | analysis plan, metadata locks, limitations register |
| ✅ | `reproducibility/session_info/` | environment of record |
| ✅ | `README.md`, `DATA_SOURCES.md`, `PROVENANCE.md` | documentation |
| ❌ | manuscript text (`.docx`, `.html`) | excluded — plagiarism-corpus risk |
| ❌ | figures (`.pdf`, `.tif`, `.png`) | excluded — figure-library risk |
| ❌ | `data/` | excluded — third-party; re-download, see `DATA_SOURCES.md` |
| ❌ | `results/` | excluded under the code-only release decision |
| ❌ | internal work snapshots | excluded — contain process artefacts |

The exclusion list is not only a `.gitignore` convention:
`scripts/00_repo_hygiene/check_hygiene.py` asserts against the **git index**, so
a stray `git add -f` or a renamed pattern is caught rather than silently
committed. Run it before every push.

> **Consequence worth stating plainly.** Because `data/` and `results/` are not
> shipped, **this repository cannot be executed end-to-end as-is**. It is
> auditable as it stands — you can read every step of the analysis and check the
> assertions — but reproducing the numbers requires first downloading the source
> data (`DATA_SOURCES.md`) and then running the pipeline (`run_all.sh`). The
> result tables behind every number in the manuscript are available from the
> corresponding author on request.

---

## Data

All datasets are public; **no new sequence data were generated**. Accessions and
per-file provenance (including recorded byte sizes) are in
[`DATA_SOURCES.md`](DATA_SOURCES.md).

Primary cohort for the validation module in `scripts/07_bulk_validation/`:

| | |
|---|---|
| Accession | **GSE221492** (bulk RNA-seq atlas of liposarcoma) |
| Samples | 53 samples / 41 patients |
| File used | `GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz` |

Place the downloaded file at:

```
data/GSE221492_bulk/raw/GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz
```

The pipeline locates the repository root via the `LPS_ROOT` environment variable
and falls back to inferring it from the script location, so the working
directory does not matter.

---

## Requirements

Verified on the environment of record (`reproducibility/session_info/`):

| | Version |
|---|---|
| R | 4.6.0 (2026-04-24 ucrt), x86_64-w64-mingw32 |
| R packages | `data.table`, `ggplot2`, `dplyr`, `tidyr` |
| Python | 3.13 |
| Python packages | `numpy`, `scipy` |

The R implementation deliberately uses only those four packages, so it runs on a
stock CRAN setup.

---

## Quick start

```sh
git clone https://github.com/wqxrmyy/LipidScore_DDLPS.git
cd LipidScore_DDLPS

# 1) fetch the source data (see DATA_SOURCES.md for the accessions)
#    -> data/GSE221492_bulk/raw/GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz

# 2) run both implementations and the cross-check
sh run_all.sh
```

> **Branch of record: `master`.** The code described here is the repository
> default branch, so a plain `git clone` is all that is needed — no `-b` flag.
>
> **A note on this repository's history.** An earlier, unrelated project — a
> different manuscript and pipeline, with no shared commit ancestry —
> previously occupied this repository. It has been withdrawn and is not
> reachable from any ref.

`run_all.sh` runs, in order: the R implementation, the Python implementation,
the cross-validation assertions, and the repository hygiene check. Any failing
assertion exits non-zero.

To run a single step:

```sh
Rscript         scripts/07_bulk_validation/01_lipid_score_bulk.R
python          scripts/07_bulk_validation/01_lipid_score_gse221492.py
python          scripts/09_qc/assert_manuscript_numbers.py
python          scripts/00_repo_hygiene/check_hygiene.py
```

---

## What is reproduced

`scripts/07_bulk_validation/` contains **two independent implementations** of the
same analysis. They are kept side by side on purpose: the Python implementation
is the origin of the numbers in the manuscript, the R implementation is an
independent re-derivation, and `scripts/09_qc/assert_manuscript_numbers.py`
requires them to agree digit-for-digit before either is trusted.

Reported in the manuscript vs. produced by the code (GSE221492, paired
components):

| Quantity | Manuscript | Source |
|---|---|---|
| Patients with strictly paired components | 8 | `*_summary.csv` |
| Paired mean Δ | **−0.587** | `*_paired_delta.csv` |
| Pairs with DD < WD | 6 of 8 | `*_summary.csv` |
| Wilcoxon signed-rank P | **0.078125** | `*_summary.csv` |
| Unpaired DD mean | −0.305 | `*_summary.csv` |
| Unpaired WD mean | +0.164 | `*_summary.csv` |
| Unpaired Δ / Welch P | −0.469 / 0.025 | `*_summary.csv` |
| AUC, DD vs WD / P | 0.249 / 0.013 | `*_summary.csv` |

Note that the manuscript reports the **mean** paired Δ. The **median** paired Δ
(−0.689) is also produced but is not a reported value; the two were once
conflated in the R implementation, which is why both are now emitted under
explicit labels.

---

## Cross-validation

`scripts/09_qc/assert_manuscript_numbers.py` checks three layers:

- **L1 — implementations agree**: per-patient Δ, n, DD<WD count, P, mean and
  median Δ must match between R and Python.
- **L2 — reported values**: every number in the table above is asserted against
  both implementations.
- **L3 — manuscript text** (optional, `--manuscript PATH`): the manuscript
  really prints those numbers — guarding against the failure mode where the code
  is right and the text is wrong.

The two implementations agree to within `4.5e-5` on Δ, which is the rounding
floor of the JSON object (values stored to 4 decimals).

---

## Two implementation pitfalls, documented

Both were found and fixed during preparation of this repository; they are
recorded here because either one silently shifts the reported Δ in its third
decimal.

**1. Library size must be computed before collapsing duplicate gene symbols.**
Collapsing duplicate `hgnc_symbol` rows discards 12,873 of 45,473 rows (28.3%)
of counts and shrinks the size factor by 0.47%, inflating CPM uniformly. The
six score genes each have exactly one row, so gene values are unaffected — the
denominator is the *only* downstream effect. It changes the paired mean Δ from
−0.586858 (correct) to −0.586407.

**2. `DDLPS_sWD` is not `DDLPS_WD`.** The strict paired set requires a patient to
have both a genuine `DDLPS_WD` and a `DDLPS_DD` sample → n = 8. Folding
superficial well-differentiated components into the WD group yields n = 10, a
result **the manuscript never reports**. The R implementation previously
labelled the n=10 analysis "Main" and the n=8 analysis "Sensitivity"; the
priority is now reversed to match the manuscript.

---

## Repository layout

```
scripts/
  00_repo_hygiene/check_hygiene.py            no manuscript/figures/snapshots tracked
  07_bulk_validation/01_lipid_score_bulk.R    independent implementation (R)
  07_bulk_validation/01_lipid_score_gse221492.py   reference implementation (Python)
  09_qc/assert_manuscript_numbers.py          three-layer cross-validation
governance/
  analysis_plan/analysis_plan_v1.md           pre-registered analysis plan
  metadata_lock/                              frozen gene sets, sample manifests, annotations
  limitations/                                limitations register per analysis
reproducibility/session_info/                 environment of record
run_all.sh                                    one-shot pipeline
DATA_SOURCES.md                               dataset accessions and per-file provenance
PROVENANCE.md                                 how the commit history was reconstructed
```

---

## History

The commit history is **reconstructed from file timestamps** — the analysis
predates any commit, so the history orders artefacts rather than recording the
work as it happened. This is disclosed in full, including what could not be
recovered, in **[`PROVENANCE.md`](PROVENANCE.md)**. Please read it before drawing
conclusions from commit dates.

---

## A note on language

The analysis code carries its comments in **Chinese**, the working language of
the project — this includes `scripts/07_bulk_validation/` and `scripts/09_qc/`.
The documentation you are reading, `DATA_SOURCES.md`, `PROVENANCE.md`, all
assertion output and every commit message are in English.

If you are reviewing the code and a comment is unclear, please open an issue.
The analysis logic itself is language-independent, and the QC script prints its
assertions in English so the numbers can be checked without reading the comments.

---

## Citation

If you use this code, please cite the manuscript (full citation to be added upon
publication) and this repository.

## Contact

Corresponding author: **Weihua Xiao** — hua2577@126.com
