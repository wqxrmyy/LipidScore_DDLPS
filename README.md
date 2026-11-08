# LipidScore_DDLPS — analysis code

Reproducibility code accompanying the manuscript:

> **A six-gene lipid metabolism score tracks copy-number clusters independently
> of grade in dedifferentiated liposarcoma**

The repository holds the analysis code and the governance locks behind the
reported statistics. It is intentionally **not** a copy of the manuscript.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23274039.svg)](https://doi.org/10.5281/zenodo.23274039)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](LICENSE)

Every published version is archived on Zenodo. **Cite the concept DOI**
`10.5281/zenodo.23274039` — it always resolves to the newest version, so a
citation stays correct when the code is updated. To pin the exact snapshot a
result came from, use that version's own DOI instead (listed on the Zenodo page).

---

## Scope: what is here, and what is deliberately not

| | Item | Note |
|---|---|---|
| ✅ | `scripts/01_analysis/` | analysis code behind the figures and tables — single-cell, TCGA-SARC, GSE cohorts, score definition |
| ✅ | `scripts/02_figures/` | one script per figure panel, plus the shared plotting library and the assembly kit |
| ✅ | `scripts/07_bulk_validation/` | two independent implementations of the GSE221492 validation + QC assertions |
| ✅ | `scripts/09_qc/` | assertions that tie the code output to the reported values |
| ✅ | `governance/` | analysis plan, metadata locks, limitations register |
| ✅ | `reproducibility/session_info/` | environment of record, plus `required_packages.txt` derived from the code |
| ✅ | `README.md`, `DATA_SOURCES.md`, `PROVENANCE.md` | documentation |
| ❌ | manuscript text (`.docx`, `.html`) | excluded — plagiarism-corpus risk |
| ❌ | figures (`.pdf`, `.tif`, `.png`) | excluded — figure-library risk |
| ❌ | `data/` | excluded — third-party; re-download, see `DATA_SOURCES.md` |
| ❌ | `results/` | excluded under the code-only release decision |
| ❌ | internal work snapshots | excluded — contain process artefacts |
| ❌ | scripts that wrote the manuscript, the caption or the limitations register | excluded — they are provenance of the document, not code behind the paper |

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
stock CRAN setup. The analysis and figure code needs more — the full set that the
scripts actually load is listed in
[`reproducibility/session_info/required_packages.txt`](reproducibility/session_info/required_packages.txt)
(derived from `library()` / `import` across the tree, because the archived
sessionInfo files record only base R).

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

## The analysis and figure code

`scripts/01_analysis/` and `scripts/02_figures/` carry the code that produced the
figures and tables. They are published as they were executed, one script per
analytic step or per panel, so a reader can trace any reported number back to the
step that computed it.

| Module | Directory | What it covers |
|---|---|---|
| Single-cell (GSE221493) | `scripts/01_analysis/scrna/` | 10x matrices → QC/normalisation → UMAP → Harmony → L1/L2 annotation → score projection → Kitagawa decomposition |
| TCGA-SARC | `scripts/01_analysis/tcga/` | copy-number clusters, expression, purity, robustness → score revisions → Harrell C and bootstrap |
| GSE cohorts | `scripts/01_analysis/validation/` | HR, C-index and direction checks for GSE30929 / GSE21122 / GSE159659 |
| Acquisition | `scripts/01_analysis/fetch/` | GEO/TCGA download and unpacking, cohort tiering |
| Score definition | `scripts/01_analysis/caliber/` | the six-gene formula, Ensembl checks, four endpoints, subgroup re-checks |
| Figures | `scripts/02_figures/` | one script per panel; shared library in `utils/`, assembly in `assemble/` |

**Paths.** The scripts read their roots from the environment instead of hard-coding
a developer machine:

| Variable | Meaning | Default |
|---|---|---|
| `LPS_DATA_ROOT` | project root holding `data/` and `results/` | `.` |
| `LPS_FIG_ROOT` | figure working root — `scripts/02_figures/` | `.` |
| `LPS_WORK_ROOT` | internal work tree referred to by a few analysis steps | `.` |
| `R_BIN_DIR` | directory holding `Rscript` (only where one script shells out) | `.` |

```sh
export LPS_DATA_ROOT=$PWD
export LPS_FIG_ROOT=$PWD/scripts/02_figures
```

**What is deliberately not here.** The 2026-09 archive contained 141 analysis
scripts and 79 figure scripts. Forty-one of the analysis scripts and twelve of
the figure helpers did not compute anything — they wrote the manuscript, the
figure captions, the limitations register or the project memory, and several of
them carried absolute paths into a personal notes directory. Those are
**provenance of the document, not code behind the paper**, and are not published.
The selection is recorded file by file in
[`scripts/SELECTION.json`](scripts/SELECTION.json), not applied by pattern.

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
  01_analysis/                                scrna/ tcga/ validation/ fetch/ caliber/
  02_figures/                                 fig1…fig5/ figS1…figS3/ utils/ assemble/
  07_bulk_validation/01_lipid_score_bulk.R    independent implementation (R)
  07_bulk_validation/01_lipid_score_gse221492.py   reference implementation (Python)
  09_qc/assert_manuscript_numbers.py          three-layer cross-validation
  SELECTION.json                              what the 2026-11 archive held and what was dropped
governance/
  analysis_plan/analysis_plan_v1.md           pre-registered analysis plan
  metadata_lock/                              frozen gene sets, sample manifests, annotations
  limitations/                                limitations register per analysis
reproducibility/session_info/                 environment of record + required_packages.txt
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
publication) and this archive. The archive DOI is the **concept DOI**, which
resolves to the newest published version:

```
Xiao, W. (2026). LipidScore-DDLPS: analysis code for a six-gene
lipid-metabolism score in dedifferentiated liposarcoma (v1.1.0).
Zenodo. https://doi.org/10.5281/zenodo.23274039
```

BibTeX:

```bibtex
@software{xiao_lipidscore_ddlps_2026,
  author    = {Xiao, Weihua},
  title     = {LipidScore-DDLPS: analysis code for a six-gene
               lipid-metabolism score in dedifferentiated liposarcoma},
  year      = {2026},
  version   = {v1.1.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23274039},
  url       = {https://doi.org/10.5281/zenodo.23274039}
}
```

`CITATION.cff` at the repository root carries the same metadata and can be read
directly by GitHub, Zenodo and reference managers.

## Contact

Corresponding author: **Weihua Xiao** — hua7708@126.com
