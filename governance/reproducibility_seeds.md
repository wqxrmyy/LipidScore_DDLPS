# Random seeds used across the study

The manuscript states that "all analysis scripts, session information and
random seeds are available". This file records the seeds so that the
statement is verifiable.

## Scope: this repository vs. the rest of the study

**The analysis released in this repository (`scripts/07_bulk_validation/`,
`, `scripts/09_qc/`) is deterministic and uses no random seed at all.** The
GSE221492 paired-component statistics are computed by exact enumeration
(Wilcoxon signed-rank), by z-score averaging, and by ranking (AUC) — none of
these draws random numbers. Re-running the module reproduces the reported
values bit-for-bit without any seed being set.

The seeds below belong to the **other modules of the study** (TCGA-SARC
survival / discrimination analysis, mediation and power analyses, and the
permutation-test figures). Those scripts were developed in the project
working area and are **not** part of this repository; the seed values are
recorded here so that the corresponding analyses remain reproducible.

## How this table was produced

```sh
# scan every script in the project working area for a seed call
grep -rnoE 'set\.seed\([0-9]+\)|SEED *= *[0-9]+|default_rng\([0-9]+\)|random_state *= *[0-9]+' <project-working-area>
```

Regenerate with `_work/s103_seeds_inventory.py --scan <working-area>`.

## Seeds

| Script (project working area) | Line | Seed | Purpose (from adjacent comment) |
|---|---:|---:|---|
| `d12_p4_robustness.py` | 13 | `20260917` | — |
| `p5_figures/fig2/fig2c.R` | 12 | `20260917` | 观测 r=0.5354, P_perm=0.0001 ============================================================ |
| `p5_figures/figS2/figS2d.R` | 19 | `20260918` | — |
| `s22_harrell_A_baseline.R` | 7 | `20260917` | S22 · D18 后：用 A 评分（6 基因 z 分均值）重算 TCGA-SARC 的 Harrell C 对照：同时输出 C 评分（02_feature_scores.csv）下的 Harrell C，确认 A 与 C 的差异 |
| `s22b_delta_boot.R` | 5 | `20260918` | S22b · A vs C 的 Harrell C 差异：配对 bootstrap（同一批 50 例重抽，保持两评分相关性） 目的：回答「A 比 C 高 0.072，这个差是真的吗」；为 Fig4D 提供 ΔC 与其 95% CI |
| `s3_mediation.py` | 37 | `20260917` | — |
| `s3_power.py` | 32 | `20260917` | — |
| `s3_robustness.py` | 27 | `20260917` | — |
| `s4_addendum.R` | 14 | `1` | 用 concordancefit 的 reverse 参数直接暴露两种约定 |
| `s4_cindex.R` | 13 | `20260917` | 队列：TCGA-SARC (OS, n=50) / GSE30929 all (DFS, n=140) / GSE30929 DDLPS (DFS, n=40) ============================================================ |

## Distinct seed values

| Seed | Set in | Note |
|---|---:|---|
| `20260917` | 7 scripts | project-wide default (the date the analysis phase opened) |
| `20260918` | 2 scripts | paired bootstrap / permutation figure using the later data cut |
| `1` | 1 script | only this script; value recorded as found, meaning not inferred |

Scanning also surfaced `set.seed(2024)` inside a plain-text note
(`_cprov.txt`). That is a quotation, not an executable seed call, so it is
excluded from the table and the count above.

Total: **10 seed calls across 10 scripts**, distinct values: `1`, `20260917`, `20260918`.
