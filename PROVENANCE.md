# Provenance and history disclosure

## Why this file exists

This repository is the **release snapshot** prepared for manuscript submission.
The analysis behind it started on **2026-05-30** and ran through **2026-09-25**,
but **no git commit was made during that entire period** — `git init` had been
run, the object store was empty.

The commit history you see here is therefore **reconstructed from file-system
modification times (`mtime`)**. It orders the artefacts faithfully. It is *not*
a day-by-day log of the work, and it should not be read as one.

This is disclosed explicitly because **a reconstructed history that is not
disclosed is indistinguishable from a fabricated one**. If you are a reviewer or
an editor and you need to judge the authenticity of this codebase, read
"Verify it yourself" at the bottom — the check is reproducible.

## Evidence

All dates below are `mtime` values of files still present on the authoring
machine, measured on **2026-09-25**.

### 1. Analysis scripts across the whole project

```sh
find <project-root> -type f \( -name '*.py' -o -name '*.R' \) \
     -printf '%TY-%Tm-%Td\n' | sort | uniq -c
```

→ **356 scripts**, spanning 2026-05-30 → 2026-09-25:

| Date | scripts | | Date | scripts |
|---|---:|---|---|---:|
| 2026-05-30 | 9 | | 2026-07-28 | 12 |
| 2026-05-31 | 1 | | 2026-08-01 | 6 |
| 2026-06-01 | 5 | | 2026-08-05 | 4 |
| 2026-06-06 | 4 | | 2026-08-06 | 4 |
| 2026-06-21 | 78 | | 2026-09-06 | 3 |
| 2026-06-22 | 3 | | 2026-09-20 | 186 |
| 2026-06-23 | 10 | | 2026-09-22 | 3 |
| 2026-07-04 | 2 | | 2026-09-25 | 5 |
| 2026-07-18 | 16 | | | |
| 2026-07-27 | 5 | | **total** | **356** |

The earliest surviving script is dated **2026-05-30**. Nothing on this machine
carries an earlier analysis timestamp.

### 2. Files inside *this* repository (internal snapshots excluded)

| Date | files | note |
|---|---:|---|
| 2026-09-20 | 42 | data + metadata locks + analysis plan + first result objects |
| 2026-09-21 | 12 | **all** are `.Rproj.user/` session state (ignored, not committed) |
| 2026-09-25 | 17 | the two-implementation refactor and the library-size fix |

Note the shape of this: the repository content was **imported in one pass on
2026-09-20**, not grown commit-by-commit. That import is why the reconstructed
history has only a handful of dates.

### 3. Upstream work area (not published)

An internal working area holds the 597 analysis scripts produced between
**2026-09-17 13:25** and **2026-09-24 23:03**. It is deliberately excluded from
this repository (it contains internal process artefacts), so it is cited here as
evidence only and cannot be offered as a public artefact.

## How commits were dated

Each commit is stamped with the `mtime` of the artefacts it contains — no date
in this history is earlier than the file evidence supporting it. Where the
evidence supports only a single bulk date, a single bulk commit was made rather
than splitting it into invented daily commits.

Concretely: the 2026-09-20 import is split by artefact class (repository
scaffolding → governance locks → analysis code), all stamped 2026-09-20; the
2026-09-25 work is stamped 2026-09-25.

## What could NOT be recovered

- **Ordering finer than one day for 2026-09-20.** Several distinct pieces of
  work (metadata locks, analysis plan, first results) were imported the same
  evening. Their internal order is inferred from `mtime` to the minute but is
  not independently corroborated.
- **Anything before 2026-05-30.** No file on this machine carries an earlier
  timestamp. If analysis was performed before that date, **it is not evidenced
  here and no commit claims it**.
- **Intermediate states.** Because nothing was committed during development, no
  intermediate version of any script survives except the fix backup retained
  outside this repository. The single published draft of each script is its
  final state.

## Verify it yourself

```sh
git log --format='%h  %ad  %s' --date=short     # read the reconstructed timeline
git rev-list --all --count                       # number of commits
python scripts/00_repo_hygiene/check_hygiene.py  # no manuscript/figures/snapshots tracked
python scripts/09_qc/assert_manuscript_numbers.py  # every reported number is reproducible
```

If your reading of this repository's history matters for your decision, the
last two commands are the ones that carry the actual weight: the first only
tells you *when* the files were stamped, the last two tell you whether the
science in them holds up.
