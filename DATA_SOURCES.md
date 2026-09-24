# Data sources

All datasets analysed are public. **No new sequence data were generated in this
study.** This file records the accessions and the per-file provenance needed to
rebuild the `data/` directory, which is excluded from version control.

Downloaded artefacts are placed under `data/<accession>/raw/`. The recorded byte
sizes below are the sizes of the files actually analysed, so they can be used to
confirm a re-download is the same object.

---

## Primary cohort — validation module

### GSE221492 — bulk RNA-seq atlas of liposarcoma

Used by `scripts/07_bulk_validation/`.

| | |
|---|---|
| Accession | GSE221492 |
| Portal | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221492> |
| Composition | 53 samples / 41 patients — 17 dedifferentiated and 17 well-differentiated DDLPS components, 2 superficial well-differentiated components, 6 WDLPS, 5 lipomas, 5 peritumoural fat, 1 hibernoma |
| Strictly paired patients | 8 (the paired analysis set) |

Required files:

| File | Bytes | Used for |
|---|---:|---|
| `GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz` | 2,556,624 | **all reported statistics** |
| `GSE221492_bulk_LPS_atlas_log2_FPKM_counts_matrix.tsv.gz` | 2,728,365 | cross-checking only |

Expected layout:

```
data/GSE221492_bulk/raw/GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz
data/GSE221492_bulk/raw/GSE221492_bulk_LPS_atlas_log2_FPKM_counts_matrix.tsv.gz
```

The matrix is tab-separated with the sample identifiers on the first line
(column 3 onward) and, per gene row, column 1 = Ensembl gene ID,
column 2 = `hgnc_symbol`.

> **Retrieval note.** These files are distributed as GEO supplementary files.
> The direct FTP host was not reachable from the authoring machine, so the
> download route is not asserted here; retrieve them through the GEO accession
> page above and verify against the recorded byte sizes.

---

## Cohorts referenced in the manuscript

Listed for completeness. Only GSE221492 feeds the module released in this
repository; the others belong to the wider study and were analysed in the
project's own workspace.

| Accession | Role in the study | n |
|---|---|---|
| TCGA-SARC | discovery cohort, copy-number clusters | 50 DDLPS patients |
| GSE21122 | discovery cohort | 158 (149 sarcoma + 9 normal fat) |
| GSE159659 | discovery cohort, matched DD / WD / normal adipose | 45 (15 patients) |
| GSE30929 | survival validation (DFS) | 140 |
| GSE221493 | single-cell atlas, composition-vs-state analysis | 101,890 QC-passed cells / 17 patients |

Portal: <https://www.ncbi.nlm.nih.gov/geo/>

TCGA-SARC data were obtained from the Genomic Data Commons
(<https://portal.gdc.cancer.gov/>). Copy-number cluster and tumour-purity
annotations were adopted from the published annotation accompanying the source
TCGA-SARC study.

---

## Annotation and metadata locks

These are **not** third-party downloads — they are frozen by this project and
are committed under `governance/metadata_lock/` so that the analysis inputs
cannot drift:

| File | Contents |
|---|---|
| `score_genes_v1.tsv` | the six score genes |
| `sample_manifest_v1.tsv` | per-sample cohort and component labels |
| `cell_type_annotation_v1.tsv` | single-cell annotation used by the atlas analysis |
| `pathways_v1.tsv` | pathway definitions |

The `.backup.tsv` counterpart of the sample manifest is retained for
traceability of a manual correction; see `governance/limitations/`.
