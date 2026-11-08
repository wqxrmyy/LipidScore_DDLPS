# 01_analysis — 分析代码 / analysis code

> 本目录是**公开发布集**，不是工作档案。原始归档共 141 个分析脚本，其中 41 个的职责是写稿件/记忆/局限登记/图注（属**溯源过程文件**），已剔除；剔除清单见 `../SELECTION.json` 的 `dropped_analysis`。

| 组 / group | 目录 | 脚本数 | 覆盖 |
|---|---|---:|---|
| 单细胞 / scRNA-seq (GSE221493) | `scrna/` | 51 | 10x 矩阵 → QC/标准化 → UMAP → Harmony → L1/L2 注释 → 评分投影 → Kitagawa 分解 |
| TCGA-SARC | `tcga/` | 31 | 拷贝数簇 / 表达 / 纯度 / 稳健性 → 评分修订与溯源 → Harrell C 与 bootstrap |
| GSE 队列验证 / GSE cohort validation | `validation/` | 5 | GSE30929 / GSE21122 / GSE159659 的 HR、C-index 与方向核对 |
| 数据获取 / Data acquisition | `fetch/` | 6 | GEO/TCGA 下载与解包、队列分层、文献交叉核对 |
| 评分口径与终点 / Score caliber & endpoints | `caliber/` | 7 | 评分公式、Ensembl 校验、四终点（OS/DSS/RFS/NTE）、GSE30929 亚组复核 |

## 路径口径 / how paths work

脚本中的项目根不再写死为开发者机器上的绝对路径，而是运行时读取环境变量：

| 变量 | 含义 | 未设置时 |
|---|---|---|
| `LPS_DATA_ROOT` | 项目数据根（`data/`、`results/` 所在） | `.` |
| `LPS_WORK_ROOT` | 内部工作树（历史 `_work/`） | `.` |
| `LPS_FIG_ROOT` | 出图工作根，对应 `scripts/02_figures/` | `.` |
| `R_BIN_DIR` | R 可执行文件目录（少数脚本调用 `Rscript`） | `.` |

```sh
export LPS_DATA_ROOT=$PWD        # repo root also holds data/ after download
export LPS_FIG_ROOT=$PWD/scripts/02_figures
```

## 关键链 / chains reported in the manuscript

```
单细胞  s39_s1_matrix_dl → s3_read10x → s3b1_qc_norm → s3b2_umap → s4_batch_harmony
        → s5_keyclass → s6_score → s7_kitagawa → s9a_prep
TCGA    d7_cncluster → d12_p4_robustness → d15*_score_revision
        → s22_harrell_A_baseline → s22b_delta_boot → s149_tcga_endpoints_analysis
```

> 与仓库既有模块的关系：`scripts/07_bulk_validation/` 是 GSE221492 的**双实现互证**模块，`scripts/09_qc/` 断言稿件数值；本目录是产生图与表的分析代码。
