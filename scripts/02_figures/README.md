# 02_figures — 出图代码 / figure code

每个子图一个脚本，按图分目录；输出按**最终印刷尺寸 1:1** 生成 PNG + 矢量 PDF。

| 图 | 面板 | 目录 | 脚本数 |
|---|---:|---|---:|
| Fig 1 | 7 | `fig1/` | 6 |
| Fig 2 | 8 | `fig2/` | 7 |
| Fig 3 | 6 | `fig3/` | 6 |
| Fig 4 | 7 | `fig4/` | 7 |
| Fig 5 | 12 | `fig5/` | 12 |
| S1 Fig | 5 | `figS1/` | 5 |
| S2 Fig | 4 | `figS2/` | 4 |
| S3 Fig | 9 | `figS3/` | 9 |

> 脚本与面板不总是一一对应：`fig1b.R` 产 `fig1b1`+`fig1b2`，`fig2g.R` 产 `fig2g1`+`fig2g2`；`fig5e1–e4`、`fig5f1–f4` 各自独立成脚本。逐面板对应以 `panel_manifest.csv` 的 `stem` 为准。

## 共享库 `utils/`

| 文件 | 作用 |
|---|---|
| `colors.R` | 全图配色唯一真相源 |
| `themes.R` | `theme_publication()`、`panel_box()`、`export_fig()` |
| `assemble_lib.R` | 拼版公共库 |
| `prep_data.R` / `prep_gse221492.R` | TCGA / GSE221492 侧数据准备 |
| `fig5_common.R` / `figS6_common.R` | Fig 5 / Fig S6 面板辅助函数 |

## 装配 `assemble/`

| 文件 | 作用 |
|---|---|
| `make_manual_kit.py` | 主装配器：按 `panel_manifest.csv` 把面板贴成整图 |
| `_render_all.py` | 按最终印刷尺寸批量重渲染全部面板 |
| `assemble_fig_main_v3.R` / `assemble_figS_v3.R` | 拼版实现（v3 合规版） |

> 未收录：8 个 `_asm_*_reference.R`（对账用参考装配）、`assemble_fig_main.R`/`assemble_figS.R`（被 v3 取代的早期实现）、`_dim_scan.py`/`_capacity_calc.py`（版面尺寸试算）。

## 运行方式

路径相对 `LPS_FIG_ROOT`（默认 `.`）：

```sh
cd scripts/02_figures
export LPS_FIG_ROOT=$PWD
Rscript fig1/fig1a.R
```
