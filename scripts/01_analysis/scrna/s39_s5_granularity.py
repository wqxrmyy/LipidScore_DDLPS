## Step 5 粒度诊断：四级粒度（L1 10 / L2 23 / L3 36 / subtype_code 43）
## 核心问题：哪一级粒度能在「不掩盖关键类」与「分解可稳定估计」之间取得平衡
import os, json
import numpy as np
import pandas as pd

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
OUT = []


def A(s=""):
    print(s, flush=True)
    OUT.append(str(s))


df = pd.read_csv(os.path.join(W, "s39_s5_meta.csv"), encoding="utf-8-sig", low_memory=False)
A("=== 载入 s39_s5_meta.csv ===")
A("shape = %s   cols = %s" % (df.shape, list(df.columns)))
A("")

PAIR = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
WC = ["WDLPS", "DDLPS-WD"]     # WD 成分
DC = ["DDLPS-DD"]              # DD 成分

# ---------- [A] 配对样本结构 ----------
A("=" * 78)
A("[A] 配对 7 例的 tissue 构成（QC 后）")
A("=" * 78)
A(pd.crosstab(df.patient, df.tissue).to_string())
A("")

A("-- 每患者 WD/DD 成分细胞数 --")
rows = []
for p in PAIR:
    s = df[df.patient == p]
    nw = int(s.tissue.isin(WC).sum())
    nd = int(s.tissue.isin(DC).sum())
    rows.append({"patient": p, "WD_tissue": "/".join(sorted(set(s.tissue) & set(WC))),
                 "n_WD": nw, "n_DD": nd, "ratio": round(nw / max(nd, 1), 2)})
A(pd.DataFrame(rows).to_string(index=False))
A("")

# ---------- [B] 四级粒度类别清单 ----------
GRANS = ["L1", "L2", "L3", "subtype_code"]
A("=" * 78)
A("[B] 四级粒度的类别数与计数分布（全量 101,890 细胞）")
A("=" * 78)
A("%-14s %8s %10s %10s %10s %10s" % ("粒度", "类别数", "最小类", "p25", "中位", "最大类"))
for g in GRANS:
    vc = df[g].value_counts()
    A("%-14s %8d %10d %10.1f %10.1f %10d" % (g, len(vc), vc.min(), vc.quantile(.25), vc.median(), vc.max()))
A("")

A("-- L2 全类别计数（23 类）--")
A(df.L2.value_counts().to_string())
A("")

A("-- L1 全类别计数（10 类）--")
A(df.L1.value_counts().to_string())
A("")

# ---------- [C] 关键类定位 ----------
A("=" * 78)
A("[C] 关键类（脂肪谱系）逐级定位")
A("=" * 78)
for kw in ["Adipo", "adipo", "Normal"]:
    for g in GRANS:
        hit = [c for c in df[g].unique() if kw in str(c)]
        if hit:
            A("  [%s] 关键词 '%s' 命中: %s" % (g, kw, hit))
A("")
A("-- L3 中含 'Adipo' 的类别计数 --")
A(df[df.L3.astype(str).str.contains("Adipo", case=False, na=False)].L3.value_counts().to_string())
A("")
A("-- subtype_code 的取值形态（前 12 个）--")
A(sorted(df.subtype_code.astype(str).unique())[:12])
A("")
A("-- L2='Adipo diff' 的细胞在 subtype_code 上的分布（该列取值为编码，不含基因名）--")
A(df[df.L2 == "Adipo diff"].subtype_code.value_counts().to_string())
A("")

A("-- 层级一致性检查：下级粒度是否被上级唯一包含 --")
for lo, hi in [("L2", "L1"), ("L3", "L2"), ("subtype_code", "L3")]:
    c = pd.crosstab(df[lo], df[hi])
    n_multi = int(((c > 0).sum(axis=1) > 1).sum())
    pur = float((c.values.max(axis=1) / c.values.sum(axis=1)).mean())
    A("  %-13s -> %-12s  n_%s=%d, 跨多个 %s 的类数=%d, 平均归并纯度=%.6f"
      % (lo, hi, lo, len(c), hi, n_multi, pur))
A("")

# ---------- [D] 稀疏诊断（核心） ----------
A("=" * 78)
A("[D] ★ 稀疏诊断：配对样本 × 类别 的细胞计数矩阵")
A("=" * 78)
A("样本集 = 7 患者 × {WD, DD} = 14 个样本")
A("")
A("%-14s %7s %7s %9s %9s %9s %9s %9s %9s" %
  ("粒度", "类别", "格数", "空格数", "空格率", "N<5率", "N<10率", "N<20率", "可估率"))
summary = {}
for g in GRANS:
    recs = []
    # 每个样本的类别计数
    for p in PAIR:
        s = df[df.patient == p]
        for tag, tl in [("WD", WC), ("DD", DC)]:
            sub = s[s.tissue.isin(tl)]
            vc = sub[g].value_counts()
            vc.name = "%s_%s" % (p, tag)
            recs.append(vc)
    mat = pd.concat(recs, axis=1).fillna(0).astype(int)   # 类别 × 14
    N = mat.values
    cells = N.size
    zero = int((N == 0).sum())
    lt5 = float((N < 5).mean())
    lt10 = float((N < 10).mean())
    lt20 = float((N < 20).mean())
    ok = float((N >= 10).mean())
    A("%-14s %7d %7d %9d %8.1f%% %8.1f%% %8.1f%% %8.1f%% %8.1f%%" %
      (g, len(mat), cells, zero, 100 * zero / cells, 100 * lt5, 100 * lt10, 100 * lt20, 100 * ok))
    summary[g] = {"n_class": int(len(mat)), "n_cell_slots": int(cells), "n_zero": zero,
                  "rate_zero": round(zero / cells, 4), "rate_lt5": round(lt5, 4),
                  "rate_lt10": round(lt10, 4), "rate_lt20": round(lt20, 4),
                  "rate_est_ok": round(ok, 4)}
A("")
A("注：'可估率' = N>=10 的格占比（m_{s,c} 可稳定估计的比例）")
A("")

# 保存计数矩阵供 Step 7 使用
os.makedirs(os.path.join(W, "s5"), exist_ok=True)
for g in GRANS:
    recs = []
    for p in PAIR:
        s = df[df.patient == p]
        for tag, tl in [("WD", WC), ("DD", DC)]:
            sub = s[s.tissue.isin(tl)]
            vc = sub[g].value_counts()
            vc.name = "%s_%s" % (p, tag)
            recs.append(vc)
    mat = pd.concat(recs, axis=1).fillna(0).astype(int)
    mat.index.name = "celltype"
    mat.to_csv(os.path.join(W, "s5", "count_%s.csv" % g), encoding="utf-8-sig")
A("→ 计数矩阵已落盘 _work/s5/count_<粒度>.csv（共 %d 个）" % len(GRANS))
A("")

# ---------- [E] 最细粒度的"尾部"问题 ----------
A("=" * 78)
A("[E] 最细粒度 subtype_code(43) 在两成分间的平衡性")
A("=" * 78)
cnt = []
for p in PAIR:
    s = df[df.patient == p]
    for tag, tl in [("WD", WC), ("DD", DC)]:
        v = s[s.tissue.isin(tl)].subtype_code.value_counts()
        cnt.append(v)
M = pd.concat(cnt, axis=1).fillna(0).astype(int)
present = (M >= 10).sum(axis=1)
A("subtype_code 总数 = %d" % len(M))
A("在 14 个样本中「N>=10 的样本数」分布：")
A(present.value_counts().sort_index(ascending=False).to_string())
A("")
A("仅出现在 <7 个样本（即 <半数）的 subtype_code 数 = %d (%.1f%%)" %
  (int((present < 7).sum()), 100 * (present < 7).mean()))
A("")

A("DONE_S5_GRANULARITY")

rep = os.path.join(W, "s39_s5_report.txt")
with open(rep, "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
with open(os.path.join(W, "s39_s5_summary.json"), "w", encoding="utf-8") as f:
    json.dump({"pair_patients": PAIR, "granularity": summary}, f, ensure_ascii=False, indent=1)
print("\n[written]", rep)
