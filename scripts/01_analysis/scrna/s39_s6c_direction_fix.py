# =============================================================================
#  DDLPS 重构 · 单细胞 #39  Step 6c：§7 描述性方向判读（标签笔误修正重算）
#  ---------------------------------------------------------------------------
#  背景：s39_s6_score.R 原 §7 把配对 WD 的组织标签写成 "WDLPS"（笔误），
#        真实标签为 "DDLPS-WD" → WD 侧全取不到 → §7 输出全 NA。
#        该 bug 仅影响 §7，不影响 §2–§6 的评分与 pseudobulk 产物（均全样本聚合）。
#  修正：口径对齐 Step 3/Step 4 —— 配对 WD = DDLPS-WD（肿瘤内高分化成分），
#        DD = DDLPS-DD。WDLPS 为独立纯高分化病例，不入配对。
#  输入：_work/s6/{pseudobulk_sample,pseudobulk_sample_x_L2,pseudobulk_sample_x_bucket}.csv
#  输出：s39_s6c_direction.txt / s6/pair_wide_{sample,L2,bucket}.csv
#  ★ 本文件全部为描述性判读，判决权归 Step 7 的 Kitagawa 分解。
# =============================================================================
import os
import numpy as np
import pandas as pd
from scipy import stats

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
OUT = os.path.join(W, "s6")
LOG = os.path.join(W, "s39_s6c_direction.txt")

PAIR = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
WD_T, DD_T = "DDLPS-WD", "DDLPS-DD"

lines = []


def A(s=""):
    print(s)
    lines.append(s)


def split_sample(s):
    """'P1_DDLPS-WD' -> ('P1', 'DDLPS-WD')"""
    p, _, t = s.partition("_")
    return p, t


A("=" * 79)
A("Step 6c  §7 描述性方向判读 —— 标签笔误修正重算")
A("=" * 79)
A("")
A("原 §7 失效原因：WD_T 误写为 'WDLPS'（独立纯高分化病例标签）")
A("真实配对标签：WD = 'DDLPS-WD'（DDLPS 肿瘤内高分化成分）｜DD = 'DDLPS-DD'")
A("")

# ---------------------------------------------------------------- 载入
pb_all = pd.read_csv(os.path.join(OUT, "pseudobulk_sample.csv"))
pb_l2 = pd.read_csv(os.path.join(OUT, "pseudobulk_sample_x_L2.csv"))
pb_bk = pd.read_csv(os.path.join(OUT, "pseudobulk_sample_x_bucket.csv"))

for df in (pb_all, pb_l2, pb_bk):
    if "patient" not in df.columns:
        df["patient"] = [split_sample(s)[0] for s in df["sample"]]
        df["tissue"] = [split_sample(s)[1] for s in df["sample"]]

A("-- [0] 样本清单核对 --")
A(f"pseudobulk 样本数 = {pb_all['sample'].nunique()}（应 31）")
A(f"tissue 取值：{sorted(pb_all['tissue'].unique())}")
A("")
chk = pb_all[pb_all.patient.isin(PAIR) & pb_all.tissue.isin([WD_T, DD_T])]
cnt = chk.groupby("patient")["tissue"].apply(lambda s: tuple(sorted(s)))
A("-- [1] 配对完整性（liveness check）--")
ok = 0
for p in PAIR:
    got = cnt.get(p, ())
    full = (WD_T in got) and (DD_T in got)
    ok += full
    A(f"  {p:>4}  WD={'✔' if WD_T in got else '✘'}  DD={'✔' if DD_T in got else '✘'}"
      f"   {'齐全' if full else '★缺失'}")
A(f"完整配对：{ok} / {len(PAIR)}")
A("")

# ---------------------------------------------------------------- (a) 样本级
A("=" * 79)
A("-- [2] (a) 样本级（全体细胞）评分 WD → DD —— 配对 7 例 --")
A("=" * 79)


def paired_table(df, value_col, tag):
    rows = []
    for p in PAIR:
        w = df[(df.patient == p) & (df.tissue == WD_T)][value_col]
        d = df[(df.patient == p) & (df.tissue == DD_T)][value_col]
        rows.append({
            "patient": p,
            "WD": float(w.mean()) if len(w) else np.nan,
            "DD": float(d.mean()) if len(d) else np.nan,
        })
    t = pd.DataFrame(rows)
    t["DD_minus_WD"] = t["DD"] - t["WD"]
    if tag:
        A(f"--- {tag} ---")
        A(f"{'pat':>5} {'WD':>12} {'DD':>12} {'DD-WD':>12}")
        for _, r in t.iterrows():
            A(f"{r['patient']:>5} {r['WD']:>12.4f} {r['DD']:>12.4f} {r['DD_minus_WD']:>12.4f}")
        dd = t["DD_minus_WD"].dropna()
        nneg = int((dd < 0).sum())
        A(f"DD < WD 的配对: {nneg} / {len(dd)}")
        if len(dd) >= 5:
            try:
                _, wp = stats.wilcoxon(t["DD"].to_numpy(), t["WD"].to_numpy())
            except Exception:
                wp = np.nan
            A(f"配对 Wilcoxon（描述性，非判决）: p = {wp:.4f}")
        A("")
    return t


t_ams = paired_table(pb_all, "mean_score_ams", "AddModuleScore")
t_auc = paired_table(pb_all, "mean_score_aucell", "AUCell")

# ---------------------------------------------------------------- (b) 桶级
A("=" * 79)
A("-- [3] (b) 谱系桶级评分（配对）--")
A("=" * 79)
A(f"{'bucket':<10} {'WD':>12} {'DD':>12} {'DD-WD':>12}   {'nWD':>4} {'nDD':>4}")
buckets = ["Adipose", "Tumor", "Immune", "Stroma"]
rows_b = []
for b in buckets:
    s = pb_bk[pb_bk.bucket == b]
    w = s[s.patient.isin(PAIR) & (s.tissue == WD_T)]["mean_score_ams"]
    d = s[s.patient.isin(PAIR) & (s.tissue == DD_T)]["mean_score_ams"]
    A(f"{b:<10} {w.mean():>12.4f} {d.mean():>12.4f} {d.mean() - w.mean():>12.4f}"
      f"   {len(w):>4} {len(d):>4}")


def bucket_paired(df, value_col, b):
    s = df[df.bucket == b]
    rows = []
    for p in PAIR:
        w = s[(s.patient == p) & (s.tissue == WD_T)][value_col]
        d = s[(s.patient == p) & (s.tissue == DD_T)][value_col]
        rows.append({"bucket": b, "patient": p,
                     "WD": float(w.mean()) if len(w) else np.nan,
                     "DD": float(d.mean()) if len(d) else np.nan})
    return pd.DataFrame(rows)


bw = pd.concat([bucket_paired(pb_bk, "mean_score_ams", b) for b in buckets])
bw.to_csv(os.path.join(OUT, "pair_wide_bucket.csv"), index=False, encoding="utf-8-sig")
A("")

# ---------------------------------------------------------------- (c) Adipo diff
A("=" * 79)
A("-- [4] (c) 关键类 Adipo diff（样本级，配对 7 例）--")
A("=" * 79)


def l2_paired(df, value_col, cls):
    s = df[df.L2 == cls]
    rows = []
    for p in PAIR:
        w = s[(s.patient == p) & (s.tissue == WD_T)][value_col]
        d = s[(s.patient == p) & (s.tissue == DD_T)][value_col]
        rows.append({"L2": cls, "patient": p,
                     "WD": float(w.mean()) if len(w) else np.nan,
                     "DD": float(d.mean()) if len(d) else np.nan,
                     "nW": int(s[(s.patient == p) & (s.tissue == WD_T)]["n_cells"].sum()),
                     "nD": int(s[(s.patient == p) & (s.tissue == DD_T)]["n_cells"].sum())})
    return pd.DataFrame(rows)


for cls in ["Adipo diff", "Normal adipocytes", "Invasion", "ECM remodeling"]:
    t = l2_paired(pb_l2, "mean_score_ams", cls)
    t["DD_minus_WD"] = t["DD"] - t["WD"]
    A(f"--- {cls} (AddModuleScore) ---")
    A(f"{'pat':>5} {'WD':>10} {'DD':>10} {'DD-WD':>10} {'nW':>7} {'nD':>7}")
    for _, r in t.iterrows():
        wv = "NA" if pd.isna(r["WD"]) else f"{r['WD']:.4f}"
        dv = "NA" if pd.isna(r["DD"]) else f"{r['DD']:.4f}"
        dlt = "-" if pd.isna(r["DD_minus_WD"]) else f"{r['DD_minus_WD']:.4f}"
        A(f"{r['patient']:>5} {wv:>10} {dv:>10} {dlt:>10} {r['nW']:>7} {r['nD']:>7}")
    dd = t["DD_minus_WD"].dropna()
    A(f"DD < WD: {int((dd < 0).sum())} / {len(dd)}")
    A("")

# 汇总宽表落盘
l2_all = pd.concat([l2_paired(pb_l2, "mean_score_ams", c)
                    for c in sorted(pb_l2.L2.unique())])
l2_all.to_csv(os.path.join(OUT, "pair_wide_L2.csv"), index=False, encoding="utf-8-sig")

t_ams.to_csv(os.path.join(OUT, "pair_wide_sample.csv"), index=False, encoding="utf-8-sig")

A("=" * 79)
A("-- [5] 落盘 --")
A("=" * 79)
A("pair_wide_sample.csv / pair_wide_L2.csv / pair_wide_bucket.csv  →  _work/s6/")
A("★ 以上均为描述性判读；组成 vs 状态的判决归 Step 7 Kitagawa 分解。")
A("")
A("=" * 79)
A("DONE_STEP6C")
A("=" * 79)

with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n[wrote] {LOG}")
