# =============================================================================
#  DDLPS 重构 · 单细胞 #39  Step 7：Kitagawa 分解（★ 判决层）
#  ---------------------------------------------------------------------------
#  依据：plan_scrna_composition_vs_state.md §2 Q2 + §1.10⑦（规则已冻结）
#    对称分解  ΔS_p = Σ_c [ Δf_{p,c}·m̄_{p,c} + f̄_{p,c}·Δm_{p,c} ]
#                      └──── 组成项 C_p ────┘   └──── 状态项 T_p ────┘
#    Δf = f_DD,c − f_WD,c ;  Δm = m_DD,c − m_WD,c ;  m̄/f̄ = 两成分均值
#
#  ★ 预设判决线（分析前登记，禁止事后调整）
#    |C| / (|C|+|T|) ≥ 0.8  → H_comp 主导（支柱1 = 组织学组成的读数）
#                     ≤ 0.5  → H_state 主导（升级：克隆内脂质身份可塑性）
#                     0.5–0.8 → 混合（定量报告比例）
#
#  冻结的粒度与缺失规则（§1.10⑦）：
#    主粒度 = L2(23)；敏感性粒度 = L1 / L3
#    m 缺失填充：主分析 = pooled 均值（14 配对样本、按细胞数加权）并报填充规模
#                敏感性A = counterpart（Δm=0，保守偏 H_comp）
#                敏感性B = complete（整类剔除 + f 重归一化）
#  输入：_work/s6/cell_scores.csv.gz, _work/s6/pseudobulk_sample_x_{L2,L3}.csv
#  输出：s39_s7_report.txt / s7/*.csv
# =============================================================================
import os
import json
import numpy as np
import pandas as pd
from scipy import stats

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
OUT = os.path.join(W, "s7")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(W, "s39_s7_report.txt")

PAIR = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
WDT, DDT = "DDLPS-WD", "DDLPS-DD"
lines = []


def A(s=""):
    print(s)
    lines.append(str(s))


A("=" * 79)
A("Step 7  Kitagawa 分解（判决层）")
A("=" * 79)
A("")

# ------------------------------------------------------------------ 载入
cells = pd.read_csv(os.path.join(W, "s6", "cell_scores.csv.gz"))
A("-- [0] 载入 --")
A(f"细胞级表: {len(cells):,} 细胞 × {cells.shape[1]} 列")
A(f"tissue 取值: {sorted(cells['tissue'].unique())}")
A("")


def agg_from_cells(cls_col):
    g = (cells[["sample", "patient", "tissue", cls_col, "score_ams", "score_aucell"]]
         .groupby(["sample", "patient", "tissue", cls_col], as_index=False)
         .agg(n=("score_ams", "size"), m=("score_ams", "mean"),
              m_auc=("score_aucell", "mean")))
    return g.rename(columns={cls_col: "cls"})


def read_l3():
    d = pd.read_csv(os.path.join(W, "s6", "pseudobulk_sample_x_L3.csv"))
    d = d.rename(columns={"n_cells": "n", "mean_score_ams": "m",
                          "mean_score_aucell": "m_auc", "L3": "cls"})
    d["patient"] = [s.split("_")[0] for s in d["sample"]]
    d["tissue"] = ["_".join(s.split("_")[1:]) for s in d["sample"]]
    return d[["sample", "patient", "tissue", "cls", "n", "m", "m_auc"]]


TABLES = {
    "L1": agg_from_cells("L1"),
    "L2": agg_from_cells("L2"),
    "L3": read_l3(),
}
for k, v in TABLES.items():
    A(f"  粒度 {k}: {v['cls'].nunique()} 类 / {len(v)} 行")

# ------------------------------------------------------- 配对完整性与样本集
A("")
A("-- [1] 配对完整性断言（★ 分组标签实测，沿用 Step 6c 教训）--")
ok = 0
for p in PAIR:
    s = TABLES["L2"]
    hw = ((s.patient == p) & (s.tissue == WDT)).any()
    hd = ((s.patient == p) & (s.tissue == DDT)).any()
    ok += (hw and hd)
    A(f"  {p:>4}  WD={'ok' if hw else 'MISS'}  DD={'ok' if hd else 'MISS'}")
A(f"完整配对: {ok} / {len(PAIR)}")
assert ok == len(PAIR), "配对不完整"
A("")

# ------------------------------------------------------------ 分解核心函数
def decompose(g, fill="pooled"):
    """g: long 表(sample,patient,tissue,cls,n,m). 返回 (per-patient list, class detail)"""
    res, detail = [], []
    # pooled 类均值（跨 14 个配对样本、按细胞数加权）
    gp = g[g.patient.isin(PAIR) & g.tissue.isin([WDT, DDT]) & (g.n > 0)]
    pooled = ((gp.assign(w=gp.n * gp.m)
               .groupby("cls")[["w", "n"]].sum()
               .eval("w / n")).to_dict())
    filled_cnt = 0
    for p in PAIR:
        w = g[(g.patient == p) & (g.tissue == WDT)].set_index("cls")
        d = g[(g.patient == p) & (g.tissue == DDT)].set_index("cls")
        K = sorted(set(w.index) | set(d.index))
        nw_tot = float(w["n"].sum())
        nd_tot = float(d["n"].sum())
        fw = (w["n"] / nw_tot).reindex(K).fillna(0.0)
        fd = (d["n"] / nd_tot).reindex(K).fillna(0.0)
        mw = w["m"].reindex(K).astype(float)
        md = d["m"].reindex(K).astype(float)

        if fill == "pooled":
            ps = pd.Series(pooled).reindex(K)
            filled_cnt += int(mw.isna().sum() + md.isna().sum())
            mw = mw.fillna(ps)
            md = md.fillna(ps)
        elif fill == "counterpart":
            mw2 = mw.where(mw.notna(), md)
            md2 = md.where(md.notna(), mw)
            mw, md = mw2, md2
        elif fill == "complete":
            keep = [c for c in K if (c in w.index) and (c in d.index)]
            if len(keep) == 0:
                continue
            fw, fd = fw[keep], fd[keep]
            mw, md = mw[keep], md[keep]
            fw = fw / fw.sum()
            fd = fd / fd.sum()
            K = keep
        else:
            raise ValueError(fill)

        df_ = fd - fw
        dm_ = md - mw
        mbar = (mw + md) / 2.0
        fbar = (fw + fd) / 2.0
        cC = df_ * mbar
        cT = fbar * dm_
        C = float(cC.sum())
        T = float(cT.sum())
        dS = float((fd * md - fw * mw).sum())
        res.append({"patient": p, "C": C, "T": T, "dS": dS,
                    "C_plus_T": C + T, "identity_resid": abs(C + T - dS)})
        for c in K:
            detail.append({"patient": p, "cls": c, "f_WD": fw[c], "f_DD": fd[c],
                           "m_WD": mw[c], "m_DD": md[c], "df": df_[c], "dm": dm_[c],
                           "contrib_C": cC[c], "contrib_T": cT[c]})
    return res, pd.DataFrame(detail), filled_cnt


def verdict(ratio):
    if ratio >= 0.8:
        return "H_comp 主导（支柱1 = 组织学组成的读数）"
    if ratio <= 0.5:
        return "H_state 主导（升级：克隆内脂质身份可塑性）"
    return "混合（定量报告比例）"


# --------------------------------------------------------------- 主分析 L2
A("=" * 79)
A("-- [2] 主分析：L2（23 类）粒度 + pooled 填充 --")
A("=" * 79)
res, det, filled = decompose(TABLES["L2"], "pooled")
A(f"填充单元数（m 缺失）: {filled}")
A("")
A(f"{'pat':>5} {'C_p':>11} {'T_p':>11} {'ΔS_p':>11} {'C+T':>11} {'|C+T-ΔS|':>10}")
for r in res:
    A(f"{r['patient']:>5} {r['C']:>11.5f} {r['T']:>11.5f} {r['dS']:>11.5f} "
      f"{r['C_plus_T']:>11.5f} {r['identity_resid']:>10.2e}")
C_tot = sum(r["C"] for r in res)
T_tot = sum(r["T"] for r in res)
dS_tot = sum(r["dS"] for r in res)
A("")
A(f"合计:  C = {C_tot:+.5f}   T = {T_tot:+.5f}   ΔS = {dS_tot:+.5f}")
A(f"恒等式校验: |C+T-ΔS| = {abs(C_tot + T_tot - dS_tot):.2e}  （应为 ~0）")
A("")

ratio = abs(C_tot) / (abs(C_tot) + abs(T_tot))
signed = C_tot / (C_tot + T_tot)
per_pat = [abs(r["C"]) / (abs(r["C"]) + abs(r["T"])) for r in res]
A("-- ★ 判决 --")
A(f"主判据  |C|/(|C|+|T|) = {ratio:.4f}   （方案书 §2 Q2 预设口径）")
A(f"带符号贡献占比  C/(C+T) = {signed:.4f}")
A(f"逐患者比值  中位数 {np.median(per_pat):.4f}  范围 [{min(per_pat):.4f}, {max(per_pat):.4f}]")
A(f"             全部 7 例: " + "  ".join(f"{v:.3f}" for v in per_pat))
A(f"")
A(f"判决线 ≥0.8 = H_comp ｜ ≤0.5 = H_state ｜ 0.5–0.8 = 混合")
A(f"→ ★ 判决结果：{verdict(ratio)}")
A("")

A("-- 配对检验（Wilcoxon signed-rank, n = 7）--")
try:
    s1, p1 = stats.wilcoxon([r["C"] for r in res], [r["T"] for r in res])
    A(f"  C_p vs T_p : W = {s1:.1f},  p = {p1:.4f}")
except Exception as e:
    A(f"  C_p vs T_p : 无法计算（{e}）")
try:
    s2, p2 = stats.wilcoxon([r["dS"] for r in res])
    A(f"  ΔS_p vs 0  : W = {s2:.1f},  p = {p2:.4f}   （总方向是否显著）")
except Exception as e:
    A(f"  ΔS_p vs 0  : 无法计算（{e}）")
n_neg = sum(1 for r in res if r["dS"] < 0)
A(f"  ΔS < 0 的配对: {n_neg} / {len(res)}   （bulk 侧为 6/8 负）")
A("")

# --------------------------------------------------------- 类级贡献明细
A("-- L2 类级贡献（7 例合计，按 |组成项| 降序前 12）--")
agg = (det.groupby("cls", as_index=False)
       .agg(contrib_C=("contrib_C", "sum"), contrib_T=("contrib_T", "sum"),
            f_WD=("f_WD", "mean"), f_DD=("f_DD", "mean"),
            m_WD=("m_WD", "mean"), m_DD=("m_DD", "mean")))
agg["df"] = agg["f_DD"] - agg["f_WD"]
agg["dm"] = agg["m_DD"] - agg["m_WD"]
agg["total"] = agg["contrib_C"] + agg["contrib_T"]
agg = agg.sort_values("contrib_C", key=lambda s: s.abs(), ascending=False)
A(f"{'L2 类':<24} {'f_WD':>7} {'f_DD':>7} {'Δf(pp)':>8} {'Δm':>8} "
  f"{'组成项':>10} {'状态项':>10}")
for _, r in agg.head(12).iterrows():
    A(f"{r['cls']:<24} {r['f_WD']:>7.4f} {r['f_DD']:>7.4f} {r['df']*100:>8.2f} "
      f"{r['dm']:>8.4f} {r['contrib_C']:>10.5f} {r['contrib_T']:>10.5f}")
A("")

# ------------------------------------------------------------- 敏感性分析
A("=" * 79)
A("-- [3] 敏感性分析（填充规则 × 粒度，均已预先登记）--")
A("=" * 79)
A(f"{'粒度':<6} {'填充':<12} {'C':>11} {'T':>11} {'ratio':>8} {'判决':<10}")
sens_rows = []
for gr in ["L2", "L1", "L3"]:
    for fl in ["pooled", "counterpart", "complete"]:
        try:
            r2, _, fc = decompose(TABLES[gr], fl)
            c = sum(x["C"] for x in r2)
            t = sum(x["T"] for x in r2)
            rt = abs(c) / (abs(c) + abs(t))
            vd = ("H_comp" if rt >= 0.8 else "H_state" if rt <= 0.5 else "混合")
            mark = " ←主" if (gr == "L2" and fl == "pooled") else ""
            A(f"{gr:<6} {fl:<12} {c:>11.5f} {t:>11.5f} {rt:>8.4f} {vd:<10}{mark}")
            sens_rows.append({"granularity": gr, "fill": fl, "C": c, "T": t,
                              "ratio": rt, "verdict": vd, "n_filled": fc})
        except Exception as e:
            A(f"{gr:<6} {fl:<12}  失败: {e}")
A("")

# ------------------------------------------------------ nl 参照（LIM-008）
A("=" * 79)
A("-- [4] 纯脂肪参照（LIM-P3-008：Adipo diff 是否只是脂肪谱系读数）--")
A("=" * 79)
ref = cells[cells["tissue"].isin(["nl", "lipoma"])]
A(f"{'样本':<16} {'类':<20} {'占比%':>8} {'评分均值':>10} {'n':>7}")
for s in sorted(ref["sample"].unique()):
    sub = ref[ref["sample"] == s]
    tot = len(sub)
    for cls in ["Adipo diff", "Normal adipocytes"]:
        k = sub[sub["L2"] == cls]
        if len(k) == 0:
            A(f"{s:<16} {cls:<20} {'0.0':>8} {'NA':>10} {0:>7}")
            continue
        A(f"{s:<16} {cls:<20} {len(k)/tot*100:>8.2f} {k['score_ams'].mean():>10.4f} {len(k):>7}")

pair_cells = cells[cells.patient.isin(PAIR) & cells.tissue.isin([WDT, DDT])]
A("")
A(f"{'配对成分':<20} {'类':<20} {'占比%':>8} {'评分均值':>10}")
for tt, lab in [(WDT, "对配 WD"), (DDT, "对配 DD")]:
    sub = pair_cells[pair_cells.tissue == tt]
    tot = len(sub)
    for cls in ["Adipo diff", "Normal adipocytes"]:
        k = sub[sub["L2"] == cls]
        A(f"{lab:<20} {cls:<20} {len(k)/tot*100:>8.2f} {k['score_ams'].mean():>10.4f}")
A("")

# ------------------------------------------------------------------ 落盘
pd.DataFrame(res).to_csv(os.path.join(OUT, "kitagawa_patient_L2.csv"),
                         index=False, encoding="utf-8-sig")
agg.to_csv(os.path.join(OUT, "kitagawa_class_L2.csv"),
           index=False, encoding="utf-8-sig")
det.to_csv(os.path.join(OUT, "kitagawa_detail_L2.csv"),
           index=False, encoding="utf-8-sig")
pd.DataFrame(sens_rows).to_csv(os.path.join(OUT, "kitagawa_sensitivity.csv"),
                               index=False, encoding="utf-8-sig")

summary = {
    "date": "2026-09-22", "granularity_main": "L2", "fill_main": "pooled",
    "n_pairs": len(res), "n_filled": filled,
    "C_total": C_tot, "T_total": T_tot, "dS_total": dS_tot,
    "identity_resid": abs(C_tot + T_tot - dS_tot),
    "ratio_abs": ratio, "ratio_signed": signed,
    "per_patient_ratio": per_pat, "per_patient_ratio_median": float(np.median(per_pat)),
    "verdict": verdict(ratio),
    "dS_negative_pairs": n_neg,
}
with open(os.path.join(W, "s39_s7_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

A("=" * 79)
A("-- [5] 落盘 --")
A("=" * 79)
A("s7/kitagawa_{patient,class,detail}_L2.csv ｜ s7/kitagawa_sensitivity.csv")
A("s39_s7_summary.json")
A("")
A("=" * 79)
A("DONE_STEP7")
A("=" * 79)

with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print(f"\n[wrote] {LOG}")
