# -*- coding: utf-8 -*-
# ============================================================
# Step 9a · Fig5 十二图件的数据层准备
#   数据源（全部为已冻结产物，本脚本只做「搬运 + 聚合 + 交叉核对」，不改口径）：
#     ① 官方嵌入坐标 + 注释：gse221493/GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz
#        列序 = barcode | patient | tissue | location | subtype_code | L1 | L2 | L3 | x | y
#     ② 逐细胞评分（QC 通过 101,890 细胞）：s6/cell_scores.csv.gz
#     ③ 样本×谱系桶 pseudobulk（Step 6 冻结）：s6/pseudobulk_sample_x_bucket.csv
#     ④ Kitagawa 分解（Step 7 冻结）：s7/kitagawa_patient_L2.csv、s24_evidence_pack/T1_main_pooled.csv
#   输出：p5_figures/fig5_data/*.csv（供 fig5a…fig5f4.R 读取；图内数值一律不硬编码）
# ============================================================
import io, os, gzip, json
import numpy as np
import pandas as pd

W   = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
OUT = os.path.join(W, "p5_figures", "fig5_data")
os.makedirs(OUT, exist_ok=True)
REP = []            # 报告行
def A(*a):
    line = " ".join(str(x) for x in a)
    print(line); REP.append(line)

A("=" * 78)
A("Step 9a · Fig5 数据层准备")
A("=" * 78)

PAIR = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
BUCKETS = ["Adipose", "Tumor", "Immune", "Stroma"]

# ---------------------------------------------------------------- ① 官方坐标
bc_path = os.path.join(W, "gse221493",
                       "GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz")
cols = ["barcode", "patient", "tissue", "location", "subtype_code",
        "L1", "L2", "L3", "umap1", "umap2"]
bar = pd.read_csv(bc_path, sep="\t", header=None, names=cols,
                  compression="gzip", quoting=3)
A(f"[1] 官方 barcodes: {bar.shape[0]} 行 × {bar.shape[1]} 列")
assert bar.shape[1] == 10, "barcodes 列数不是 10 —— 结构已变，须人工复核"
A(f"    坐标值域 umap1 [{bar.umap1.min():.3f}, {bar.umap1.max():.3f}] | "
  f"umap2 [{bar.umap2.min():.3f}, {bar.umap2.max():.3f}]")

# ---------------------------------------------------------------- ② 逐细胞评分
cs = pd.read_csv(os.path.join(W, "s6", "cell_scores.csv.gz"))
A(f"[2] cell_scores: {cs.shape[0]} 行（QC 通过）")
assert cs.shape[0] == 101890

# 键一致性：barcode 唯一且能全量对齐
assert bar.barcode.is_unique and cs.barcode.is_unique
inter = set(cs.barcode) & set(bar.barcode)
A(f"    barcode 交集 = {len(inter)} / cell_scores {len(cs)} / barcodes {len(bar)}")
assert len(inter) == len(cs), "存在 QC 通过但官方坐标缺失的细胞"

# 标签一致性：L1/L2 两处必须逐细胞相同（防两套注释漂移）
m = cs.merge(bar[["barcode", "L1", "L2", "patient", "tissue"]],
             on="barcode", suffixes=("_cs", "_bar"))
assert len(m) == len(cs)
bad1 = int((m.L1_cs != m.L1_bar).sum()); bad2 = int((m.L2_cs != m.L2_bar).sum())
badp = int((m.patient_cs != m.patient_bar).sum()); badt = int((m.tissue_cs != m.tissue_bar).sum())
A(f"    ★ 注释一致性: L1 不符 {bad1} | L2 不符 {bad2} | patient {badp} | tissue {badt}")
assert bad1 == 0 and bad2 == 0 and badp == 0 and badt == 0, "两套注释不一致"

# ---------------------------------------------------------------- ①′ UMAP 图件数据
um = bar[["barcode", "umap1", "umap2", "L1", "L2", "patient", "tissue"]].merge(
    cs[["barcode", "bucket", "score_ams", "score_aucell"]], on="barcode", how="inner")
um = um.rename(columns={"L1": "celltype"})
# L1 类按平均评分降序 → 同一顺序贯穿 A 图例 / C 点图 / D 箱线
l1_order = (um.groupby("celltype").score_ams.mean()
              .sort_values(ascending=False).index.tolist())
l1_n = um.celltype.value_counts()
A(f"[3] UMAP 数据: {um.shape[0]} 细胞 | L1 共 {len(l1_order)} 类")
for c in l1_order:
    A(f"      {c:22} n={l1_n[c]:>6}  mean score A = {um.loc[um.celltype==c,'score_ams'].mean():+.4f}")
um.to_csv(os.path.join(OUT, "umap_cells.csv.gz"), index=False,
          compression="gzip", float_format="%.6f")
A(f"    → umap_cells.csv.gz")

# ---------------------------------------------------------------- ③ 组成层（占比）
cnt = cs.groupby(["sample", "patient", "tissue", "bucket"]).size().rename("n_cells").reset_index()
tot = cnt.groupby("sample").n_cells.sum().rename("n_total").reset_index()
cnt = cnt.merge(tot, on="sample")
cnt["frac"] = cnt.n_cells / cnt.n_total
cnt.to_csv(os.path.join(OUT, "comp_sample_bucket.csv"), index=False,
           encoding="utf-8-sig", float_format="%.6f")
A(f"[4] 组成层: {len(cnt)} 行（{cnt.tissue.nunique()} 组织 × {cnt['bucket'].nunique()} 桶）")

# 配对完整度
pc = cnt[cnt.patient.isin(PAIR) & cnt.tissue.isin(["DDLPS-WD", "DDLPS-DD"])]
n_ok = sum(1 for p in PAIR
           if set(pc.loc[pc.patient == p, "tissue"]) == {"DDLPS-WD", "DDLPS-DD"})
A(f"    配对完整度 = {n_ok}/7")
assert n_ok == 7

# ---------------------------------------------------------------- ④ 状态层（谱系内均值评分）
# 读 Step 6 冻结的 pseudobulk，并与逐细胞重算交叉核对（必须逐值相同）
pb = pd.read_csv(os.path.join(W, "s6", "pseudobulk_sample_x_bucket.csv"))
chk = (cs.groupby(["sample", "bucket"])
         .agg(mean_score_recalc=("score_ams", "mean"), n_recalc=("score_ams", "size"))
         .reset_index())
pb2 = pb.merge(chk, on=["sample", "bucket"], how="outer")
assert pb2.n_cells.notna().all() and pb2.n_recalc.notna().all(), "pseudobulk 与逐细胞计数对不上"
dmax = float(np.abs(pb2.mean_score_ams - pb2.mean_score_recalc).max())
dn = int((pb2.n_cells.astype(int) != pb2.n_recalc.astype(int)).sum())
A(f"[5] 状态层: pseudobulk {len(pb)} 行 | ★ 与逐细胞重算：最大偏差 {dmax:.3e} | 计数不符 {dn} 行")
assert dmax < 1e-12 and dn == 0, "pseudobulk 与逐细胞重算不一致"
st = cnt[["sample", "patient", "tissue", "bucket", "n_cells", "n_total", "frac"]].merge(
    pb[["sample", "bucket", "mean_score_ams", "mean_score_aucell"]], on=["sample", "bucket"])
st.to_csv(os.path.join(OUT, "state_sample_bucket.csv"), index=False,
          encoding="utf-8-sig", float_format="%.6f")
A(f"    → state_sample_bucket.csv")

# 与 Step 6 报告「桶级 Δ(DD−WD)」对照（7 对患者上 DD−WD 的均值）
rep6 = {"Adipose": -0.0395, "Tumor": -0.0006, "Immune": -0.0086, "Stroma": -0.0622}
w = st[st.patient.isin(PAIR) & st.tissue.isin(["DDLPS-WD", "DDLPS-DD"])]
piv = w.pivot_table(index=["bucket", "patient"], columns="tissue",
                    values="mean_score_ams").reset_index()
piv["d"] = piv["DDLPS-DD"] - piv["DDLPS-WD"]
A("    ★ 桶级 Δ(DD−WD) 与 Step 6 报告对照：")
for b in BUCKETS:
    got = piv.loc[piv.bucket == b, "d"].mean()
    A(f"      {b:8} 本脚本 {got:+.4f} | Step 6 报告 {rep6[b]:+.4f} | 差 {abs(got-rep6[b]):.5f}")
    assert abs(got - rep6[b]) < 5e-5, f"{b} 桶级 Δ 与 Step 6 报告不一致"

# ---------------------------------------------------------------- ⑤ E4 脂肪占比梯度
# 参照序列：nl / lipoma（纯脂肪上限）→ WDLPS → DDLPS-sWD → DDLPS-WD → DDLPS-DD
GRP_ORDER = ["nl", "lipoma", "WDLPS", "DDLPS-sWD", "DDLPS-WD", "DDLPS-DD"]
ad = cnt[cnt["bucket"] == "Adipose"].copy()
excluded = sorted(set(cnt.tissue) - set(GRP_ORDER))
A(f"[6] E4 梯度：保留组织 {GRP_ORDER}")
A(f"    ⚠️ 未纳入本面板的组织（须显式披露，非静默剔除）：{excluded}")
ad["grp"] = pd.Categorical(ad.tissue, categories=GRP_ORDER, ordered=True)
ad = ad.dropna(subset=["grp"]).sort_values(["grp", "frac"])
ad[["sample", "patient", "tissue", "grp", "n_cells", "n_total", "frac"]].to_csv(
    os.path.join(OUT, "adipose_gradient.csv"), index=False,
    encoding="utf-8-sig", float_format="%.6f")
A("    各组脂肪系占比（中位）：")
for g in GRP_ORDER:
    s = ad.loc[ad.grp == g, "frac"]
    A(f"      {g:12} n={len(s):>2}  中位 {s.median():.4f}  范围 [{s.min():.4f}, {s.max():.4f}]")

# ★ 诊断：脂肪系桶由哪两个 L2 类构成（决定 E4 的判读口径，防止把「Adipo diff」误读为正常脂肪细胞）
A("    ★ 脂肪系桶的 L2 构成（占该样本总细胞的比，按组中位）：")
sub = cs[cs["bucket"] == "Adipose"]
l2cnt = (sub.groupby(["sample", "L2"]).size().rename("n").reset_index()
            .merge(tot, on="sample"))
l2cnt["frac"] = l2cnt.n / l2cnt.n_total
l2cnt["tissue"] = l2cnt["sample"].map(dict(zip(cs["sample"], cs["tissue"])))
l2cnt = l2cnt[l2cnt.tissue.isin(GRP_ORDER)]
l2cnt["grp"] = pd.Categorical(l2cnt.tissue, categories=GRP_ORDER, ordered=True)
pv = l2cnt.pivot_table(index="grp", columns="L2", values="frac", aggfunc="median")
pv = pv.reindex(GRP_ORDER)
for g in GRP_ORDER:
    cells = " | ".join(f"{c}={pv.loc[g, c]:.4f}" for c in pv.columns if pd.notna(pv.loc[g, c]))
    A(f"      {g:12} {cells}")
l2cnt[["sample", "tissue", "grp", "L2", "n", "n_total", "frac"]].to_csv(
    os.path.join(OUT, "adipose_L2_breakdown.csv"), index=False,
    encoding="utf-8-sig", float_format="%.6f")

# ---------------------------------------------------------------- ⑥ F4 Kitagawa
kp = pd.read_csv(os.path.join(W, "s7", "kitagawa_patient_L2.csv"))
t1 = pd.read_csv(os.path.join(W, "s24_evidence_pack", "T1_main_pooled.csv"))
A(f"[7] Kitagawa 逐患者表列: {list(kp.columns)}")
A(f"    合计表列: {list(t1.columns)}")
A(t1.to_string(index=False))
# 注意：必须用 kp["T"] —— kp.T 是 DataFrame 转置属性，会被静默当成转置（abs 报 str 错）
cC, cT = kp["C"].abs(), kp["T"].abs()
kp["ratio"] = cC / (cC + cT)
kp["dS_check"] = kp["C"] + kp["T"]
assert float((kp["dS_check"] - kp["dS"]).abs().max()) < 1e-12, "dS 与 C+T 不符"
kp = kp.sort_values("dS").reset_index(drop=True)
kp.to_csv(os.path.join(OUT, "kitagawa_patient_F4.csv"), index=False,
          encoding="utf-8-sig", float_format="%.6f")
A("    逐患者（按 ΔS 升序）：")
for _, r in kp.iterrows():
    A(f"      {r.patient:4} C={r['C']:+.5f} T={r['T']:+.5f} dS={r['dS']:+.5f} ratio={r['ratio']:.4f}")

# 主判据与裁定值对照（防漂移）
rat = abs(t1["C_total"].iloc[0]) / (abs(t1["C_total"].iloc[0]) + abs(t1["T_total"].iloc[0]))
A(f"[8] ★ 复核主判据 ratio = {rat:.4f}（D23 记录 0.8663）| 逐患者中位 {kp.ratio.median():.4f}（D23 记录 0.5370）")
assert abs(rat - 0.8663) < 5e-4 and abs(kp["ratio"].median() - 0.5370) < 5e-4

# L1 顺序（供 C/D 面板复用）
with io.open(os.path.join(OUT, "l1_order.json"), "w", encoding="utf-8") as f:
    json.dump(l1_order, f, ensure_ascii=False, indent=1)

with io.open(os.path.join(W, "s39_s9a_prep_report.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(REP) + "\n")
A("")
A("DONE · 产物 → p5_figures/fig5_data/")
