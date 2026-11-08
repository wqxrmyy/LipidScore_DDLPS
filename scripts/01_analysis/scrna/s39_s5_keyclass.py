## Step 5 补充：关键类逐样本细胞数 + 谱系桶映射（Fig5 E/F 行必需）
import os, json
import numpy as np
import pandas as pd

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
OUT = []
def A(s=""):
    print(s, flush=True); OUT.append(str(s))

df = pd.read_csv(os.path.join(W, "s39_s5_meta.csv"), encoding="utf-8-sig", low_memory=False)
PAIR = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
WC, DC = ["WDLPS", "DDLPS-WD"], ["DDLPS-DD"]

samples = []
for p in PAIR:
    for tag, tl in [("WD", WC), ("DD", DC)]:
        s = df[(df.patient == p) & (df.tissue.isin(tl))]
        samples.append((p, tag, s))
labels = ["%s_%s" % (p, t) for p, t, _ in samples]
FMT9 = "%-22s" + "".join(["%9s"] * len(labels))
FMT8 = "%-18s" + "".join(["%9s"] * len(labels))
FMTI = "%-22s" + "".join(["%9d"] * len(labels))

# ---------- [F] 关键脂肪类逐样本细胞数 ----------
A("=" * 78)
A("[F] ★ 评分承载类（脂肪谱系）在 14 个配对样本中的细胞数")
A("=" * 78)
A(FMT9 % ("细胞类型", *labels))
for ct in ["Adipo diff", "Normal adipocytes"]:
    row = [int((s.L2 == ct).sum()) for _, _, s in samples]
    A(FMTI % (ct, *row))
A("")
A("-- 各样本总细胞数 --")
A(FMTI % ("n_cells", *[len(s) for _, _, s in samples]))
A("")

# ---------- [G] L3 的 Adipo diff 三子类逐样本 ----------
A("=" * 78)
A("[G] L3 把 Adipo diff 拆成 3 个子类后的逐样本细胞数（为何 L3 不可作主粒度）")
A("=" * 78)
A(FMT8 % ("L3 子类", *labels))
for ct in ["Adipo diff_1", "Adipo diff_2", "Adipo diff_3"]:
    row = [int((s.L3 == ct).sum()) for _, _, s in samples]
    A(FMT8 % (ct, *row))
A("")
A("→ 若用 L3，状态项 Δm 需对 3 个子类分别估计；Adipo diff_3（全库仅 154 细胞）在多数样本中 N<10，")
A("  其 m_{s,c} 方差极大 → 分解结果被小样本类主导。")
A("")

# ---------- [H] 谱系桶映射（L2 -> 4 桶） ----------
A("=" * 78)
A("[H] ★ L2(23 类) → 谱系桶映射（拟冻结，供 Fig5 E/F 行与分解使用）")
A("=" * 78)
BUCKETS = {
    "Adipose":  ["Adipo diff", "Normal adipocytes"],
    "Tumor":    ["Invasion", "ECM remodeling", "Stemness", "Hypoxia", "Angiogenesis", "Cycling cells"],
    "Immune":   ["Monocytes", "T cells", "NK cells", "B cells", "Plasmocytes", "NKT cells",
                 "DC", "Macrophages", "MDSC", "TAN", "Mast cells"],
    "Stroma":   ["Endothelial cells", "Pericytes", "Mesothelial cells", "RBC"],
}
allL2 = set(df.L2.unique())
mapped = [c for v in BUCKETS.values() for c in v]
A("L2 全部类数 = %d；已映射 = %d；未映射 = %s" % (len(allL2), len(mapped), sorted(allL2 - set(mapped))))
A("重复映射 = %s" % sorted([c for c in mapped if mapped.count(c) > 1]))
A("")
df["_bucket"] = df.L2.map({c: b for b, v in BUCKETS.items() for c in v})
# 重新切片（samples 是在加 _bucket 之前切的，不含该列）
samples2 = []
for p in PAIR:
    for tag, tl in [("WD", WC), ("DD", DC)]:
        samples2.append((p, tag, df[(df.patient == p) & (df.tissue.isin(tl))]))
A("%-10s %10s %8s   %s" % ("桶", "细胞数", "占比", "成员类"))
for b, v in BUCKETS.items():
    n = int((df._bucket == b).sum())
    A("%-10s %10d %7.2f%%   %s" % (b, n, 100 * n / len(df), ", ".join(v)))
A("")
A("-- 各桶在两成分的占比（7 对配对，均值）--")
A("%-10s %10s %10s %10s" % ("桶", "WD均", "DD均", "差(DD-WD)"))
for b in BUCKETS:
    ws, ds = [], []
    for p in PAIR:
        sw = df[(df.patient == p) & (df.tissue.isin(WC))]
        sd = df[(df.patient == p) & (df.tissue.isin(DC))]
        ws.append((sw._bucket == b).mean()); ds.append((sd._bucket == b).mean())
    A("%-10s %9.2f%% %9.2f%% %+9.2f%%" % (b, 100*np.mean(ws), 100*np.mean(ds), 100*(np.mean(ds)-np.mean(ws))))
A("")

# 桶级稀疏诊断
A("-- 桶级稀疏诊断（14 样本 × 4 桶）--")
mat = pd.DataFrame({l: [int((s._bucket == b).sum()) for b in BUCKETS]
                    for l, (_, _, s) in zip(labels, samples2)},
                   index=list(BUCKETS))
A(mat.to_string())
A("空格数 = %d / %d" % (int((mat.values == 0).sum()), mat.size))
A("")

# ---------- [I] 分解类别集预判 ----------
A("=" * 78)
A("[I] 分解可行性预判（L2 粒度，逐类别的可估样本数）")
A("=" * 78)
rows = []
for ct in df.L2.value_counts().index:
    n_ok = sum(1 for _, _, s in samples if (s.L2 == ct).sum() >= 10)
    n_tot = sum(int((s.L2 == ct).sum()) for _, _, s in samples)
    rows.append({"L2": ct, "n_ok_samples(N>=10)": n_ok, "n_in_pairs": n_tot})
r = pd.DataFrame(rows).sort_values("n_ok_samples(N>=10)", ascending=False)
A(r.to_string(index=False))
A("")
A("→ 「两成分均可估」的严格门槛 = 对每个患者在 WD 与 DD 两侧都 N>=10。")
strict = []
for ct in df.L2.value_counts().index:
    ok = all((s.L2 == ct).sum() >= 10 for _, _, s in samples)
    strict.append((ct, ok))
n_ok = sum(1 for _, o in strict if o)
A("严格全可估（14/14 样本均 N>=10）的 L2 类数 = %d / 23" % n_ok)
A("  可估类: %s" % ", ".join([c for c, o in strict if o]))
A("  不可估: %s" % ", ".join([c for c, o in strict if not o]))
A("")

A("DONE_S5_KEYCLASS")
rep = os.path.join(W, "s39_s5_keyclass.txt")
with open(rep, "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
with open(os.path.join(W, "s5", "bucket_map.json"), "w", encoding="utf-8") as f:
    json.dump(BUCKETS, f, ensure_ascii=False, indent=1)
print("[written]", rep)
