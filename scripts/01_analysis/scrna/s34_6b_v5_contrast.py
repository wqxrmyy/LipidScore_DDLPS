# -*- coding: utf-8 -*-
"""S34 6B v5 对照：rho(A) vs rho(A-PPARG)，从落盘 CSV 重算（可追溯）"""
import csv, os, math

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
SRC = os.path.join(W, "p5_figures", "fig6_data", "pathway_vs_score.csv")
OUT = os.path.join(W, "s34_6b_v5_contrast.txt")
assert os.path.exists(SRC), SRC

rows = []
with open(SRC, "r", encoding="utf-8-sig", newline="") as f:
    for r in csv.DictReader(f):
        for k in ("rho_tum", "rho_tum_noPP", "rho_tum_FULL", "rho_own",
                  "P_tum", "P_tum_noPP", "P_tum_FULL", "P_own", "q_tum", "q_tum_noPP"):
            if k in r and r[k] not in (None, ""):
                r[k] = float(r[k])
        r["n"] = int(float(r["n"]))
        r["d_noPP"] = r["rho_tum_noPP"] - r["rho_tum"]
        r["d_clean"] = r["rho_tum"] - r["rho_tum_FULL"]
        rows.append(r)

assert len(rows) == 260, f"expect 260 rows, got {len(rows)}"

def line(s=""):
    print(s)

L = []
def A(s=""):
    L.append(s)

A("S34 6B v5 contrast  rho(A) vs rho(A-PPARG)  [recomputed from pathway_vs_score.csv]")
A(f"rows={len(rows)}  pools={sorted(set(r['pool'] for r in rows))}")
A("")
A("== 1. 循环确证：去 PPARG 后 |delta| 最大的 25 项 ==")
A(f"{'pool':<14} {'pathway':<46} {'n':>3} {'rho(A)':>9} {'rho(A-PPARG)':>13} {'delta':>9} {'P_noPP':>10}")
for r in sorted(rows, key=lambda z: -abs(z["d_noPP"]))[:25]:
    A(f"{r['pool']:<14} {r['pathway']:<46} {r['n']:>3} {r['rho_tum']:>+9.4f} {r['rho_tum_noPP']:>+13.4f} {r['d_noPP']:>+9.4f} {r['P_tum_noPP']:>10.3g}")

A("")
A("== 2. mTOR / PPAR / 脂肪酸 专项对照（全部 4 池）==")
sel = [r for r in rows if any(t in r["pathway"] for t in
       ("MTOR", "PPAR", "FATTY_ACID", "UNSATURATED", "ADIPOCYTO", "INSULIN"))]
A(f"{'pool':<14} {'pathway':<46} {'n':>3} {'rho(A)':>9} {'rho(A-PPARG)':>13} {'delta':>9}")
for r in sorted(sel, key=lambda z: (z["pathway"], z["pool"])):
    A(f"{r['pool']:<14} {r['pathway']:<46} {r['n']:>3} {r['rho_tum']:>+9.4f} {r['rho_tum_noPP']:>+13.4f} {r['d_noPP']:>+9.4f}")

A("")
A("== 3. 炎症 / 免疫通路（对照 N3 的『脂质高 ↔ 炎症低』）==")
sel2 = [r for r in rows if any(t in r["pathway"] for t in
        ("INFLAMMATORY", "INTERFERON", "COMPLEMENT", "ALLOGRAFT", "IL6_JAK", "KRAS_SIGNALING_UP",
         "TGF_BETA", "TNFA_SIGNALING", "ANTIGEN"))]
A(f"{'pool':<14} {'pathway':<46} {'n':>3} {'rho(A)':>9} {'rho(A-PPARG)':>13} {'delta':>9} {'q(A)':>9}")
for r in sorted(sel2, key=lambda z: (z["pathway"], z["pool"])):
    A(f"{r['pool']:<14} {r['pathway']:<46} {r['n']:>3} {r['rho_tum']:>+9.4f} {r['rho_tum_noPP']:>+13.4f} {r['d_noPP']:>+9.4f} {r['q_tum']:>9.3g}")

A("")
A("== 4. 关键判读锚点（逐项取值，供正文引用）==")
def get(pool, path):
    for r in rows:
        if r["pool"] == pool and r["pathway"] == path:
            return r
    return None
anchors = [
    ("P1_tumor_all", "KEGG_MTOR_SIGNALING_PATHWAY"),
    ("P1_tumor_all", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
    ("P3_adipo_diff", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
    ("P3_adipo_diff", "HALLMARK_MTORC1_SIGNALING"),
    ("P1_tumor_all", "KEGG_PPAR_SIGNALING_PATHWAY"),
    ("P3_adipo_diff", "KEGG_PPAR_SIGNALING_PATHWAY"),
    ("P1_tumor_all", "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS"),
    ("P2_tumor_prog", "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS"),
    ("P1_tumor_all", "HALLMARK_FATTY_ACID_METABOLISM"),
    ("P1_tumor_all", "KEGG_FATTY_ACID_METABOLISM"),
    ("P1_tumor_all", "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),
    ("P2_tumor_prog", "HALLMARK_KRAS_SIGNALING_UP"),
    ("P1_tumor_all", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    ("P3_adipo_diff", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    ("P4_normal_adi", "HALLMARK_COMPLEMENT"),
    ("P4_normal_adi", "HALLMARK_KRAS_SIGNALING_UP"),
]
for pool, path in anchors:
    r = get(pool, path)
    if r is None:
        A(f"   !! MISSING {pool} / {path}")
        continue
    A(f"   {pool:<14} {path:<46} n={r['n']:>2}  rhoA={r['rho_tum']:+.4f} (P={r['P_tum']:.3g}, q={r['q_tum']:.3g})"
      f"  rhoA-noPP={r['rho_tum_noPP']:+.4f} (P={r['P_tum_noPP']:.3g})"
      f"  full={r['rho_tum_FULL']:+.4f}  d_clean={r['d_clean']:+.4f}")

A("")
A("== 5. full / clean 差异总览 ==")
dc = [r["d_clean"] for r in rows]
A(f"   |d_clean| max={max(abs(x) for x in dc):.4f}  median={sorted(abs(x) for x in dc)[len(dc)//2]:.4f}"
  f"  mean_abs_full={sum(abs(r['rho_tum_FULL']) for r in rows)/len(rows):.4f}"
  f"  mean_abs_clean={sum(abs(r['rho_tum']) for r in rows)/len(rows):.4f}")

txt = "\n".join(L)
with open(OUT, "w", encoding="utf-8") as f:
    f.write(txt + "\n")
print(txt)
print("\nwrote:", OUT)
