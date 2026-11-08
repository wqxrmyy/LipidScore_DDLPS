# -*- coding: utf-8 -*-
"""6B 出图数据：把 pathway_vs_score.csv 转成 通路 x 池 宽表（rho / P / q），供 P5 热图直接用。
   避免 P5 再从长表二次聚合导致口径分叉。"""
import csv, os

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
FIG = os.path.join(W, "p5_figures", "fig6_data")
SRC = os.path.join(FIG, "pathway_vs_score.csv")
assert os.path.exists(SRC)

POOLS = ["P1_tumor_all", "P2_tumor_prog", "P3_adipo_diff", "P4_normal_adi"]
rows = []
with open(SRC, "r", encoding="utf-8-sig", newline="") as f:
    for r in csv.DictReader(f):
        rows.append(r)
assert len(rows) == 260, len(rows)

paths = sorted({r["pathway"] for r in rows})
assert len(paths) == 65, len(paths)

def pivot(val, fmt="%.4f"):
    out = []
    for p in paths:
        d = {r["pool"]: r for r in rows if r["pathway"] == p}
        assert set(d) == set(POOLS), (p, sorted(d))
        out.append([p] + [fmt % float(d[q][val]) for q in POOLS])
    return out

def write(name, val, fmt, header):
    p = os.path.join(FIG, name)
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pathway"] + POOLS)
        w.writerows(pivot(val, fmt))
    print("wrote", name, os.path.getsize(p), "bytes")

write("fig6b_wide_rho_clean.csv", "rho_tum", "%.4f", None)
write("fig6b_wide_P_clean.csv", "P_tum", "%.3g", None)
write("fig6b_wide_q_clean.csv", "q_tum", "%.3g", None)
write("fig6b_wide_rho_full.csv", "rho_tum_FULL", "%.4f", None)
write("fig6b_wide_rho_noPPARG.csv", "rho_tum_noPP", "%.4f", None)

# 推荐子集（mTOR/PPAR/脂肪酸/炎症/关键对照）
KEY = ("MTOR|PPAR|FATTY_ACID|UNSATURATED|ADIPOCYTO|INSULIN|OXIDATIVE|CHOLESTEROL|"
       "INFLAMMATORY|INTERFERON|COMPLEMENT|ALLOGRAFT|IL6_JAK|KRAS_SIGNALING_UP|TGF_BETA|"
       "TNFA_SIGNALING|ADIPOGENESIS|MYC_TARGETS|DNA_REPAIR|GLYCOLYSIS|PEROXISOME")
sel = [p for p in paths if any(k in p for k in KEY.split("|"))]
print("recommended subset:", len(sel), "of", len(paths))
with open(os.path.join(FIG, "fig6b_recommended_pathways.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(sel) + "\n")

# 自检：宽表行数
for n in ("fig6b_wide_rho_clean.csv",):
    with open(os.path.join(FIG, n), "r", encoding="utf-8-sig", newline="") as f:
        rr = list(csv.reader(f))
    assert len(rr) == 66 and len(rr[0]) == 5, (len(rr), len(rr[0]))
print("OK: wide tables 66 rows (65 pathways + header) x 5 cols")
