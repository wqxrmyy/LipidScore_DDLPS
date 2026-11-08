# -*- coding: utf-8 -*-
"""
S34 · Fig6 六面板出图数据准备（冻结件）
--------------------------------------------------------------
原则：
  1) 只读既有产物，不重算任何统计量（ρ/OR/r 一律取自原始输出）
  2) 唯一例外：95%CI 由 ρ 与 n 用 Fisher z 变换现算（公式公开、可复算）
     —— 并与交付件 §5.1 已列值交叉核对，不一致即 assert 失败
  3) 输出到 fig6_data/（冻结目录），CSV 一律 utf-8-sig
来源：
  6A 层读数 : s33_evidence_pack/fig6_stage1_results.md（阶段一，gene symbol 口径）
              + 交付件 §2（bulk 三重校正，D27 修正 ID 后）
  6C LR 对  : s34_6c_report.txt（F/G/E/H 段解析）
  6D 样本级 : s33_6d2_sample_coupling.csv / s33_6d2_coupling_corr.csv（复制入冻结目录）
  6E 配对   : fig6_data/slingshot_pseudotime_V1.csv + cell_meta_tumor.csv（患者映射）
"""
import io, os, csv, re, shutil, math, hashlib, json

W   = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
GD  = os.path.join(W, "p5_figures", "fig6_data")
EP  = os.path.join(W, "s34_evidence_pack")
os.makedirs(GD, exist_ok=True)

log = []
def A(s):
    print(s); log.append(str(s))

def wcsv(path, header, rows):
    with io.open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    A("  wrote %-38s rows=%d" % (os.path.basename(path), len(rows)))

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]

A("=" * 78)
A("S34 Fig6 outplot-data prep")
A("=" * 78)

# ---------------------------------------------------------------- 1. 6D 复制入冻结目录
A("\n[1] 6D sample-level coupling -> fig6_data/ (copy, 冻结)")
for f in ["s33_6d2_sample_coupling.csv", "s33_6d2_coupling_corr.csv"]:
    src = os.path.join(W, f)
    dst = os.path.join(GD, "fig6d_" + f.replace("s33_6d2_", ""))
    shutil.copyfile(src, dst)
    A("  %s -> %s (%d bytes)" % (f, os.path.basename(dst), os.path.getsize(dst)))

# ---------------------------------------------------------------- 2. 6A 五层级读数
A("\n[2] 6A layer ladder (5 layers x 2 genes)")
# 来源逐格标注；metric 说明量纲，供出图分面（各层独立 x 轴，不跨层比较）
A6 = [
    # layer,       label(EN),                       gene,   est,     lo,   hi,   P,        metric
    ("L1_cell",    "single cell",                   "RICTOR", 0.0426, None, None, 3.3e-16,  "rho"),
    ("L1_cell",    "single cell",                   "CD36",   0.0015, None, None, 0.77,     "rho"),
    ("L2_fisher",  "within detected",               "RICTOR", 1.460,  None, None, 3.03e-10, "OR"),
    ("L2_fisher",  "within detected",               "CD36",   1.120,  None, None, 0.138,    "OR"),
    ("L3_scc",     "sample x program",              "RICTOR", 0.2903, None, None, 1.6e-4,   "rho"),
    ("L3_scc",     "sample x program",              "CD36",   0.0721, None, None, 0.36,     "rho"),
    ("L4_sample",  "sample (scRNA)",                "RICTOR", 0.0335, None, None, 0.86,     "rho"),
    ("L4_sample",  "sample (scRNA)",                "CD36",  -0.1405, None, None, 0.45,     "rho"),
    ("L5_bulk_raw","TCGA-SARC (raw)",               "RICTOR", 0.5370, None, None, 5.8e-5,   "r"),
    ("L5_bulk_raw","TCGA-SARC (raw)",               "CD36",   0.3574, None, None, 0.0108,   "r"),
    ("L6_bulk_adj","TCGA-SARC (purity+CN+prolif)",  "RICTOR", 0.3230, None, None, 0.0394,   "partial r"),
    ("L6_bulk_adj","TCGA-SARC (purity+CN+prolif)",  "CD36",   0.3877, None, None, 0.0123,   "partial r"),
]
# 交叉核对（防手抄错）：与交付件 §2 表逐字对照
doc = io.open(os.path.join(EP, "fig6_id_correction_and_N1N3.md"), encoding="utf-8").read()
for token in ["+0.3230 (P=0.0394)", "+0.3877 (P=0.0123)", "+0.5370", "+0.3574"]:
    assert token in doc, "交付件缺少 token: " + token
st1 = io.open(os.path.join(W, "s33_evidence_pack", "fig6_stage1_results.md"), encoding="utf-8").read()
for token in ["OR=1.460, P=3.03e-10", "rho=+0.2903, P=1.6e-4", "rho=+0.0426 (P=3.3e-16)"]:
    assert token in st1, "阶段一缺少 token: " + token
A("  cross-check vs 交付件§2 与 阶段一: OK")
wcsv(os.path.join(GD, "fig6a_layers.csv"),
     ["layer", "label", "gene", "est", "lo", "hi", "P", "metric", "sig"],
     [[a, b, c, d, e, f, g, h, ("yes" if (g is not None and float(g) < 0.05) else "no")]
      for (a, b, c, d, e, f, g, h) in A6])

# ---------------------------------------------------------------- 3. 6C LR pairs（解析报告）
A("\n[3] 6C LR pairs (parse s34_6c_report.txt F/G/E/H 段)")
rep = io.open(os.path.join(W, "s34_6c_report.txt"), encoding="utf-8", errors="replace").read()

def parse_lr(seg_text):
    out = []
    for line in seg_text.splitlines():
        m = re.match(r"^\s{2,}(\S+)\s+([0-9]*\.?[0-9]+)\s*$", line)
        if m:
            out.append((m.group(1), float(m.group(2))))
    return out

def seg(start, end):
    i = rep.find(start)
    j = rep.find(end)
    assert i > 0 and j > i, "报告段未找到: " + start
    return rep[i + len(start):j]

f_seg = seg("-- F. top LR pairs: TAM3 -> tumor programs --", "-- G. top LR pairs: tumor -> TAM3 --")
g_seg = seg("-- G. top LR pairs: tumor -> TAM3 --", "-- H. TAM1 vs TAM3")
F = parse_lr(f_seg); G = parse_lr(g_seg)
A("  parsed F rows=%d  G rows=%d" % (len(F), len(G)))
assert abs(dict(F)["SPP1_CD44"] - 0.0827) < 1e-9, "SPP1_CD44 解析异常"
assert abs(dict(G)["THBS2_CD36"] - 0.0770) < 1e-9, "THBS2_CD36 解析异常"
assert abs(dict(G)["APP_CD74"] - 0.5494) < 1e-9, "APP_CD74 解析异常"
A("  spot-check SPP1_CD44=0.0827 / THBS2_CD36=0.0770 / APP_CD74=0.5494 : OK")

# TAM3 专属判定：SPP1 家族在 TAM1 的对应权重（E 段 SPP1 m2t 中 Invasion 列 = 0.0478 属对方方向）
# 只登记「本图要标出的专属轴」——来源：交付件 §10.3（SPP1 家族 TAM1=0.0000）
EXCL = {"SPP1_CD44", "SPP1_ITGAV_ITGB1", "SPP1_ITGA5_ITGB1", "SPP1_ITGAV_ITGB5", "LGALS9_CD44", "NAMPT_ITGA5_ITGB1"}
rows = []
for name, v in F:
    rows.append(["TAM3->tumour", name, v, ("exclusive" if name in EXCL else "shared")])
for name, v in G:
    rows.append(["tumour->TAM3", name, v, "na"])
wcsv(os.path.join(GD, "fig6c_lr_pairs.csv"), ["dir", "pair", "weight", "flag"], rows)

# hub（H 段）
h_seg = seg("-- H. TAM1 vs TAM3 as SENDER", "-- I. per-pathway")
m = re.search(r"TAM1 out=([\d.]+)\s+in=([\d.]+)\s+TAM3 out=([\d.]+)\s+in=([\d.]+)", h_seg)
m2 = re.search(r"per-cell out: TAM1=([\d.]+)\s+TAM3=([\d.]+)", h_seg)
assert m and m2, "H 段解析失败"
hub = [["TAM1", float(m.group(1)), float(m.group(2)), float(m2.group(1))],
       ["TAM3", float(m.group(3)), float(m.group(4)), float(m2.group(2))]]
assert abs(hub[1][1] - 3.165) < 1e-9 and abs(hub[1][2] - 7.980) < 1e-9, "hub 解析异常"
wcsv(os.path.join(GD, "fig6c_hub.csv"), ["subtype", "outgoing", "incoming", "per_cell_out"], hub)

# 方向对照（E 段 focus pathways，t2m vs m2t，5 个肿瘤程序池）
# ★ 报告排版：通路名只印在 t2m 行，m2t 行首列为空白 → 须沿用上一通路名
e_seg = seg("-- E. tumor <-> TAM : focus pathways", "-- F. top LR pairs:")
pools = ["Adipo diff", "Invasion", "ECM remodeling", "Stemness", "Hypoxia"]
drows = []
last_pw = None
for line in e_seg.splitlines():
    m = re.match(r"^(?:(\S+)\s+)?(t2m|m2t)\s+([0-9.\s]+)$", line.strip())
    if not m:
        continue
    if m.group(1):
        last_pw = m.group(1)
    vals = [float(x) for x in m.group(3).split()]
    if len(vals) != 5:
        raise AssertionError("E 段列数 != 5: %r -> %s" % (line, vals))
    drows.append([last_pw, m.group(2)] + vals)
A("  E-seg parsed rows=%d (通路 %d)" % (len(drows), len(set(r[0] for r in drows))))
assert len(drows) == 34, "E 段解析行数异常: %d" % len(drows)
assert len(set(r[0] for r in drows)) == 17
# 交叉核对报告原文三个关键值（SPP1 m2t Invasion / COLLAGEN t2m ECM remod / TGFB 全 0）
dd = {(r[0], r[1]): r[2:] for r in drows}
assert abs(dd[("SPP1", "m2t")][1] - 0.0478) < 1e-9, "SPP1 m2t 解析异常"
assert abs(dd[("COLLAGEN", "t2m")][2] - 0.7012) < 1e-9, "COLLAGEN t2m 解析异常"
A("  spot-check SPP1|m2t[Invasion]=0.0478 / COLLAGEN|t2m[ECM]=0.7012 : OK")
wcsv(os.path.join(GD, "fig6c_direction.csv"), ["pathway", "dir"] + pools, drows)

# ---------------------------------------------------------------- 4. 6D forest（ρ + CI 现算 + MDE）
A("\n[4] 6D forest (rho + Fisher-z 95%CI 现算 + MDE)")
corr = list(csv.DictReader(io.open(os.path.join(W, "s33_6d2_coupling_corr.csv"), encoding="utf-8-sig")))
def ci(rho, n):
    z = math.atanh(rho); s = 1.0 / math.sqrt(n - 3.0)
    return math.tanh(z - 1.96 * s), math.tanh(z + 1.96 * s)
MDE27, MDE19 = 0.517, 0.605
NAME = {"tumor_score|mac_PPARG":   ("Macrophage PPARG", "ref"),
        "tumor_score|mac_m2like":  ("M2-like score",   "neg"),
        "tumor_score|mac_CD36":    ("Macrophage CD36", "neg"),
        "tumor_score|mac_RICTOR":  ("Macrophage RICTOR", "NS"),
        "tumor_score_noPP|mac_m2like": ("A-PPARG vs M2-like", "sens"),
        "tumor_score_noPP|mac_CD36":   ("A-PPARG vs CD36", "sens")}
frows = []
for r in corr:
    key = r["x"] + "|" + r["y"]
    lab, cls = NAME[key]
    rho = float(r["rho"]); n = int(r["n"]); P = float(r["p"])
    lo, hi = ci(rho, n)
    frows.append([lab, key, round(rho, 4), round(lo, 4), round(hi, 4), P, n, cls, MDE27])
    A("  %-22s rho=%+.4f CI=[%+.3f,%+.3f] P=%.4g n=%d" % (lab, rho, lo, hi, P, n))
# 交叉核对交付件 §5.1 的 CI（防公式/口径漂移）
# 注：交付件排版用 Unicode 减号 U+2212，故断言避开符号位以免误报
assert "[+0.004, +0.666]" in doc and "0.308, +0.448" in doc
chk = [r for r in frows if r[0] == "Macrophage PPARG"][0]
assert abs(chk[2] - 0.3836) < 0.001 and abs(chk[3] - 0.004) < 0.002 and abs(chk[4] - 0.666) < 0.002, \
    "6D CI 现算与交付件不一致: %s" % chk
A("  CI formula cross-check vs 交付件§5.1: OK (Fisher z, n=27)")
wcsv(os.path.join(GD, "fig6d_forest.csv"),
     ["label", "key", "rho", "lo", "hi", "P", "n", "class", "MDE"], frows)

# ---------------------------------------------------------------- 5. 6E 患者级配对伪时序
A("\n[5] 6E per-patient pseudotime (V1-L1)")
meta = {}
with io.open(os.path.join(GD, "cell_meta_tumor.csv"), encoding="utf-8-sig") as f:
    rd = csv.reader(f); hdr = next(rd)
    ib, ip, it = hdr.index("barcode"), hdr.index("patient"), hdr.index("tissue")
    for row in rd:
        meta[row[ib]] = (row[ip], row[it])
A("  cell_meta rows=%d" % len(meta))
agg = {}   # (patient, tissue) -> [pt...]
with io.open(os.path.join(GD, "slingshot_pseudotime_V1.csv"), encoding="utf-8-sig") as f:
    rd = csv.reader(f); hdr = next(rd)
    ic, ipt = hdr.index("cell"), hdr.index("pt_L1")
    ncell = 0
    for row in rd:
        v = row[ipt]
        if v in ("", "NA"):
            continue
        pm = meta.get(row[ic])
        if pm is None:
            continue
        agg.setdefault(pm, []).append(float(v))
        ncell += 1
A("  cells with pt_L1 & patient: %d" % ncell)
# ★ 聚合口径必须与 s34_6e_patientlevel.R 一致：per patient x tissue 用 **mean**（非 median）
erows = []
for (pat, tis), v in sorted(agg.items()):
    erows.append([pat, tis, round(sum(v) / len(v), 4), len(v)])
wcsv(os.path.join(GD, "fig6e_pair_pt.csv"), ["patient", "tissue", "mean_pt", "n_cells"], erows)

# 7 对配对差值（DD − WD）
PAIRS = ["P1", "P4", "P7", "P10", "P14", "P15", "P16"]
d = {}
for pat, tis, m, n in erows:
    if pat in PAIRS:
        d.setdefault(pat, {})[tis] = m
prows = []
for pat in PAIRS:
    if "DDLPS-WD" in d.get(pat, {}) and "DDLPS-DD" in d.get(pat, {}):
        prows.append([pat, d[pat]["DDLPS-WD"], d[pat]["DDLPS-DD"],
                      round(d[pat]["DDLPS-DD"] - d[pat]["DDLPS-WD"], 4)])
A("  paired (DD-WD) n=%d" % len(prows))
for r in prows:
    A("    %-5s WD=%.3f DD=%.3f diff=%+.3f" % tuple(r))
mean_d = sum(r[3] for r in prows) / len(prows)

# ★ Wilcoxon signed-rank，exact（n=7，枚举 2^7 符号组合，无 ties）—— 与交付件 P=0.0313 对账
import itertools
def wilcox_exact(diffs):
    dd = [x for x in diffs if x != 0]
    n = len(dd)
    idx = sorted(range(n), key=lambda i: abs(dd[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(dd[idx[j + 1]]) == abs(dd[idx[i]]):
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[idx[k]] = avg
        i = j + 1
    Wp = sum(ranks[i] for i in range(n) if dd[i] > 0)
    base = sum(ranks)
    cnt = tot = 0
    for signs in itertools.product([0, 1], repeat=n):
        w = sum(ranks[i] for i in range(n) if signs[i])
        tot += 1
        if abs(w - base / 2.0) >= abs(Wp - base / 2.0) - 1e-9:
            cnt += 1
    return Wp, cnt / float(tot)

Wp, Pw = wilcox_exact([r[3] for r in prows])
A("  mean delta = %+.3f | Wilcoxon W+=%.1f exact two-sided P=%.4f" % (mean_d, Wp, Pw))
assert abs(mean_d - 3.573) < 0.02, "6E 配对差值均值 %.3f 与交付件 +3.573 不符" % mean_d
assert abs(Pw - 0.0313) < 0.005, "6E Wilcoxon P=%.4f 与交付件 0.0313 不符" % Pw
A("  cross-check vs 交付件§11.3 (+3.573 / P=0.0313): OK")
for r in prows:
    r.append(round(Pw, 4))
wcsv(os.path.join(GD, "fig6e_pair_delta.csv"), ["patient", "WD", "DD", "delta", "P_wilcox"], prows)
with io.open(os.path.join(GD, "fig6e_pair_meta.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    w.writerow(["mean_delta_DD_minus_WD", round(mean_d, 4)])
    w.writerow(["wilcoxon_P_exact", round(Pw, 4)])
    w.writerow(["n_pairs", len(prows)])
    w.writerow(["lineage", "V1-L1"])

# ---------------------------------------------------------------- 6. 哈希留档
A("\n[6] freeze hash (fig6_data 本次新增/更新件)")
files = sorted(set([f for f in os.listdir(GD) if f.startswith(("fig6a_", "fig6c_", "fig6d_", "fig6e_"))]))
hs = {}
for f in files:
    p = os.path.join(GD, f)
    hs[f] = {"bytes": os.path.getsize(p), "sha256_16": sha(p)}
    A("  %-30s %8d  %s" % (f, hs[f]["bytes"], hs[f]["sha256_16"]))
with io.open(os.path.join(GD, "fig6_outplot_data_hash.json"), "w", encoding="utf-8") as f:
    json.dump(hs, f, ensure_ascii=False, indent=2)
A("\nDONE")

with io.open(os.path.join(W, "s34_f6_prep_report.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log) + "\n")
