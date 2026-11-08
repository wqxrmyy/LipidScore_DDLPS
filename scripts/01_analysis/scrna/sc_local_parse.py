# -*- coding: utf-8 -*-
"""本地解析已下载文件（不联网）：features 基因核验 + raw barcodes 后缀分布。"""
import io, os, gzip

D = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "gse221493")
out = []


def A(s=""):
    out.append(str(s))


def gz_ok(p):
    try:
        with gzip.open(p, "rb") as g:
            while g.read(1 << 20):
                pass
        return True
    except Exception as e:
        return "BAD: %s" % type(e).__name__


def lines(p):
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        return [l.rstrip("\n") for l in f if l.strip()]


for fn in sorted(os.listdir(D)):
    p = os.path.join(D, fn)
    A("%-64s %10d B  完整性=%s" % (fn, os.path.getsize(p), gz_ok(p)))

# ---- features ----
fp = os.path.join(D, "GSE221493_sc_LPS_atlas_raw_counts_features.tsv.gz")
if os.path.exists(fp):
    A("")
    A("=" * 70)
    A("features 基因核验")
    A("=" * 70)
    fl = lines(fp)
    A("行数: %d | 首行: %s" % (len(fl), fl[0]))
    syms = set(l.split("\t")[0] for l in fl)
    A("唯一基因符号: %d" % len(syms))
    groups = {
        "评分 6 基因": ["SREBF1", "SCD", "HMGCR", "ACLY", "PPARG", "FASN"],
        "脂肪细胞": ["ADIPOQ", "PLIN1", "LEP", "FABP4", "LPL", "CIDEC", "PLIN4"],
        "脂质摄取": ["CD36", "LPL", "VLDLR", "FABP3"],
        "间充质/前体": ["PDGFRA", "DCN", "LUM", "THY1", "PDPN", "NT5E"],
        "免疫": ["PTPRC", "CD3E", "CD8A", "CD68", "CD163", "LYZ", "CD19", "NKG7"],
        "内皮/平滑肌": ["PECAM1", "VWF", "ACTA2", "RGS5", "MCAM"],
        "增殖": ["MKI67", "TOP2A"],
        "12q13-15": ["MDM2", "CDK4", "HMGA2", "FRS2", "YEATS4"],
    }
    for g, genes in groups.items():
        hit = [x for x in genes if x in syms]
        miss = [x for x in genes if x not in syms]
        A("  %-12s 命中 %-52s %s" % (g, ",".join(hit), ("缺: " + ",".join(miss)) if miss else ""))

# ---- raw barcodes ----
bp = os.path.join(D, "GSE221493_sc_LPS_atlas_raw_counts_barcodes.tsv.gz")
if os.path.exists(bp):
    A("")
    A("=" * 70)
    A("raw barcodes 后缀分布")
    A("=" * 70)
    bl = lines(bp)
    cnt = {}
    for l in bl:
        if "-" in l:
            k = int(l.rsplit("-", 1)[1]) if l.rsplit("-", 1)[1].isdigit() else l.rsplit("-", 1)[1]
            cnt[k] = cnt.get(k, 0) + 1
    A("barcode 总数: %d | 不同后缀: %d" % (len(bl), len(cnt)))
    A("后缀 -> 细胞数:")
    for k in sorted(cnt, key=lambda x: (isinstance(x, str), x)):
        A("   -%-4s %6d" % (k, cnt[k]))

dst = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "sc_local_parse.txt")
with io.open(dst, "w", encoding="utf-8") as f:
    f.write("\n".join(out))
print("\n".join(out))
