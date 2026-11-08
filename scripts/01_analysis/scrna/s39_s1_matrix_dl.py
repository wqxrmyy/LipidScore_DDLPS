# -*- coding: utf-8 -*-
"""
#39 Step 1（前半）：下载 raw counts 矩阵（550.9 MB）
- HEAD 取 Content-Length → Range 断点续传
- 全量 gzip 解压校验（CRC32 + ISIZE 由 gzip 自身保证）
- MatrixMarket 头部维度断言：行 == 49632（features）列 == 102753（cells）
- SHA256 + 下载日期留痕
输出：_work/s39_s1_matrix_report.json / .txt
"""
import os, io, json, gzip, time, hashlib, urllib.request, datetime

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
G = os.path.join(W, "gse221493")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221493/suppl/"
FN = "GSE221493_sc_LPS_atlas_raw_counts_matrix.mtx.gz"
UA = {"User-Agent": "Mozilla/5.0 (compatible; academic-research; contact: local)"}
TODAY = datetime.date.today().isoformat()

out = []
def A(s=""):
    out.append(str(s)); print(s, flush=True)

def flush():
    with io.open(os.path.join(W, "s39_s1_matrix_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(out))

def head_size(url, tries=6):
    for k in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers=UA, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                cl = r.headers.get("Content-Length")
                if cl:
                    return int(cl)
        except Exception as e:
            A("  [head retry %d] %s" % (k, str(e)[:60]))
        time.sleep(1 + k)
    return None

def gz_full_ok(path):
    try:
        n = 0
        with gzip.open(path, "rb") as g:
            while True:
                b = g.read(1 << 24)
                if not b:
                    break
                n += len(b)
        return True, n
    except Exception as e:
        return False, "%s: %s" % (type(e).__name__, str(e))

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

A("=" * 74)
A("Step 1 前半：raw counts 矩阵下载与校验   [%s]" % TODAY)
A("=" * 74)

url = BASE + FN
dst = os.path.join(G, FN)
tmp = dst + ".part"
expect = head_size(url)
A("  HEAD Content-Length = %s B (%.1f MB)" % (expect, expect / 1048576.0))
flush()

# 残留 .part 处理
if os.path.exists(tmp) and expect is not None:
    sz = os.path.getsize(tmp)
    if sz >= expect:
        A("  [part] 残留 .part %d B >= 期望 → 校验后认领" % sz)
    else:
        A("  [part] 残留 .part %d B → 续传" % sz)

tries = 40
for k in range(1, tries + 1):
    have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    if expect is not None and have > expect:
        os.remove(tmp); have = 0
    hdr = dict(UA)
    if have:
        hdr["Range"] = "bytes=%d-" % have
    try:
        req = urllib.request.Request(url, headers=hdr)
        with urllib.request.urlopen(req, timeout=300) as r:
            code = getattr(r, "status", 200)
            mode = "ab" if (have and code == 206) else "wb"
            if mode == "wb":
                have = 0
                if expect is None:
                    cl = r.headers.get("Content-Length")
                    if cl:
                        expect = int(cl)
            got = have
            t0 = time.time()
            with open(tmp, mode) as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b); got += len(b)
                    if got % (32 << 20) < (1 << 20):
                        el = time.time() - t0
                        A("      ... %6.1f MB / %5.1f MB  (%.2f MB/s)"
                          % (got / 1048576.0, (expect or 0) / 1048576.0,
                             (got - have) / 1048576.0 / max(el, 0.01)))
                        flush()
        if expect is not None and got != expect:
            raise ValueError("长度不符 %d != %d" % (got, expect))
        os.replace(tmp, dst)
        A("  [ok] 下载完成 %d B（尝试 %d 次）" % (got, k))
        break
    except Exception as e:
        A("  [retry %d/%d] %s: %s" % (k, tries, type(e).__name__, str(e)[:70]))
        if "416" in str(e) or "Range Not Satisfiable" in str(e):
            if os.path.exists(tmp):
                os.remove(tmp)
            A("      → 416：丢弃 .part 从头下载")
        time.sleep(min(3 + 3 * k, 30))
    flush()
else:
    A("  [FAIL] 重试耗尽")

# ---------------- 校验
if os.path.exists(dst):
    A("")
    A("=" * 74)
    A("校验")
    A("=" * 74)
    A("  文件 %s" % dst)
    A("  字节 %d（期望 %s）%s" % (os.path.getsize(dst), expect,
                                 "✓" if expect == os.path.getsize(dst) else "✗"))

    A("  ① 全量 gzip 解压中（约 550 MB → 数 GB）…")
    flush()
    t0 = time.time()
    ok, info = gz_full_ok(dst)
    A("     解压 %s：%s（%.1fs）" % ("成功" if ok else "失败",
                                   ("%d B" % info) if ok else info, time.time() - t0))
    if not ok:
        A("  !! gzip 校验失败，文件不可用")
    else:
        # ② MatrixMarket 头部
        with gzip.open(dst, "rt", encoding="utf-8", errors="replace") as g:
            hdr = []
            for _ in range(6):
                hdr.append(g.readline().rstrip("\n"))
        A("  ② MatrixMarket 头部：")
        for i, h in enumerate(hdr[:4]):
            A("      L%d: %s" % (i + 1, h[:100]))
        data_line = None
        for h in hdr:
            if h and not h.startswith("%"):
                data_line = h; break
        dims = None
        if data_line:
            p = data_line.split()
            if len(p) >= 3:
                dims = tuple(int(x) for x in p[:3])
        A("  ③ 维度 = %s   （矩阵 Market 为 genes × cells）" % (dims,))
        A("     断言 features(行) == 49632 : %s" % ("✓" if dims and dims[0] == 49632 else "✗"))
        A("     断言 cells(列)   == 102753: %s" % ("✓" if dims and dims[1] == 102753 else "✗"))

        A("  ④ SHA256 计算中…")
        flush()
        t0 = time.time()
        h = sha256(dst)
        A("     sha256 = %s" % h)
        A("     耗时 %.1fs" % (time.time() - t0))

        rep = {"date": TODAY, "file": FN, "bytes": os.path.getsize(dst),
               "expected_bytes": expect, "decompressed_bytes": info,
               "dims": dims, "sha256": h,
               "assert_features_49632": bool(dims and dims[0] == 49632),
               "assert_cells_102753": bool(dims and dims[1] == 102753)}
        with io.open(os.path.join(W, "s39_s1_matrix_report.json"), "w", encoding="utf-8") as f:
            json.dump(rep, f, ensure_ascii=False, indent=2)
        A("")
        A("STEP0/1-MATRIX DONE")

flush()
print("written s39_s1_matrix_report.txt")
