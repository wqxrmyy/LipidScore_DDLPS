# -*- coding: utf-8 -*-
"""
#39 Step 1：raw counts 矩阵（550,902,472 B）分块并行下载
理由：单连接整体 GET 会被拖到 ~10 KB/s 并静默截断；实测 Range 取数可达 ~580 KB/s。
策略：4 MB 分块 + 每块独立 Range + 多线程 + 逐块落地（断点可续）→ 顺序合并 → gzip 全量校验 → 维度断言 → SHA256
输出：_work/s39_s1_matrix_report.json / .txt；日志 _work/s39_s1_dl.log
"""
import os, io, json, gzip, time, hashlib, urllib.request, threading, datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
G = os.path.join(W, "gse221493")
CH = os.path.join(G, "_chunks")
os.makedirs(CH, exist_ok=True)
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221493/suppl/"
FN = "GSE221493_sc_LPS_atlas_raw_counts_matrix.mtx.gz"
URL = BASE + FN
DST = os.path.join(G, FN)
UA = {"User-Agent": "Mozilla/5.0 (compatible; academic-research)"}
TODAY = datetime.date.today().isoformat()

CHUNK = 4 << 20          # 4 MB
WORKERS = 4
TRIES = 12
TIMEOUT = 100

TOTAL = 550902472

_log_lock = threading.Lock()
lf = io.open(os.path.join(W, "s39_s1_dl.log"), "w", encoding="utf-8")
def L(s):
    with _log_lock:
        lf.write("%s %s\n" % (time.strftime("%H:%M:%S"), s)); lf.flush()

def head_size(url, tries=6):
    for k in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers=UA, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                cl = r.headers.get("Content-Length")
                if cl:
                    return int(cl)
        except Exception as e:
            L("HEAD retry %d: %s" % (k, str(e)[:60]))
        time.sleep(1 + k)
    return None

def fetch_window(idx):
    """下载第 idx 块；返回 (idx, 状态字符串, 耗时)"""
    start = idx * CHUNK
    end = min(start + CHUNK, TOTAL) - 1
    want = end - start + 1
    dst = os.path.join(CH, "c%05d.bin" % idx)
    if os.path.exists(dst) and os.path.getsize(dst) == want:
        return idx, "cached", 0.0
    t0 = time.time()
    for k in range(1, TRIES + 1):
        h = dict(UA); h["Range"] = "bytes=%d-%d" % (start, end)
        try:
            req = urllib.request.Request(URL, headers=h)
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                if getattr(r, "status", 200) != 206:
                    raise ValueError("status %s != 206" % getattr(r, "status", "?"))
                buf = bytearray()
                while len(buf) < want:
                    c = r.read(min(1 << 16, want - len(buf)))
                    if not c:
                        break
                    buf += c
            if len(buf) != want:
                raise ValueError("窗口长度 %d != %d" % (len(buf), want))
            with open(dst + ".t", "wb") as f:
                f.write(buf)
            os.replace(dst + ".t", dst)
            return idx, "ok(try%d)" % k, time.time() - t0
        except Exception as e:
            L("  块 %d 重试 %d/%d: %s: %s" % (idx, k, TRIES, type(e).__name__, str(e)[:60]))
            time.sleep(min(1 + k, 10))
    return idx, "FAIL", time.time() - t0

def main():
    L("=" * 70)
    L("分块并行下载 raw counts 矩阵  [%s]" % TODAY)
    L("=" * 70)
    exp = head_size(URL)
    L("HEAD Content-Length = %s" % exp)
    global TOTAL
    if exp:
        TOTAL = exp
    nchunk = (TOTAL + CHUNK - 1) // CHUNK
    L("分块数 = %d（每块 %d B，%d 线程）" % (nchunk, CHUNK, WORKERS))
    L("")

    done = 0
    t_start = time.time()
    fails = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch_window, i): i for i in range(nchunk)}
        for fu in as_completed(futs):
            idx, st, dt = fu.result()
            done += 1
            if st == "FAIL":
                fails.append(idx)
                L("[%3d/%3d] 块 %d %s" % (done, nchunk, idx, st))
            elif st == "cached":
                pass
            else:
                el = time.time() - t_start
                rate = (done * CHUNK / 1048576.0) / max(el, 1e-9)
                L("[%3d/%3d] 块 %d %s %.1fs | 累计 %.1f MB | %.0f KB/s"
                  % (done, nchunk, idx, st, dt, done * CHUNK / 1048576.0, rate * 1024))

    L("")
    L("分块完成：%d 块，失败 %d 块 %s" % (nchunk, len(fails), fails[:20]))
    if fails:
        L("!! 存在失败块，未合并。可重跑本脚本续传。")
        return 1

    # ---- 合并
    L("按序合并 -> %s" % DST)
    tmp = DST + ".part"
    with open(tmp, "wb") as out:
        for i in range(nchunk):
            p = os.path.join(CH, "c%05d.bin" % i)
            with open(p, "rb") as f:
                while True:
                    b = f.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
    if os.path.getsize(tmp) != TOTAL:
        L("!! 合并后长度 %d != %d" % (os.path.getsize(tmp), TOTAL))
        return 1
    os.replace(tmp, DST)
    L("合并完成 %d B (%.1f MB)" % (TOTAL, TOTAL / 1048576.0))

    # ---- 校验
    L("")
    L("=" * 70); L("校验"); L("=" * 70)
    L("① 全量 gzip 解压…")
    t0 = time.time(); nbytes = 0; ok = True; info = ""
    try:
        with gzip.open(DST, "rb") as g:
            while True:
                b = g.read(1 << 24)
                if not b:
                    break
                nbytes += len(b)
    except Exception as e:
        ok = False; info = "%s: %s" % (type(e).__name__, str(e))
    L("   解压 %s：%s（%.1fs）" % ("成功" if ok else "失败",
                                  ("%d B" % nbytes) if ok else info, time.time() - t0))

    dims = None
    if ok:
        with gzip.open(DST, "rt", encoding="utf-8", errors="replace") as g:
            hdr = [g.readline().rstrip("\n") for _ in range(6)]
        L("② MatrixMarket 头部：")
        for i, h in enumerate(hdr[:4]):
            L("   L%d: %s" % (i + 1, h[:110]))
        for h in hdr:
            if h and not h.startswith("%"):
                p = h.split()
                if len(p) >= 3:
                    dims = tuple(int(x) for x in p[:3])
                break
        L("③ 维度 = %s" % (dims,))
        L("   断言 features(行)==49632 : %s" % ("OK" if dims and dims[0] == 49632 else "FAIL"))
        L("   断言 cells(列)==102753   : %s" % ("OK" if dims and dims[1] == 102753 else "FAIL"))

    L("④ SHA256…")
    t0 = time.time()
    h = hashlib.sha256()
    with open(DST, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    sha = h.hexdigest()
    L("   sha256 = %s  (%.1fs)" % (sha, time.time() - t0))

    rep = {"date": TODAY, "file": FN, "bytes": os.path.getsize(DST),
           "head_content_length": TOTAL, "decompressed_bytes": nbytes if ok else None,
           "gzip_ok": ok, "dims": dims, "sha256": sha,
           "assert_features_49632": bool(dims and dims[0] == 49632),
           "assert_cells_102753": bool(dims and dims[1] == 102753),
           "chunks": nchunk, "chunk_bytes": CHUNK, "workers": WORKERS,
           "elapsed_sec": round(time.time() - t_start, 1)}
    with io.open(os.path.join(W, "s39_s1_matrix_report.json"), "w", encoding="utf-8") as f:
        json.dump(rep, f, ensure_ascii=False, indent=2)
    with io.open(os.path.join(W, "s39_s1_matrix_report.txt"), "w", encoding="utf-8") as f:
        f.write(open(os.path.join(W, "s39_s1_dl.log"), encoding="utf-8").read())
    L("")
    L("STEP1-MATRIX DONE  gzip_ok=%s dims=%s" % (ok, dims))
    return 0

rc = main()
lf.close()
print("rc =", rc)
