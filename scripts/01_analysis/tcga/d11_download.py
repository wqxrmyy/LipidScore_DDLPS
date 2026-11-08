# -*- coding: utf-8 -*-
"""D11：下载 GSE221492 bulk LPS atlas 矩阵（独立队列复现用）"""
import os, sys, urllib.request, json, time

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "gse221492")
os.makedirs(OUT, exist_ok=True)
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221492/suppl/'
FILES = [
    'GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz',
    'GSE221492_bulk_LPS_atlas_log2_FPKM_counts_matrix.tsv.gz',
]
rep = {}
for f in FILES:
    dst = os.path.join(OUT, f)
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        rep[f] = {'status': 'cached', 'bytes': os.path.getsize(dst)}
        continue
    t0 = time.time()
    try:
        req = urllib.request.Request(BASE + f, headers={'User-Agent': 'Mozilla/5.0 (research)'})
        with urllib.request.urlopen(req, timeout=180) as r, open(dst, 'wb') as w:
            n = 0
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                w.write(b); n += len(b)
        rep[f] = {'status': 'ok', 'bytes': n, 'sec': round(time.time() - t0, 1)}
    except Exception as e:
        rep[f] = {'status': 'FAIL', 'err': repr(e)[:300]}
    print(f, rep[f], flush=True)

json.dump(rep, open(os.path.join(OUT, '_download.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
