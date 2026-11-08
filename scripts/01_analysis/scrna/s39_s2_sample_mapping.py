# -*- coding: utf-8 -*-
"""Step 2 样本映射验证。
一次流式扫描 raw mtx.gz，同时得到：
  (a) 每个细胞的 total counts（用于 CP10K 归一化）
  (b) 一组 marker 基因在每个细胞中的 counts（稀疏提取）
然后按 组织 / 患者 聚合，用脂肪细胞标志（ADIPOQ/PLIN1/LEP）反查样本身份，
验证注释标签与真实生物学一致（映射无错位）。
"""
import os, gzip, json, time
import numpy as np

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
G = os.path.join(W, 'gse221493')
MTX = os.path.join(G, 'GSE221493_sc_LPS_atlas_raw_counts_matrix.mtx.gz')
FEAT = os.path.join(G, 'GSE221493_sc_LPS_atlas_raw_counts_features.tsv.gz')
ANN = os.path.join(G, 'GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz')

OUT_TXT = os.path.join(W, 's39_s2_report.txt')
OUT_JSON = os.path.join(W, 's39_s2_report.json')
OUT_NPZ = os.path.join(W, 's39_s2_marker_matrix.npz')

log = []
def A(s):
    log.append(s)
    print(s, flush=True)

A('=' * 78)
A('Step 2 样本映射验证 — %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
A('=' * 78)

# ---------- 1. features ----------
feat_sym = []
with gzip.open(FEAT, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if p and p[0]:
            feat_sym.append(p[0])
n_features = len(feat_sym)
A('features: %d 行' % n_features)

PANEL = {
    'Adipocyte identity': ['ADIPOQ', 'PLIN1', 'LEP', 'FABP4', 'CIDEC', 'PLIN4', 'PLIN2', 'CEBPA'],
    'Lipid uptake': ['CD36', 'LPL', 'VLDLR', 'FABP3', 'FABP5', 'SLC27A1'],
    'Score-6 (bulk-consistent)': ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN'],
    'Lipogenesis-others': ['INSIG1', 'FADS1', 'FADS2', 'ELOVL6', 'ACACA', 'SCD5'],
    'Mesenchymal / preadipocyte': ['PDGFRA', 'DCN', 'LUM', 'THY1', 'PDPN', 'NT5E'],
    'Immune': ['PTPRC', 'CD3E', 'CD8A', 'CD68', 'CD163', 'LYZ', 'NKG7'],
    'Endothelial / mural': ['PECAM1', 'VWF', 'ACTA2', 'RGS5'],
    'Proliferation': ['MKI67', 'TOP2A'],
    'Housekeeping (technical control)': ['ACTB', 'GAPDH', 'B2M', 'TBP'],
}
wanted = []
for v in PANEL.values():
    wanted.extend(v)
wanted = list(dict.fromkeys(wanted))
sym2rows = {}
for i, s in enumerate(feat_sym):
    sym2rows.setdefault(s, []).append(i)

picked = []      # (symbol, feature_row_0based)
missing = []
for g in wanted:
    if g in sym2rows:
        for r in sym2rows[g]:
            picked.append((g, r))
    else:
        missing.append(g)
A('目标 marker: 命中 %d 个基因（含 %d 个重复名行）；未命中 %d 个' %
  (len(set(g for g, _ in picked)), len(picked), len(missing)))
if missing:
    A('  未命中: %s' % ', '.join(missing))

row2slot = {}
slots = []
sym_list = []
for gi, (g, r) in enumerate(picked):
    if r in row2slot:      # 同一 feature 行已在（理论上不会）
        continue
    row2slot[r] = gi
    slots.append(gi)
    sym_list.append(g)
target_rows = np.array(sorted(row2slot.keys()), dtype=np.int64)
maxrow = int(target_rows.max())
A('唯一目标 feature 行: %d（最大行号 %d / %d）' % (len(target_rows), maxrow, n_features))

rowmap = np.full(n_features, -1, dtype=np.int32)
for r, gi in row2slot.items():
    rowmap[r] = gi

# ---------- 2. 注释表 ----------
barcodes, patients, tissues, ct_l1, ct_l2, ct_l3 = [], [], [], [], [], []
with gzip.open(ANN, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 10:
            continue
        barcodes.append(p[0]); patients.append(p[1]); tissues.append(p[2])
        ct_l1.append(p[5]); ct_l2.append(p[6]); ct_l3.append(p[7])
n_cells = len(barcodes)
A('注释表: %d 细胞' % n_cells)
tissue_arr = np.array(tissues)
patient_arr = np.array(patients)
l1_arr = np.array(ct_l1)
l2_arr = np.array(ct_l2)
l3_arr = np.array(ct_l3)

# ---------- 3. 一次流式扫描 mtx ----------
A('--- 流式扫描 mtx（%.1f MB）...' % (os.path.getsize(MTX) / 1048576))
expr = np.zeros((len(sym_list), n_cells), dtype=np.float32)
total = np.zeros(n_cells, dtype=np.float64)
nnz_target = 0

def consume(block):
    global nnz_target
    a = np.fromstring(block, dtype=np.float64, sep=' ')
    if a.size % 3 != 0:
        raise RuntimeError('block not multiple of 3: %d' % a.size)
    a = a.reshape(-1, 3)
    r0 = a[:, 0].astype(np.int64) - 1          # ★ mtx 行号是 1-based → 必须转 0-based
    c = a[:, 1].astype(np.int64) - 1
    v = a[:, 2].astype(np.float64)
    total[:] += np.bincount(c, weights=v, minlength=n_cells)[:n_cells]
    m = r0 <= maxrow
    if m.any():
        rr = r0[m]
        keep = rowmap[rr] >= 0
        if keep.any():
            gidx = rowmap[rr[keep]].astype(np.int64)
            cidx = c[m][keep]
            vals = v[m][keep]
            np.add.at(expr, (gidx, cidx), vals.astype(np.float32))
            nnz_target += int(keep.sum())

t0 = time.time()
with gzip.open(MTX, 'rb') as f:
    hdr1 = f.readline()
    hdr2 = f.readline()
    A('  header1: %s' % hdr1.decode('utf-8', 'replace').strip())
    A('  header2: %s' % hdr2.decode('utf-8', 'replace').strip())
    buf = b''
    nblk = 0
    while True:
        chunk = f.read(1 << 23)
        if not chunk:
            break
        buf += chunk
        nl = buf.rfind(b'\n')
        if nl < 0:
            continue
        consume(buf[:nl])
        buf = buf[nl + 1:]
        nblk += 1
        if nblk % 20 == 0:
            A('    blk %-4d  %.1f%%  %.0fs' % (nblk, 100 * f.tell() / os.path.getsize(MTX), time.time() - t0))
    if buf.strip():
        consume(buf)
A('  扫描完成: %d 块, %.0fs, 目标基因非零 %d' % (nblk, time.time() - t0, nnz_target))

# 校验：per-cell 深度是否合理
A('  per-cell total counts: median=%.0f  min=%.0f  max=%.0f  sum=%.0f' %
  (np.median(total), total.min(), total.max(), total.sum()))

# ---------- 4. 归一化 ----------
scale = 1e4 / np.maximum(total, 1.0)
logn = np.log1p(expr.astype(np.float64) * scale[None, :])   # slot × cell
sym_index = {g: i for i, g in enumerate(sym_list)}

# --- sanity：行号映射正确性自检（ACTB/B2M 应接近全细胞检出）---
A('')
A('--- sanity: marker 检出率（预期 ACTB/B2M ≈ 全细胞；ADIPOQ 仅脂肪来源高）---')
for g in ['ACTB', 'GAPDH', 'B2M', 'MALAT1', 'ADIPOQ', 'PLIN1', 'LEP', 'FABP4', 'CD36', 'PPARG']:
    if g in sym_index:
        i = sym_index[g]
        A('  %-8s nnz=%-8d detect=%6.1f%%  mean_log1p=%.3f' %
          (g, int((expr[i] > 0).sum()), 100.0 * (expr[i] > 0).mean(), float(logn[i].mean())))

def agg(mask, label):
    if mask.sum() == 0:
        return None
    sub = logn[:, mask]
    det = (expr[:, mask] > 0).mean(axis=1)
    return {'n_cells': int(mask.sum()),
            'mean_log1p_cp10k': {sym_list[i]: round(float(sub[i].mean()), 4) for i in range(len(sym_list))},
            'pct_detected': {sym_list[i]: round(float(det[i] * 100), 2) for i in range(len(sym_list))}}

# 4a. 按组织
A('')
A('=' * 78)
A('4a. 按组织层级聚合（关键 marker）')
A('=' * 78)
KEY = ['ADIPOQ', 'PLIN1', 'LEP', 'FABP4', 'CIDEC', 'CD36', 'LPL',
       'SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN',
       'PDGFRA', 'PTPRC', 'PECAM1', 'ACTB']
order = ['nl', 'lipoma', 'WDLPS', 'DDLPS-WD', 'DDLPS-sWD', 'DDLPS-DD', 'DDLPS-DDK']
present = [t for t in order if (tissue_arr == t).any()]
present += [t for t in sorted(set(tissues)) if t not in present]
A('%-12s %8s | %s' % ('tissue', 'n_cells', ' '.join('%9s' % g for g in KEY)))
tissue_stats = {}
for t in present:
    m = tissue_arr == t
    st = agg(m, t)
    tissue_stats[t] = st
    if st:
        A('%-12s %8d | %s' % (t, st['n_cells'],
                              ' '.join('%9.3f' % st['mean_log1p_cp10k'][g] for g in KEY)))
A('')
A('  检测率 pct_detected (%)')
A('%-12s %8s | %s' % ('tissue', 'n_cells', ' '.join('%9s' % g for g in KEY)))
for t in present:
    st = tissue_stats.get(t)
    if st:
        A('%-12s %8d | %s' % (t, st['n_cells'],
                              ' '.join('%9.2f' % st['pct_detected'][g] for g in KEY)))

# 4b. 按患者×组织（脂肪 marker 三联）
A('')
A('=' * 78)
A('4b. 按 患者×组织（ADIPOQ / PLIN1 / LEP 均值）')
A('=' * 78)
A('%-6s %-12s %8s | %9s %9s %9s' % ('pat', 'tissue', 'n', 'ADIPOQ', 'PLIN1', 'LEP'))
per_sample = {}
for p in sorted(set(patients), key=lambda x: int(x[1:])):
    for t in present:
        m = (patient_arr == p) & (tissue_arr == t)
        if m.sum() == 0:
            continue
        st = agg(m, t)
        per_sample['%s|%s' % (p, t)] = st
        A('%-6s %-12s %8d | %9.3f %9.3f %9.3f' % (
            p, t, st['n_cells'],
            st['mean_log1p_cp10k']['ADIPOQ'], st['mean_log1p_cp10k']['PLIN1'],
            st['mean_log1p_cp10k']['LEP']))

# 4c. 与数据自带 L3 注释交叉（第三路验证）
A('')
A('=' * 78)
A('4c. 与数据自带 L3 细胞注释交叉验证（ADIPOQ 是否落在 adipocyte 类）')
A('=' * 78)
l3_levels = sorted(set(ct_l3))
A('%-34s %8s | %9s %9s %9s' % ('L3 annotation', 'n', 'ADIPOQ', 'PLIN1', 'FABP4'))
l3_rows = []
for l3 in l3_levels:
    m = l3_arr == l3
    if m.sum() < 20:
        continue
    st = agg(m, l3)
    l3_rows.append((l3, st['n_cells'], st['mean_log1p_cp10k']['ADIPOQ'],
                    st['mean_log1p_cp10k']['PLIN1'], st['mean_log1p_cp10k']['FABP4']))
l3_rows.sort(key=lambda x: -x[2])
for l3, n, a, pl, f4 in l3_rows:
    A('%-34s %8d | %9.3f %9.3f %9.3f' % (l3, n, a, pl, f4))

# ---------- 5. 判决 ----------
A('')
A('=' * 78)
A('5. 映射判决')
A('=' * 78)
def msym(t, g):
    st = tissue_stats.get(t)
    return st['mean_log1p_cp10k'][g] if st else None

fat_marks = ['ADIPOQ', 'PLIN1', 'LEP']
verdict = {}
for g in fat_marks:
    vals = {t: msym(t, g) for t in present}
    dd = vals.get('DDLPS-DD')
    nl = vals.get('nl')
    wd = vals.get('DDLPS-WD')
    wdl = vals.get('WDLPS')
    lip = vals.get('lipoma')
    fatmax = max(v for k, v in vals.items() if k in ('nl', 'lipoma') and v is not None) if any(vals.get(k) is not None for k in ('nl', 'lipoma')) else None
    verdict[g] = {'by_tissue': vals, 'DDLPS-DD': dd, 'nl': nl, 'WDLPS': wdl, 'DDLPS-WD': wd}
    ok = (dd is not None and fatmax is not None and fatmax > dd)
    A('%-8s nl=%.3f lipoma=%.3f WDLPS=%.3f WD=%.3f DD=%.3f sWD=%s DDK=%s  → 脂肪组织 > DD ? %s' % (
        g,
        vals.get('nl') if vals.get('nl') is not None else -1,
        vals.get('lipoma') if vals.get('lipoma') is not None else -1,
        vals.get('WDLPS') if vals.get('WDLPS') is not None else -1,
        vals.get('DDLPS-WD') if vals.get('DDLPS-WD') is not None else -1,
        dd if dd is not None else -1,
        ('%.3f' % vals['DDLPS-sWD']) if vals.get('DDLPS-sWD') is not None else '-',
        ('%.3f' % vals['DDLPS-DDK']) if vals.get('DDLPS-DDK') is not None else '-',
        'YES' if ok else 'NO'))

# 技术对照：housekeeping 应跨组织平坦
A('')
A('  技术对照（housekeeping 跨组织变异，用于确认差异非测序深度伪影）：')
for g in ['ACTB', 'GAPDH', 'B2M']:
    vals = [msym(t, g) for t in present if msym(t, g) is not None]
    if vals:
        A('    %-6s 跨组织 mean=%.3f  min=%.3f  max=%.3f  极差=%.3f' %
          (g, np.mean(vals), min(vals), max(vals), max(vals) - min(vals)))

# ---------- 6. 落盘 ----------
np.savez_compressed(OUT_NPZ,
                    symbols=np.array(sym_list),
                    expr=expr,
                    total=total,
                    patients=patient_arr,
                    tissues=tissue_arr,
                    l3=l3_arr)
A('')
A('已保存 marker 矩阵 → %s (%.1f MB)' % (os.path.basename(OUT_NPZ), os.path.getsize(OUT_NPZ) / 1048576))

rep = {'date': time.strftime('%Y-%m-%d'),
       'n_cells': n_cells, 'n_genes_panel': len(sym_list),
       'panel': PANEL, 'missing_genes': missing,
       'per_cell_total_median': float(np.median(total)),
       'tissue_stats': {k: v for k, v in tissue_stats.items() if v},
       'per_sample': per_sample,
       'l3_adipocyte_ranking': [{'l3': a, 'n': b, 'ADIPOQ': c} for a, b, c, _, _ in l3_rows],
       'verdict_fat_vs_dd': verdict}
json.dump(rep, open(OUT_JSON, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
open(OUT_TXT, 'w', encoding='utf-8').write('\n'.join(log))
A('报告 → %s' % OUT_TXT)
