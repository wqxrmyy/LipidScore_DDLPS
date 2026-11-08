# -*- coding: utf-8 -*-
"""S18 修正：CD8A（被误用为 CD3E）与 FASN（被误用为 RPL31）的正确数值
输出：s18_corrected.json / s18_corrected.txt
"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {'note': 'S18：修正两处 Ensembl ID 错配后重算'}

# ---------- 三套评分 ----------
A = {r['Sample']: float(r['Score']) for r in
     csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
C = {}
for r in csv.DictReader(open(os.path.join(RES, '02_feature_scores.csv'), encoding='utf-8-sig')):
    k = r.get('') or r.get('Sample')
    try: C[str(k).strip()] = float(r['Lipid_Metabolism_Score'])
    except Exception: pass

# ---------- 临床 ----------
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
print('clinical cols(cn-relevant):',
      [c for c in cl[0].keys() if 'pur' in c.lower() or 'grade' in c.lower()
       or 'mitot' in c.lower() or 'barc' in c.lower() or 'Sample' in c])
bybar = {r.get('barcode', ''): r for r in cl}
bypat = {}
for r in cl:
    b = r.get('barcode', '') or ''
    bypat.setdefault(b[:12], r)

def lookup(sample):
    r = bybar.get(sample)
    if r: return r
    return bypat.get(sample[:12])

def num(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except Exception: return None

# ---------- 表达矩阵 ----------
big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
bidx = {str(k).split('.')[0]: str(k) for k in big.index}
cols = set(str(c) for c in big.columns)

IDS = {
    # label            : (ensembl, 是否为本脚本认定的正确 ID)
    'CD8A_correct': ('ENSG00000153563', True),
    'CD3E_mislabeled_as_CD8A': ('ENSG00000198851', False),
    'FASN_correct': ('ENSG00000169710', True),
    'RPL31_mislabeled_as_FASN': ('ENSG00000071082', False),
}

samples = sorted(set(A) & set(C) & cols)
rep['n_ABC_common'] = len(samples)
print('n samples(A&C&matrix) =', len(samples))

vals = {}
for lab, (e, ok) in IDS.items():
    row = bidx.get(e)
    if not row:
        print('!! NOT FOUND in matrix:', lab, e); continue
    v = []
    for s in samples:
        try: v.append(float(big.loc[row, s]))
        except Exception: v.append(np.nan)
    vals[lab] = np.array(v, float)
    print(f'  {lab:<28} {e}  row={row}')

a = np.array([A[s] for s in samples], float)
c = np.array([C[s] for s in samples], float)
pur = np.array([num((lookup(s) or {}).get('paper_purity')) for s in samples], float)
print('purity available:', int(np.sum(~np.isnan(pur))), '/', len(pur))

def partial(x, y, z):
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    X, Y, Z = x[m], y[m], z[m]
    if len(X) < 6: return None
    def resid(t, u):
        M = np.column_stack([np.ones_like(u), u])
        b, *_ = np.linalg.lstsq(M, t, rcond=None)
        return t - M @ b
    rx, ry = resid(X, Z), resid(Y, Z)
    n = len(X); dfree = n - 3
    if dfree < 1: return None
    r, _ = stats.pearsonr(rx, ry)
    if abs(r) >= 1: return None
    t = r * math.sqrt(dfree / (1 - r * r))
    return {'n': n, 'partial_r': round(float(r), 4),
            'p': round(float(2 * stats.t.sf(abs(t), dfree)), 6)}

print()
print('=' * 96)
print(f"{'gene(variant)':<28}{'vs A r':>10}{'p':>11}{'vs C r':>10}{'p':>11}{'partial|purity (A)':>24}")
print('=' * 96)
for lab, v in vals.items():
    ra, pa = stats.pearsonr(v, a)
    rc, pc = stats.pearsonr(v, c)
    pr = partial(v, a, pur)
    rep[lab] = {
        'ensembl': IDS[lab][0],
        'is_correct_id_for_label': IDS[lab][1],
        'vs_A': {'r': round(float(ra), 4), 'p': round(float(pa), 6), 'n': len(v)},
        'vs_C': {'r': round(float(rc), 4), 'p': round(float(pc), 6), 'n': len(v)},
        'partial_A_ctrl_purity': pr,
    }
    pstr = f"{pr['partial_r']} (P={pr['p']}, n={pr['n']})" if pr else 'n/a'
    print(f'{lab:<28}{ra:>10.4f}{pa:>11.2e}{rc:>10.4f}{pc:>11.2e}{pstr:>24}')

# ---------- FASN(正确) 与 grade / 核分裂 ----------
print()
print('--- FASN_correct 与 grade / 核分裂 ---')
gr, mi = [], []
for s in samples:
    r0 = lookup(s) or {}
    gr.append(num(r0.get('paper_FNCLCC.grade')) if r0.get('paper_FNCLCC.grade') not in (None, '') else None)
    mi.append(num(r0.get('paper_mitotic.rate')) if r0.get('paper_mitotic.rate') not in (None, '') else None)
v = vals['FASN_correct']
for name, arr in (('grade', gr), ('mitotic', mi)):
    x, y = [], []
    for i in range(len(samples)):
        if arr[i] is not None and not np.isnan(v[i]):
            x.append(arr[i]); y.append(v[i])
    if len(x) > 5:
        rs, ps = stats.spearmanr(x, y)
        rep[f'FASN_correct_vs_{name}'] = {'spearman_rho': round(float(rs), 4), 'p': round(float(ps), 6), 'n': len(x)}
        print(f'  FASN vs {name}: rho={rs:.4f}, P={ps:.6f}, n={len(x)}')

json.dump(rep, open(os.path.join(OUT, 's18_corrected.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print()
print('saved -> s18_corrected.json')
