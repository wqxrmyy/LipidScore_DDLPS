# -*- coding: utf-8 -*-
"""S18 关键补充：定义评分 A 的 6 基因矩阵，其 'FASN' 行到底是真 FASN 还是 RPL31？
   并解释 CD3E/CD8A 原始相关的 n 差异。
"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

# 6 基因矩阵（定义评分 A 的矩阵）
m6 = pyreadr.read_r(os.path.join(RAW, 'TCGA_6genes_expr.rds'))
m6 = m6[list(m6.keys())[0]]
print('=== TCGA_6genes_expr.rds ===')
print('shape', m6.shape, '| index:', list(m6.index))
print('cols[:3]', list(m6.columns)[:3])

# 全转录组 log2CPM
big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
bidx = {str(k).split('.')[0]: str(k) for k in big.index}

IDS = {'FASN_correct': 'ENSG00000169710', 'RPL31_wrong': 'ENSG00000071082',
       'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194',
       'HMGCR': 'ENSG00000113161', 'ACLY': 'ENSG00000131473',
       'PPARG': 'ENSG00000132170'}
samples = [str(c) for c in m6.columns]
print('n samples in 6gene matrix =', len(samples))

print()
print(f"{'6gene row':<10}{'vs FASN(correct)':>20}{'vs RPL31':>14}{'verdict':>22}")
print('-' * 70)
for row in m6.index:
    v = np.array([float(m6.loc[row, s]) for s in samples], float)
    out = {}
    for lab, e in (('FASN', 'ENSG00000169710'), ('RPL31', 'ENSG00000071082')):
        r2 = bidx.get(e)
        if not r2: continue
        w = np.array([float(big.loc[r2, s]) for s in samples], float)
        rr, pp = stats.pearsonr(v, w)
        out[lab] = (round(float(rr), 5), round(float(pp), 8))
    rep[str(row)] = out
    f = out.get('FASN', ('-', '-'))[0]; rl = out.get('RPL31', ('-', '-'))[0]
    verdict = 'FASN ✔' if (isinstance(f, float) and f > 0.99) else ('RPL31 ✘' if (isinstance(rl, float) and rl > 0.99) else '?')
    print(f'{str(row):<10}{str(f):>20}{str(rl):>14}{verdict:>22}')

# CD3E / CD8A 的 n 差异
print()
print('=== CD8A(正确) vs CD3E(被误标) 的 n 分解 ===')
A = {r['Sample']: float(r['Score']) for r in
     csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
bybar = {r.get('barcode', ''): r for r in cl}
bypat = {}
for r in cl:
    bypat.setdefault((r.get('barcode') or '')[:12], r)

def num(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except Exception: return None

ss = sorted(set(A) & set(str(c) for c in big.columns))
def get(e):
    r = bidx[e]
    return np.array([float(big.loc[r, s]) for s in ss], float)
def pu(s):
    r0 = bybar.get(s) or bypat.get(s[:12]) or {}
    return num(r0.get('paper_purity'))

for lab, e in (('CD8A', 'ENSG00000153563'), ('CD3E', 'ENSG00000198851')):
    v = get(e); a = np.array([A[s] for s in ss], float)
    r50, p50 = stats.pearsonr(v, a)
    m = np.array([pu(s) is not None for s in ss])
    r44, p44 = stats.pearsonr(v[m], a[m])
    rep[f'{lab}_n50'] = {'r': round(float(r50), 4), 'p': round(float(p50), 6), 'n': 50}
    rep[f'{lab}_n44'] = {'r': round(float(r44), 4), 'p': round(float(p44), 6), 'n': int(m.sum())}
    print(f'  {lab}: n=50 r={r50:.4f} P={p50:.6f} | n=44(有纯度) r={r44:.4f} P={p44:.6f}')

json.dump(rep, open(os.path.join(OUT, 's18_matrix_check.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print()
print('saved -> s18_matrix_check.json')
