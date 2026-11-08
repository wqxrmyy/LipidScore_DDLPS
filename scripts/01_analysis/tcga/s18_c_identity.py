# -*- coding: utf-8 -*-
"""S18 补充二：用【正确】FASN 重做 d15e 的 C-identity 判定
   问题：C 是否等于"6 基因跨样本 z 分均值"？
"""
import os, csv, json
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

C = {}
for r in csv.DictReader(open(os.path.join(RES, '02_feature_scores.csv'), encoding='utf-8-sig')):
    k = r.get('') or r.get('Sample')
    try: C[str(k).strip()] = float(r['Lipid_Metabolism_Score'])
    except Exception: pass
A = {r['Sample']: float(r['Score']) for r in
     csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}

big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
bidx = {str(k).split('.')[0]: str(k) for k in big.index}

GENES = ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN']
IDS_OK = {'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
          'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710'}
IDS_BAD = dict(IDS_OK); IDS_BAD['FASN'] = 'ENSG00000071082'   # 原脚本误用的 RPL31

samples = sorted(set(A) & set(C) & set(str(c) for c in big.columns))
rep['n_samples'] = len(samples)
raw = {}
for g, e in IDS_OK.items():
    row = bidx[e]
    raw[g] = np.array([float(big.loc[row, s]) for s in samples], float)

a = np.array([A[s] for s in samples], float)
c = np.array([C[s] for s in samples], float)

def zmean(genedict, ref=None):
    """跨样本 z 分后取基因均值。ref=None 表示在本样本集内标准化"""
    M = np.column_stack([genedict[g] for g in GENES])
    if ref is None:
        Z = (M - M.mean(0)) / M.std(0, ddof=1)
    else:
        Z = (M - ref[0]) / ref[1]
    return Z.mean(1)

print('=== C 与「6 基因等权 z 分均值」的一致性 ===')
for lab, d in (('FASN 正确(ENSG00000169710)', IDS_OK), ('FASN 错误=RPL31(ENSG00000071082)', IDS_BAD)):
    gd = {}
    for g, e in d.items():
        row = bidx[e]
        gd[g] = np.array([float(big.loc[row, s]) for s in samples], float)
    s6 = zmean(gd)
    r, p = stats.pearsonr(s6, c)
    rho, _ = stats.spearmanr(s6, c)
    ra, _ = stats.pearsonr(s6, a)
    rep[lab] = {'C_vs_6genezmean_pearson': round(float(r), 6), 'p': float(p),
                'spearman': round(float(rho), 6), 'A_vs_6genezmean_pearson': round(float(ra), 8),
                'max_abs_diff_vs_C': round(float(np.max(np.abs(s6 - c))), 6)}
    print(f'  {lab:<36} C~zmean r={r:.6f} (rho={rho:.6f}) | A~zmean r={ra:.8f}')

# 直接用 C 与各基因的相关（正确 FASN）作为"C 的构造指纹"
print()
print('=== C 与各基因的相关（正确 ID）===')
for g in GENES + ['CD36', 'FABP4']:
    e = {'CD36': 'ENSG00000135218', 'FABP4': 'ENSG00000170323'}.get(g, IDS_OK.get(g))
    row = bidx[e]
    v = np.array([float(big.loc[row, s]) for s in samples], float)
    r, p = stats.pearsonr(v, c)
    ra, _ = stats.pearsonr(v, a)
    rep[f'C~{g}'] = {'r': round(float(r), 4), 'p': round(float(p), 6), 'r_vs_A': round(float(ra), 4)}
    print(f'  C~{g:<7} r={r:+.4f} (P={p:.2e})   |  A~{g} r={ra:+.4f}')

# 用正确 6 基因做"C 对分量的回归"，看权重
M = np.column_stack([raw[g] for g in GENES])
M1 = np.column_stack([np.ones(len(samples)), M])
beta, *_ = np.linalg.lstsq(M1, c, rcond=None)
pred = M1 @ beta
r2 = 1 - np.sum((c - pred) ** 2) / np.sum((c - c.mean()) ** 2)
rep['C_on_6genes_regression'] = {'r2': round(float(r2), 4),
                                 'intercept': round(float(beta[0]), 4),
                                 'coef': {g: round(float(b), 4) for g, b in zip(GENES, beta[1:])}}
print()
print('=== C 对 6 基因的多元回归（看权重是否近似等权）===')
print(f"  R2 = {r2:.4f}")
for g, b in zip(GENES, beta[1:]):
    print(f'    {g:<8} coef = {b:+.4f}')

json.dump(rep, open(os.path.join(OUT, 's18_c_identity.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print()
print('saved -> s18_c_identity.json')
