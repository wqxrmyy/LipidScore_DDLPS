# -*- coding: utf-8 -*-
"""卡 D 是否稳健于「评分口径」？
把同一批端点分别在两个评分上算：
  A = TCGA_Score.csv（6 基因 z 分均值，公式已验证；审计链 d7/d11/d12 使用）
  C = 02_feature_scores.csv 的 Lipid_Metabolism_Score（手稿全部图脚本使用；生成脚本缺失）
"""
import os, csv, json, math
import numpy as np
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}


def load(path, keycol, valcol):
    d = {}
    for r in csv.DictReader(open(path, encoding='utf-8-sig')):
        k = r.get(keycol) or r.get('')
        try: d[str(k).strip()] = float(r[valcol])
        except Exception: pass
    return d


def num(x):
    try:
        v = float(str(x).strip()); return None if math.isnan(v) else v
    except Exception: return None


def _f(v):
    out = []
    for x in v:
        try: out.append(float(x))
        except Exception: out.append(np.nan)
    return np.asarray(out, float)


def pear(x, y):
    x, y = _f(x), _f(y); m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5 or len(set(x)) < 3: return None
    r, p = stats.pearsonr(x, y); return {'r': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def spear(x, y):
    x, y = _f(x), _f(y); m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5 or len(set(x)) < 3: return None
    r, p = stats.spearmanr(x, y); return {'rho': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def partial(x, y, *zs):
    X, Y = _f(x), _f(y); Z = np.column_stack([_f(z) for z in zs])
    keep = ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z).any(axis=1))
    X, Y, Z = X[keep], Y[keep], Z[keep]
    if len(X) < 8 + Z.shape[1]: return None
    A_ = np.column_stack([np.ones(len(Z)), Z])
    bx, *_ = np.linalg.lstsq(A_, X, rcond=None); by, *_ = np.linalg.lstsq(A_, Y, rcond=None)
    rx, ry = X - A_ @ bx, Y - A_ @ by
    n = len(X); r, _ = stats.pearsonr(rx, ry); d = n - 2 - Z.shape[1]
    t = r * math.sqrt(d / max(1e-9, 1 - r * r)); p = float(2 * stats.t.sf(abs(t), d))
    return {'r': round(float(r), 4), 'p': round(p, 5), 'n': n, 'n_covars': Z.shape[1]}


A = load(os.path.join(RAW, 'TCGA_Score.csv'), 'Sample', 'Score')
C = load(os.path.join(RES, '02_feature_scores.csv'), '', 'Lipid_Metabolism_Score')
cli = {r['barcode']: r for r in
       csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}

common = sorted(set(A) & set(C) & set(cli))
rep['n_common'] = len(common)

W = []
for s in common:
    r0 = cli[s]
    W.append({'s': s, 'A': A[s], 'C': C[s],
              'cn': num(r0.get('paper_CN.cluster')),
              'grade': num(r0.get('paper_FNCLCC.grade')),
              'mitotic': num(r0.get('paper_mitotic.rate')),
              'purity': num(r0.get('paper_purity')),
              'os_time': num(r0.get('OS_time')) or num(r0.get('paper_OS.days')),
              'os_status': num(r0.get('OS_status'))})

res = {}
for nm in ['A', 'C']:
    d = {}
    d['CN_raw'] = pear([w['cn'] for w in W], [w[nm] for w in W])
    d['CN_ctrl_purity'] = partial([w['cn'] for w in W], [w[nm] for w in W], [w['purity'] for w in W])
    d['CN_ctrl_purity_grade_mitotic'] = partial(
        [w['cn'] for w in W], [w[nm] for w in W],
        [w['purity'] for w in W], [w['grade'] for w in W], [w['mitotic'] for w in W])
    d['grade'] = spear([w['grade'] for w in W], [w[nm] for w in W])
    d['mitotic'] = spear([w['mitotic'] for w in W], [w[nm] for w in W])
    res[nm] = d

# CD36 / CD8A 表达
try:
    import pyreadr
    big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
    big = big[list(big.keys())[0]]
    bidx = {str(k).split('.')[0]: str(k) for k in big.index}
    for g, e in [('CD36', 'ENSG00000135218'), ('CD8A', 'ENSG00000153563'), ('PPARG', 'ENSG00000132170')]:
        row = bidx.get(e)
        if not row: continue
        vals = []
        for w in W:
            try: vals.append(float(big.loc[row, w['s']]))
            except Exception: vals.append(np.nan)
        for nm in ['A', 'C']:
            res[nm][f'expr_{g}'] = pear(vals, [w[nm] for w in W])
except Exception as e:
    rep['expr_error'] = repr(e)[:150]

# Cox + Harrell C
try:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    import pandas as pd
    for nm in ['A', 'C']:
        rows = [{'T': w['os_time'], 'E': int(w['os_status']), 'x': w[nm],
                 'grade': w['grade'], 'mitotic': w['mitotic']}
                for w in W if w['os_time'] and w['os_status'] is not None]
        df = pd.DataFrame(rows)
        df['z'] = (df['x'] - df['x'].mean()) / df['x'].std(ddof=1)
        s1 = CoxPHFitter().fit(df[['T', 'E', 'z']], 'T', 'E').summary.loc['z']
        res[nm]['Cox_HR_perSD'] = {'n': int(len(df)), 'events': int(df['E'].sum()),
                                   'HR': round(float(np.exp(s1['coef'])), 4),
                                   'CI': [round(float(np.exp(s1['coef lower 95%'])), 4),
                                          round(float(np.exp(s1['coef upper 95%'])), 4)],
                                   'p': round(float(s1['p']), 5)}
        res[nm]['harrell_C'] = round(float(concordance_index(df['T'], -df['x'], df['E'])), 4)
except Exception as e:
    rep['cox_error'] = repr(e)[:250]

rep['A_TCGA_Score'] = res['A']
rep['C_feature_scores'] = res['C']

json.dump(rep, open(os.path.join(OUT, 'd15f_cardD_two_scores.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)

L = []
L.append('===== 卡 D 在两个评分口径下的对照（同一批 50 例）=====')
L.append(json.dumps({'n_common': rep['n_common']}, ensure_ascii=False))
L.append('')
L.append(f"{'端点':32} {'A (审计链用)':>26} {'C (手稿图用)':>26}")
def fmt(x):
    if not x: return '-'
    if 'r' in x: return f"r={x['r']}, p={x['p']} (n={x['n']})"
    if 'rho' in x: return f"rho={x['rho']}, p={x['p']}"
    return str(x)
keys = ['CN_raw', 'CN_ctrl_purity', 'CN_ctrl_purity_grade_mitotic', 'grade', 'mitotic',
        'expr_CD36', 'expr_CD8A', 'expr_PPARG']
for k in keys:
    L.append(f"{k:32} {fmt(rep['A_TCGA_Score'].get(k)):>26} {fmt(rep['C_feature_scores'].get(k)):>26}")
L.append('')
for nm, lab in [('A_TCGA_Score', 'A'), ('C_feature_scores', 'C')]:
    L.append(f"{lab}: Cox {json.dumps(rep[nm].get('Cox_HR_perSD'), ensure_ascii=False)}  HarrellC={rep[nm].get('harrell_C')}")
open(os.path.join(OUT, 'd15f_report.txt'), 'w', encoding='utf-8').write('\n'.join(L))
print('OK')
