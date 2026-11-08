# -*- coding: utf-8 -*-
"""
S33 · S32-a''（追加，不降维）：增殖混杂专项核查

动机：S32-a' 的「混杂源画像」实测出 score vs MKI67(增殖) r=+0.5707 (P=1.5e-05)，
      且 score vs paper_CN.cluster r=+0.5354。因此「评分 ↔ RICTOR/MTOR」必须再排除一层：
      **是否是通盘增殖/代谢活性驱动的假关联**。
做法：控制 MKI67（及 MKI67+纯度+CN 三重）后重算偏相关；并做「评分 vs RICTOR 在增殖分层内」检验。
输出：_work/s33_a_prime2.json / .txt
"""
import os, csv, json, math
import numpy as np
from scipy import stats
import pyreadr

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
lines = []
def L(s=''):
    print(s); lines.append(s)

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cli = {r['barcode']: r for r in
       csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}

def num(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except Exception:
        return None

ENS = {'RICTOR': 'ENSG00000133030', 'MTOR': 'ENSG00000198793', 'CD36': 'ENSG00000135218',
       'MKI67': 'ENSG00000148773', 'PCNA': 'ENSG00000132646', 'TOP2A': 'ENSG00000131747',
       'MCM2': 'ENSG00000117748', 'CCNB1': 'ENSG00000134057',
       'ADIPOQ': 'ENSG00000181092', 'GAPDH': 'ENSG00000111640'}
df = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
df = df[list(df.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in df.index}
cols = [str(c) for c in df.columns]

E = {}
for g, e in ENS.items():
    row = idx.get(e)
    if not row: continue
    E[g] = {c: float(df.loc[row, c]) for c in cols if c in score}

# 增殖指数 = 5 个增殖标志的 z-score 均值（更稳，抗单基因噪声）
prol_genes = [g for g in ['MKI67', 'PCNA', 'TOP2A', 'MCM2', 'CCNB1'] if g in E]
samples = sorted(set(score) & set.intersection(*[set(E[g]) for g in prol_genes]))
Z = np.array([[E[g][c] for g in prol_genes] for c in samples], float)
Zz = (Z - Z.mean(0)) / Z.std(0, ddof=1)
prol = {c: float(Zz[i].mean()) for i, c in enumerate(samples)}

L("S33 S32-a''  增殖混杂专项核查")
L('=' * 76)
L('增殖指数 = %s 的 z 均值（n=%d 样本）' % ('+'.join(prol_genes), len(samples)))
L('')

def partial(x, y, zs):
    X = np.asarray(x, float); Y = np.asarray(y, float)
    M = np.column_stack([np.ones(len(X))] + [np.asarray(z, float) for z in zs])
    k = M.shape[1] - 1
    if len(X) < k + 4: return None
    def rs(a):
        b, *_ = np.linalg.lstsq(M, a, rcond=None); return a - M @ b
    rx, ry = rs(X), rs(Y)
    if rx.std() == 0 or ry.std() == 0:
        return {'n': len(X), 'partial_r': None, 'p': None, 'note': '零变异'}
    r, _ = stats.pearsonr(rx, ry)
    d = len(X) - 2 - k
    if d < 1 or abs(r) >= 1: return None
    t = r * math.sqrt(d / (1 - r * r))
    return {'n': len(X), 'k_cov': k, 'partial_r': round(float(r), 4),
            'p': round(float(2 * stats.t.sf(abs(t), d)), 5)}

L('【核心】控制「增殖指数」后，评分 ↔ 轴基因是否仍相关')
L('%-9s %10s %18s %18s %24s' % ('gene', 'raw r', 'adj 增殖', 'adj 增殖+纯度', 'adj 增殖+纯度+CN'))
L('-' * 86)
def collect(g, extra=()):
    keys = [g] + list(extra)
    xs, ys, zs = [], [], {k: [] for k in extra}
    for c in samples:
        if g not in E or c not in E[g]: continue
        ok = True
        vals = {}
        for k in extra:
            v = prol.get(c) if k == 'PROL' else (num(cli.get(c, {}).get('paper_purity')) if k == 'PUR' else num(cli.get(c, {}).get('paper_CN.cluster')))
            if v is None: ok = False; break
            vals[k] = v
        if not ok: continue
        xs.append(E[g][c]); ys.append(score[c])
        for k in extra: zs[k].append(vals[k])
    return xs, ys, zs

res = {}
for g in ['RICTOR', 'MTOR', 'CD36']:
    x, y, z = collect(g, ('PROL',))
    r0, p0 = stats.pearsonr(x, y)
    a1 = partial(x, y, [z['PROL']])
    x2, y2, z2 = collect(g, ('PROL', 'PUR'))
    a2 = partial(x2, y2, [z2['PROL'], z2['PUR']])
    x3, y3, z3 = collect(g, ('PROL', 'PUR', 'CN')) if '_CNok' else (None, None, None)
    x3, y3, z3 = collect(g, ('PROL', 'PUR', 'CN'))
    a3 = partial(x3, y3, [z3['PROL'], z3['PUR'], z3['CN']]) if x3 else None
    def f(d): return '(不可估计)' if not d or d.get('partial_r') is None else 'r=%+.4f P=%.3g' % (d['partial_r'], d['p'])
    L('%-9s %10s %18s %18s %24s' % (g, 'r=%+.4f' % r0, f(a1), f(a2), f(a3)))
    res[g] = {'raw_r': round(float(r0), 4), 'adj_proliferation': a1,
              'adj_prol_purity': a2, 'adj_prol_purity_CN': a3}

L('')
L('【增殖分层内】把样本按增殖指数中位切两半，各层内重测「评分 ↔ RICTOR/MTOR」')
med = np.median([prol[c] for c in samples])
L('%-9s %-10s %5s %18s' % ('gene', 'stratum', 'n', 'r (P)'))
L('-' * 48)
strata = {}
for g in ['RICTOR', 'MTOR']:
    for tag, sel in [('prol_low', lambda c: prol[c] <= med), ('prol_high', lambda c: prol[c] > med)]:
        xs = [E[g][c] for c in samples if c in E.get(g, {}) and sel(c)]
        ys = [score[c] for c in samples if c in E.get(g, {}) and sel(c)]
        if len(xs) < 8:
            L('%-9s %-10s %5d  (过小)' % (g, tag, len(xs))); continue
        r, p = stats.pearsonr(xs, ys)
        strata.setdefault(g, {})[tag] = {'n': len(xs), 'r': round(float(r), 4), 'p': round(float(p), 5)}
        L('%-9s %-10s %5d  r=%+.4f  P=%.3g' % (g, tag, len(xs), r, p))
res['within_proliferation_stratum'] = strata
L('')

L('【对照】增殖指数自身与评分/轴基因（用于判读混杂方向）')
for lab, arr in [('score', [score[c] for c in samples]),
                 ('RICTOR', [E['RICTOR'][c] for c in samples]),
                 ('MTOR', [E['MTOR'][c] for c in samples]),
                 ('CD36', [E['CD36'][c] for c in samples]),
                 ('ADIPOQ', [E['ADIPOQ'][c] for c in samples if c in E['ADIPOQ']])]:
    pv = [prol[c] for c in samples] if lab != 'ADIPOQ' else [prol[c] for c in samples if c in E['ADIPOQ']]
    if len(arr) != len(pv): continue
    r, p = stats.pearsonr(arr, pv)
    L('  %-8s vs 增殖指数  r=%+.4f  P=%.3g' % (lab, r, p))
L('')

json.dump(res, open(os.path.join(OUT, 's33_a_prime2.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
open(os.path.join(OUT, 's33_a_prime2.txt'), 'w', encoding='utf-8').write('\n'.join(lines))
print('\n[写出] s33_a_prime2.json / s33_a_prime2.txt')
