# -*- coding: utf-8 -*-
"""D7/D9 第三步：CN 复杂度信号稳健性 + 评分公式定位（为独立队列复现做准备）"""
import os, csv, json, math, re
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
cli = {r['barcode']: r for r in cl}

def num(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except Exception: return None

# --- 1. CN cluster vs 评分：控制纯度 ---
def partial(x, y, z):
    X = np.asarray(x, float); Y = np.asarray(y, float); Z = np.asarray(z, float)
    A = np.column_stack([np.ones_like(Z), Z])
    bx, *_ = np.linalg.lstsq(A, X, rcond=None); by, *_ = np.linalg.lstsq(A, Y, rcond=None)
    rx, ry = X - A @ bx, Y - A @ by
    n = len(X); r, _ = stats.pearsonr(rx, ry); d = n - 3
    t = r * math.sqrt(d / (1 - r * r)); p = float(2 * stats.t.sf(abs(t), d))
    return {'n': n, 'partial_r': round(float(r), 4), 'p': round(p, 5)}

rec = []
for c, s in score.items():
    r0 = cli.get(c)
    if not r0: continue
    rec.append({'s': s, 'cn': num(r0.get('paper_CN.cluster')), 'pu': num(r0.get('paper_purity')),
                'pl': num(r0.get('paper_ploidy'))})
ok = [r for r in rec if r['cn'] is not None and r['pu'] is not None]
rep['score_vs_CNcluster_adj_purity'] = partial([r['cn'] for r in ok], [r['s'] for r in ok], [r['pu'] for r in ok])
rep['CNcluster_vs_purity'] = None
if len(ok) >= 8:
    r, p = stats.pearsonr([r['cn'] for r in ok], [r['pu'] for r in ok])
    rep['CNcluster_vs_purity'] = {'r': round(float(r), 4), 'p': round(float(p), 5), 'n': len(ok)}
# 按 CN cluster 分组的评分
grp = {}
for r in rec:
    if r['cn'] is None: continue
    grp.setdefault(int(r['cn']), []).append(r['s'])
rep['score_by_CNcluster'] = {k: {'n': len(v), 'mean': round(float(np.mean(v)), 4), 'sd': round(float(np.std(v, ddof=1)) if len(v) > 1 else 0, 4)} for k, v in sorted(grp.items())}
ks = [v for v in grp.values() if len(v) >= 2]
if len(ks) >= 2:
    H, pk = stats.kruskal(*ks)
    rep['kruskal_CNcluster'] = {'H': round(float(H), 4), 'p': round(float(pk), 6)}

# --- 2. CN cluster vs MDM2 表达 ---
df = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
df = df[list(df.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in df.index}
ENS = {'MDM2': 'ENSG00000135679', 'FRS2': 'ENSG00000166225', 'CDK4': 'ENSG00000135446', 'HMGA2': 'ENSG00000149948'}
t = {}
for g, e in ENS.items():
    row = idx.get(e); x, y = [], []
    for c in df.columns:
        c = str(c)
        if c not in score: continue
        r0 = cli.get(c); cn = num(r0.get('paper_CN.cluster')) if r0 else None
        if cn is None: continue
        try: v = float(df.loc[row, c])
        except Exception: continue
        x.append(cn); y.append(v)
    if len(x) >= 8:
        r, p = stats.pearsonr(x, y)
        t[g] = {'r_vs_CNcluster': round(float(r), 4), 'p': round(float(p), 5), 'n': len(x)}
rep['CNcluster_vs_amplicon_expression'] = t

# --- 3. 定位脂代谢评分的计算脚本（SREBF1+SCD+HMGCR+ACLY 同现）---
cands = []
root0 = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1")
for root, dirs, files in os.walk(root0):
    dirs[:] = [d for d in dirs if d not in ('.Rproj.user',)]
    for f in files:
        if not f.lower().endswith(('.r', '.py', '.rmd')): continue
        p = os.path.join(root, f)
        try: txt = open(p, encoding='utf-8', errors='ignore').read()
        except Exception: continue
        n = sum(txt.count(g) for g in ('SREBF1', 'SCD', 'HMGCR', 'ACLY'))
        if n >= 3 and re.search(r'[Ss]core', txt):
            hits = [l.strip()[:200] for l in txt.splitlines()
                    if re.search(r'[Ss]core\s*(<-|=)', l) and not l.strip().startswith('#')]
            cands.append({'file': os.path.relpath(p, root0), 'nhits': n, 'scorelines': hits[:6]})
rep['score_formula_candidates'] = sorted(cands, key=lambda d: -d['nhits'])[:8]

json.dump(rep, open(os.path.join(OUT, 'd7_cncluster.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
