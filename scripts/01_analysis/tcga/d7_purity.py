# -*- coding: utf-8 -*-
"""D7/D9 支撑分析：控制肿瘤纯度后的偏相关；以及评分定义的定位。"""
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

df = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
df = df[list(df.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in df.index}
ENS = {'MDM2': 'ENSG00000135679', 'CD8A': 'ENSG00000198851', 'FRS2': 'ENSG00000166225',
       'CD36': 'ENSG00000135218', 'PPARG': 'ENSG00000132170', 'SCD': 'ENSG00000099194',
       'FASN': 'ENSG00000169710', 'HMGCR': 'ENSG00000113161', 'ACLY': 'ENSG00000131473',
       'ADIPOQ': 'ENSG00000181092', 'PLIN1': 'ENSG00000166819', 'CPT1A': 'ENSG00000110090'}

cols = [str(c) for c in df.columns]
data = {}
for g, e in ENS.items():
    row = idx.get(e)
    if not row: continue
    vals = {}
    for c in cols:
        if c in score:
            try: vals[c] = float(df.loc[row, c])
            except Exception: pass
    data[g] = vals

def partial_corr(x, y, z):
    """控制 z 后 x,y 的偏相关 + n"""
    X = np.asarray(x, float); Y = np.asarray(y, float); Z = np.asarray(z, float)
    if len(X) < 6: return None
    # 残差化
    def resid(a, b):
        A = np.column_stack([np.ones_like(b), b])
        beta, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ beta
    rx, ry = resid(X, Z), resid(Y, Z)
    n = len(X); k = 1
    r, _ = stats.pearsonr(rx, ry)
    dfree = n - 2 - k
    if dfree < 1 or abs(r) >= 1: return None
    t = r * math.sqrt(dfree / (1 - r * r))
    p = float(2 * stats.t.sf(abs(t), dfree))
    return {'n': n, 'partial_r': round(float(r), 4), 'p': round(p, 5)}

pairs = [('MDM2', None), ('CD8A', None), ('FRS2', None), ('CD36', None), ('PPARG', None),
         ('SCD', None), ('ADIPOQ', None), ('CPT1A', None)]
tab = {}
for g, _ in pairs:
    if g not in data: continue
    x, y, z, cns = [], [], [], []
    for c, s in score.items():
        if c not in data[g]: continue
        r0 = cli.get(c)
        if not r0: continue
        pu = num(r0.get('paper_purity'))
        cn = num(r0.get('paper_CN.cluster'))
        x.append(data[g][c]); y.append(s); z.append(pu); cns.append(cn)
    m = [i for i, v in enumerate(z) if v is not None]
    if len(m) >= 8:
        xa = [x[i] for i in m]; ya = [y[i] for i in m]; za = [z[i] for i in m]
        raw_r, raw_p = stats.pearsonr(xa, ya)
        tab[g] = {
            'raw_r': round(float(raw_r), 4), 'raw_p': round(float(raw_p), 5), 'n_purity_ok': len(m),
            'partial_adj_purity': partial_corr(xa, ya, za),
        }
        # 也控制 CN 簇
        m2 = [i for i, v in enumerate(cns) if v is not None]
        if len(m2) >= 8:
            xb = [x[i] for i in m2]; yb = [y[i] for i in m2]; zb = [cns[i] for i in m2]
            tab[g]['partial_adj_CNcluster'] = partial_corr(xb, yb, zb)
rep['purity_adjusted'] = tab

# 评分本身 vs 纯度 / CN簇（已在 d7_mdm2_feasibility 有，此处重列以对照）
m = [(score[c], num(cli[c]['paper_purity'])) for c in score if c in cli and num(cli[c].get('paper_purity')) is not None]
if len(m) >= 8:
    xs = [a for a, b in m]; ys = [b for a, b in m]
    r, p = stats.pearsonr(ys, xs)
    rep['score_vs_purity'] = {'r': round(float(r), 4), 'p': round(float(p), 5), 'n': len(m)}
m = [(score[c], num(cli[c]['paper_CN.cluster'])) for c in score if c in cli and num(cli[c].get('paper_CN.cluster')) is not None]
if len(m) >= 8:
    xs = [a for a, b in m]; ys = [b for a, b in m]
    r, p = stats.pearsonr(ys, xs); rs, ps = stats.spearmanr(ys, xs)
    rep['score_vs_CNcluster'] = {'pearson': [round(float(r), 4), round(float(p), 5)],
                                 'spearman': [round(float(rs), 4), round(float(ps), 5)], 'n': len(m)}

# ---- 定位评分定义脚本 ----
defs = []
for root, dirs, files in os.walk(os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1")):
    dirs[:] = [d for d in dirs if d not in ('.Rproj.user',)]
    for f in files:
        if not f.lower().endswith(('.r', '.py', '.rmd')): continue
        p = os.path.join(root, f)
        try: t = open(p, encoding='utf-8', errors='ignore').read()
        except Exception: continue
        for i, l in enumerate(t.splitlines(), 1):
            if re.search(r'(Score\s*<-\s*|score\s*=\s*|Lipid_?Score\s*<-)', l) and not re.search(r'^\s*#', l):
                defs.append({'file': os.path.relpath(p, os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1")), 'line': i, 'text': l.strip()[:240]})
            elif re.search(r'zscore|scale\(|rowMeans|colMeans', l) and re.search(r'[Ss]core', l) and not re.search(r'^\s*#', l):
                defs.append({'file': os.path.relpath(p, os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1")), 'line': i, 'text': l.strip()[:240]})
rep['score_definition_hits'] = defs[:40]

json.dump(rep, open(os.path.join(OUT, 'd7_purity.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
