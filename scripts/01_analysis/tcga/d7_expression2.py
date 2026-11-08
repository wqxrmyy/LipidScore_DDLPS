# -*- coding: utf-8 -*-
"""D7 第二步（修正版）：把 MDM2 及相关基因按 Ensembl ID 定位，检验与评分的关联。"""
import os, json, csv, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")

ENS = {
    'MDM2': 'ENSG00000135679', 'CDK4': 'ENSG00000135446', 'HMGA2': 'ENSG00000149948',
    'DDIT3': 'ENSG00000175197', 'FRS2': 'ENSG00000166225', 'PPARG': 'ENSG00000132170',
    'CD36': 'ENSG00000135218', 'FASN': 'ENSG00000169710', 'SCD': 'ENSG00000099194',
    'SREBF1': 'ENSG00000072310', 'HMGCR': 'ENSG00000113161', 'ACLY': 'ENSG00000131473',
    'PLIN1': 'ENSG00000166819', 'ADIPOQ': 'ENSG00000181092', 'LPL': 'ENSG00000175445',
    'FABP4': 'ENSG00000170323', 'CPT1A': 'ENSG00000110090', 'CD8A': 'ENSG00000198851',
    'PDCD1': 'ENSG00000188389', 'LIPE': 'ENSG00000079435',
}

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}

rep = {'n_score': len(score)}

def cortab(x, y):
    if len(x) < 6: return None
    r, p = stats.pearsonr(x, y); rs, ps = stats.spearmanr(x, y)
    return {'n': len(x), 'pearson': [round(float(r), 4), round(float(p), 5)],
            'spearman': [round(float(rs), 4), round(float(ps), 5)],
            'x_range': [round(float(min(x)), 3), round(float(max(x)), 3)],
            'x_sd': round(float(np.std(x, ddof=1)), 4)}

for fn in ['TCGA_LPS_real_log2cpm.rds', 'TCGA_6genes_expr.rds']:
    res = pyreadr.read_r(os.path.join(RAW, fn))
    df = res[list(res.keys())[0]]
    idx = [str(i) for i in df.index]
    cols = [str(c) for c in df.columns]
    out = {'shape': list(df.shape), 'n_col_samples_matched': sum(1 for c in cols if c in score)}
    # 基因定位
    loc = {}
    for g, e in ENS.items():
        hit = None
        for i in idx:
            if i == g:
                hit = i; break
        if hit is None:
            for i in idx:
                if i.split('.')[0] == e:
                    hit = i; break
        if hit: loc[g] = hit
    out['genes_found'] = sorted(loc.keys())
    out['genes_missing'] = sorted(set(ENS) - set(loc))
    tab = {}
    for g, i in loc.items():
        x, y = [], []
        for c in cols:
            s = score.get(c)
            if s is None: continue
            try: v = float(df.loc[i, c])
            except Exception: continue
            if math.isnan(v): continue
            x.append(v); y.append(s)
        tab[g] = cortab(x, y)
    out['score_vs_gene'] = tab
    rep[fn] = out

json.dump(rep, open(os.path.join(OUT, 'd7_expression.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
