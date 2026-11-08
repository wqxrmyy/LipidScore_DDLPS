# -*- coding: utf-8 -*-
"""反推并验证评分公式：以 TCGA_6genes_expr.rds 复现 TCGA_Score.csv"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}

res = pyreadr.read_r(os.path.join(RAW, 'TCGA_6genes_expr.rds'))
g6 = res[list(res.keys())[0]]
genes = [str(g) for g in g6.index]
cols = [str(c) for c in g6.columns]
mat = g6.values.astype(float)
rep['genes'] = genes
rep['shape'] = list(g6.shape)

def zs(m):
    return (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, ddof=1, keepdims=True)

z = zs(mat)
cands = {
    'mean_z_all6': z.mean(axis=0),
    'mean_z_no_SREBF1': z[[i for i, g in enumerate(genes) if g != 'SREBF1']].mean(axis=0),
    'mean_z_no_PPARG': z[[i for i, g in enumerate(genes) if g != 'PPARG']].mean(axis=0),
    'mean_raw_all6': mat.mean(axis=0),
    'sum_z_all6': z.sum(axis=0),
}
y = np.array([score.get(c, np.nan) for c in cols], float)
ok = ~np.isnan(y)
out = {}
for k, v in cands.items():
    r, p = stats.pearsonr(v[ok], y[ok])
    out[k] = {'pearson_vs_TCGA_Score': round(float(r), 5), 'max_abs_diff_from_unit': round(float(np.max(np.abs(v[ok] - y[ok]))), 4)}
rep['candidate_formulas'] = out

# 最佳候选的尺度关系
best = cands['mean_z_all6']
ratio = y[ok] / best[ok]
rep['ratio_Score_over_meanZ'] = {'mean': round(float(np.mean(ratio)), 5), 'sd': round(float(np.std(ratio, ddof=1)), 5),
                                 'min': round(float(np.min(ratio)), 5), 'max': round(float(np.max(ratio)), 5)}
# 线性回归 y = a + b*z
b, a = np.polyfit(best[ok], y[ok], 1)
rep['linear_fit'] = {'slope': round(float(b), 5), 'intercept': round(float(a), 5),
                     'r2': round(float(stats.pearsonr(best[ok], y[ok])[0] ** 2), 6)}

json.dump(rep, open(os.path.join(OUT, 'score_formula.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
