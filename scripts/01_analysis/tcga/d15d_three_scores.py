# -*- coding: utf-8 -*-
"""三套「脂代谢评分」来源的相互核对
  A = data/raw/TCGA_Score.csv                     ['Score']            （审计链 d7/d11/d12 使用）
  B = data/results/02_feature_scores_37genes.csv  ['Lipid_Metabolism_Score'] (50 行)
  C = data/results/02_feature_scores.csv          ['Lipid_Metabolism_Score'] (253 行) ← 手稿全部图脚本使用
"""
import os, csv, json
import numpy as np
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")


def load(path, keycol, valcol):
    d = {}
    for r in csv.DictReader(open(path, encoding='utf-8-sig')):
        k = r.get(keycol) or r.get('')
        try:
            d[str(k).strip()] = float(r[valcol])
        except Exception:
            pass
    return d


A = load(os.path.join(RAW, 'TCGA_Score.csv'), 'Sample', 'Score')
B = load(os.path.join(RES, '02_feature_scores_37genes.csv'), 'Sample', 'Lipid_Metabolism_Score')
C = load(os.path.join(RES, '02_feature_scores.csv'), '', 'Lipid_Metabolism_Score')

rep = {'nA': len(A), 'nB': len(B), 'nC': len(C)}
rep['A_sample_rows'] = list(A.items())[:2]
rep['C_sample_rows'] = list(C.items())[:2]

common = sorted(set(A) & set(B) & set(C))
rep['n_common_ABC'] = len(common)
rep['n_common_AB'] = len(set(A) & set(B))


def cmp(x, y, nm):
    a = np.array([x[k] for k in common]); b = np.array([y[k] for k in common])
    if len(a) < 5: return
    pr, pp = stats.pearsonr(a, b); sr, sp = stats.spearmanr(a, b)
    rep[nm] = {'pearson': round(float(pr), 6), 'pearson_p': float(pp),
               'spearman': round(float(sr), 6), 'spearman_p': float(sp),
               'max_abs_diff': round(float(np.max(np.abs(a - b))), 6),
               'mean_x': round(float(a.mean()), 4), 'sd_x': round(float(a.std(ddof=1)), 4),
               'mean_y': round(float(b.mean()), 4), 'sd_y': round(float(b.std(ddof=1)), 4),
               'n_exactly_equal': int(np.sum(np.abs(a - b) < 1e-9))}


cmp(A, B, 'A_TCGA_Score__vs__B_37genes')
cmp(A, C, 'A_TCGA_Score__vs__C_feature_scores')
cmp(B, C, 'B_37genes__vs__C_feature_scores')

# 是否只是线性变换
for nm, (x, y) in {'A_vs_C': (A, C), 'A_vs_B': (A, B), 'B_vs_C': (B, C)}.items():
    a = np.array([x[k] for k in common]); b = np.array([y[k] for k in common])
    sl, ic = np.polyfit(a, b, 1)
    r2 = stats.pearsonr(a, b)[0] ** 2
    rep.setdefault('linear_fit', {})[nm] = {'slope': round(float(sl), 5),
                                            'intercept': round(float(ic), 5), 'r2': round(float(r2), 6)}

# C 的 253 行里有多少是非 TCGA / 非 DDLPS
rows = list(csv.DictReader(open(os.path.join(RES, '02_feature_scores.csv'), encoding='utf-8-sig')))
key = list(rows[0].keys())[0]
rep['C_index_col_name'] = key
rep['C_n_unique_index'] = len({r[key] for r in rows})
if 'Cohort' in rows[0]:
    from collections import Counter
    rep['C_cohort_counts'] = dict(Counter(r['Cohort'] for r in rows))
if 'Subtype' in rows[0]:
    from collections import Counter
    rep['C_subtype_counts'] = dict(Counter(r['Subtype'] for r in rows))

json.dump(rep, open(os.path.join(OUT, 'd15d_three_scores.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
