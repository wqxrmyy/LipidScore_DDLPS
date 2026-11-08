# -*- coding: utf-8 -*-
"""判定 C（02_feature_scores.csv 的 Lipid_Metabolism_Score，手稿全部图脚本使用）
究竟是：(i) 同一 6 基因复合评分在「253 样本全域」重新标准化后的版本；
        (ii) 还是一个完全不同的构造。
做法：用 TCGA_LPS_real_log2cpm.rds 在全域（feature_scores 的 253 个样本）内对 6 基因做 z 分、
     取均值，与 C 逐样本比对；同时用 6gene 矩阵（50 样本）做阳性对照。
"""
import os, csv, json
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}
ENS = {'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
       'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710'}


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
C = load(os.path.join(RES, '02_feature_scores.csv'), '', 'Lipid_Metabolism_Score')

big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
cols = [str(c) for c in big.columns]
bidx = {str(k).split('.')[0]: str(k) for k in big.index}
rep['log2cpm_ncol'] = len(cols)
rep['n_C'] = len(C)
rep['n_A'] = len(A)
rep['n_cols_overlap_C'] = len(set(cols) & set(C))

E = {}
for g, e in ENS.items():
    row = bidx.get(e)
    if not row: continue
    v = {}
    for c in cols:
        try: v[c] = float(big.loc[row, c])
        except Exception: pass
    E[g] = v
rep['genes_found_in_log2cpm'] = sorted(E.keys())


def zmean(genes, ss):
    zs = []
    for g in genes:
        v = np.array([E[g].get(s, np.nan) for s in ss], float)
        zs.append((v - np.nanmean(v)) / np.nanstd(v, ddof=1))
    return np.nanmean(np.vstack(zs), axis=0)


def cmpv(x, y, label):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    if len(x) < 8: return None
    pr, pp = stats.pearsonr(x, y); sr, sp = stats.spearmanr(x, y)
    d = {'n': len(x), 'pearson': round(float(pr), 6), 'pearson_p': float(pp),
         'spearman': round(float(sr), 6), 'spearman_p': float(sp),
         'max_abs_diff': round(float(np.max(np.abs(x - y))), 6)}
    rep[label] = d
    return d


G6 = list(ENS.keys())

# ① 全域 253 样本内标准化
ss253 = [s for s in cols if s in C]
rep['n_used_253'] = len(ss253)
v253 = zmean(G6, ss253)
cmpv(v253, [C[s] for s in ss253], 'C_vs_6gene_zmean_over_all_samples_used')

# ② 只在 C 的 DDLPS 子集内标准化（50）
rows = list(csv.DictReader(open(os.path.join(RES, '02_feature_scores.csv'), encoding='utf-8-sig')))
key = list(rows[0].keys())[0]
sub = {}
if 'Subtype' in rows[0]:
    for r in rows:
        sub[str(r[key]).strip()] = r.get('Subtype')
    from collections import Counter
    rep['C_subtype_counts'] = dict(Counter(sub.values()))
dd50 = [s for s in ss253 if sub.get(s) == 'DDLPS']
rep['n_ddlps_in_C'] = len(dd50)
if len(dd50) >= 20:
    v50 = zmean(G6, dd50)
    cmpv(v50, [C[s] for s in dd50], 'C_vs_6gene_zmean_within_DDLPS_subset')
    # 与 A 对照
    a50 = [A[s] for s in dd50 if s in A]
    if len(a50) == len(dd50):
        cmpv(a50, [C[s] for s in dd50], 'A_vs_C_within_DDLPS')
        cmpv(a50, v50, 'A_vs_6gene_zmean_from_log2cpm_within_DDLPS')

# ③ 6gene 矩阵（阳性对照）：应与 A 精确一致
res6 = pyreadr.read_r(os.path.join(RAW, 'TCGA_6genes_expr.rds'))
g6 = res6[list(res6.keys())[0]]
g6n = [str(g) for g in g6.index]
c6 = [str(c) for c in g6.columns]
M = g6.values.astype(float)
v6 = np.nanmean(np.vstack([(M[i] - np.nanmean(M[i])) / np.nanstd(M[i], ddof=1) for i in range(M.shape[0])]), axis=0)
cmpv(v6, [A[c] for c in c6], 'A_vs_6gene_matrix_zmean_positive_control')

# ④ C 与 A、B 在 DDLPS 50 例上的方向一致性（符号翻转检验）
if len(dd50) >= 20 and len([s for s in dd50 if s in A]) == len(dd50):
    a = np.array([A[s] for s in dd50]); c = np.array([C[s] for s in dd50])
    rep['A_vs_C_sign_agreement'] = round(float(np.mean(np.sign(a - a.mean()) == np.sign(c - c.mean()))), 4)

json.dump(rep, open(os.path.join(OUT, 'd15e_C_identity.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
