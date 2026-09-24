# -*- coding: utf-8 -*-
"""GSE221492（bulk LPS atlas, 53 样本/41 患者）独立队列验证 —— **规范实现**。

本脚本是 scripts/07_bulk_validation/ 下**两套并存实现**之一：
  · 本文件（Python）  <- 本文正文所报数字的**原始来源**（前身 d11_gse221492.py）
  · 01_lipid_score_bulk.R（R） <- 独立重实现，用于交叉校验
两者必须给出逐位一致的结果，由 scripts/09_qc/assert_manuscript_numbers.py 断言。

⚠️ 两处易错口径（R 侧曾在此出错，Python 侧从一开始就是对的）：
  1) **library size 必须取自完整矩阵**（合并重复 hgnc_symbol 之前）。
     合并会丢弃 12,873 行（28.3%）的计数，使分母偏小 0.47%。
  2) **DDLPS_sWD 不得并入 DDLPS_WD**。严格配对必须是「同患者同时具备真正的
     DDLPS_WD 与 DDLPS_DD 样本」-> n=8；把 sWD 当 WD 会得到本文从未报告的 n=10。

覆盖内容
--------
① 评分公式（6 基因 z 分均值，逐队列计算）
② B-P5 自检：评分对 WD/DD 的判别 AUC（配对 + 非配对）
③ A′ 复现：评分 vs 12q 扩增子成员（MDM2/FRS2/CDK4/HMGA2）表达剂量
④ 评分与 CD36 / 免疫标志 CD8A 的关联
⑤ 分量基因各自的 DD-vs-WD 效应（评分是否内部一致）
"""
import os, gzip, csv, json, math
import numpy as np
from scipy import stats

ROOT = os.environ.get('LPS_ROOT')
if not ROOT:
    ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
print('项目根:', ROOT)
OUT_OBJ = os.path.join(ROOT, 'results', 'objects')
OUT_TAB = os.path.join(ROOT, 'results', 'tables')
os.makedirs(OUT_OBJ, exist_ok=True)
os.makedirs(OUT_TAB, exist_ok=True)
GZ = os.path.join(ROOT, 'data', 'GSE221492_bulk', 'raw',
                  'GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz')
GENES = ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN']
EXTRA = ['MDM2', 'FRS2', 'CDK4', 'HMGA2', 'CD36', 'CD8A']
WANT = GENES + EXTRA
rep = {}

# ---------- 读矩阵 ----------
samples, mat = None, {}
with gzip.open(GZ, 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples = hdr[2:]
    keep = set(WANT)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 3: continue
        sym = p[1].strip()
        if sym in keep and sym not in mat:
            mat[sym] = np.array([float(x) if x not in ('', 'NA') else np.nan for x in p[2:]])
rep['n_samples'] = len(samples)
rep['genes_found'] = sorted(mat.keys())
rep['genes_missing'] = [g for g in WANT if g not in mat]

# ---------- 标签解析 ----------
def parse(s):
    b = s.split('_')
    pid = b[0].lstrip('Pp')
    lab = '_'.join(b[1:])
    return pid, lab

meta = [{'sample': s, 'patient': parse(s)[0], 'label': parse(s)[1]} for s in samples]
rep['label_counts'] = {k: sum(1 for m in meta if m['label'] == k) for k in sorted(set(m['label'] for m in meta))}

# ---------- log2CPM ----------
# library size 必须取自**完整矩阵**（合并重复基因符号之前）—— 见文件头说明 1)
colsum = np.zeros(len(samples))
with gzip.open(GZ, 'rt', encoding='utf-8', errors='ignore') as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 3: continue
        try:
            colsum += np.array([float(x) if x not in ('', 'NA') else 0.0 for x in p[2:]])
        except Exception:
            pass
l2cpm = {}
for g in mat:
    v = mat[g] / np.where(colsum > 0, colsum, np.nan) * 1e6
    l2cpm[g] = np.log2(v + 1)


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)


def score_from(genes, mask):
    zs = np.vstack([z(l2cpm[g][mask]) for g in genes if g in l2cpm])
    return np.nanmean(zs, axis=0)


def auc(x1, x0):
    """AUC + Mann-Whitney U + 精确/正态 p"""
    x1 = np.asarray(x1, float); x0 = np.asarray(x0, float)
    x1 = x1[~np.isnan(x1)]; x0 = x0[~np.isnan(x0)]
    n1, n0 = len(x1), len(x0)
    U = 0.0
    for a in x1:
        U += (x0 < a).sum() + 0.5 * (x0 == a).sum()
    a = U / (n1 * n0)
    try:
        _, p = stats.mannwhitneyu(x1, x0, alternative='two-sided')
    except Exception:
        p = float('nan')
    return {'auc': round(float(a), 4), 'U': round(float(U), 2), 'n1': n1, 'n0': n0,
            'p': round(float(p), 6)}


def comp(x1, x0):
    x1 = np.asarray(x1, float); x0 = np.asarray(x0, float)
    x1 = x1[~np.isnan(x1)]; x0 = x0[~np.isnan(x0)]
    t, p = stats.ttest_ind(x1, x0, equal_var=False)
    return {'mean1': round(float(np.mean(x1)), 4), 'mean0': round(float(np.mean(x0)), 4),
            'delta': round(float(np.mean(x1) - np.mean(x0)), 4), 'welch_t': round(float(t), 4),
            'p': round(float(p), 6)}


idx = {m['sample']: i for i, m in enumerate(meta)}
lab = np.array([m['label'] for m in meta])

# ================= ② B-P5 自检 =================
b5 = {}
# 情形 1：z 在全部 53 样本内计算
sc_all = score_from(GENES, np.ones(len(samples), bool))
rep['score_all53'] = {m['sample']: round(float(sc_all[i]), 4) for i, m in enumerate(meta)}

# 严格 DDLPS 成分对比：DDLPS_DD vs DDLPS_WD（z 在全部样本内）
mA = np.isin(lab, ['DDLPS_DD']); mB = np.isin(lab, ['DDLPS_WD'])
b5['DD_vs_WD_DDLPS_specimens'] = {'auc': auc(sc_all[mA], sc_all[mB]), 'compare': comp(sc_all[mA], sc_all[mB])}

# 广义：DD vs 所有 WD 型（DDLPS_WD + DDLPS_sWD + WDLPS）
mW = np.isin(lab, ['DDLPS_WD', 'DDLPS_sWD', 'WDLPS'])
b5['DD_vs_anyWD'] = {'auc': auc(sc_all[mA], sc_all[mW]), 'compare': comp(sc_all[mA], sc_all[mW])}

# 对照：DD vs 脂肪瘤/正常脂肪
for ctrl in ['lipoma', 'peritumoral_fat', 'Hibernoma']:
    mc = np.isin(lab, [ctrl])
    if mc.sum() >= 2:
        b5[f'DD_vs_{ctrl}'] = {'auc': auc(sc_all[mA], sc_all[mc]), 'compare': comp(sc_all[mA], sc_all[mc])}

# 情形 2：z 只在 DDLPS 标本内计算（更贴近 TCGA 口径：TCGA 也是 DDLPS-only 队列）
mD = np.isin(lab, ['DDLPS_DD', 'DDLPS_WD', 'DDLPS_sWD'])
sc_ddlps = np.full(len(samples), np.nan)
sc_ddlps[mD] = score_from(GENES, mD)
b5['z_within_DDLPS_only'] = {
    'DD_vs_WD': {'auc': auc(sc_ddlps[mA], sc_ddlps[mB]), 'compare': comp(sc_ddlps[mA], sc_ddlps[mB])},
    'DD_vs_sWD': ({'auc': auc(sc_ddlps[mA], sc_ddlps[np.isin(lab, ['DDLPS_sWD'])]),
                   'compare': comp(sc_ddlps[mA], sc_ddlps[np.isin(lab, ['DDLPS_sWD'])])}
                  if np.isin(lab, ['DDLPS_sWD']).sum() >= 2 else None),
}

# 配对分析（8 对严格 WD/DD 同患者）
pairs = {}
for m in meta:
    pairs.setdefault(m['patient'], {})[m['label']] = m['sample']
strict = [(p, d['DDLPS_WD'], d['DDLPS_DD']) for p, d in sorted(pairs.items())
          if 'DDLPS_WD' in d and 'DDLPS_DD' in d]
pa = []
for pid, w, dd in strict:
    a1 = sc_all[idx[dd]]; a0 = sc_all[idx[w]]
    b1 = sc_ddlps[idx[dd]]; b0 = sc_ddlps[idx[w]]
    pa.append({'patient': pid, 'WD': round(float(a0), 4), 'DD': round(float(a1), 4),
               'delta_all53': round(float(a1 - a0), 4), 'delta_withinDDLPS': round(float(b1 - b0), 4)})
b5['n_strict_pairs'] = len(strict)
b5['paired_detail'] = pa
if len(pa) >= 3:
    dA = [p['delta_all53'] for p in pa]
    w_, wp = stats.wilcoxon(dA)
    t_, tp = stats.ttest_rel([p['DD'] for p in pa], [p['WD'] for p in pa])
    b5['paired_test_all53'] = {'n_pairs': len(pa), 'mean_delta': round(float(np.mean(dA)), 4),
                               'median_delta': round(float(np.median(dA)), 4),
                               'n_pos': sum(1 for x in dA if x > 0), 'n_neg': sum(1 for x in dA if x < 0),
                               'wilcoxon_W': round(float(w_), 3), 'wilcoxon_p': round(float(wp), 6),
                               'paired_t': round(float(t_), 4), 'paired_t_p': round(float(tp), 6)}
    dW = [p['delta_withinDDLPS'] for p in pa]
    w2, wp2 = stats.wilcoxon(dW)
    b5['paired_test_withinDDLPS'] = {'mean_delta': round(float(np.mean(dW)), 4),
                                     'n_pos': sum(1 for x in dW if x > 0), 'n_neg': sum(1 for x in dW if x < 0),
                                     'wilcoxon_p': round(float(wp2), 6)}
rep['B5_selfcheck'] = b5

# ================= ③ A′ 复现（12q 剂量）=================
def cor(a, b, method='spearman'):
    m = ~(np.isnan(a) | np.isnan(b))
    a, b = a[m], b[m]
    if len(a) < 5: return None
    if method == 'spearman':
        r, p = stats.spearmanr(a, b)
    else:
        r, p = stats.pearsonr(a, b)
    return {'r': round(float(r), 4), 'p': round(float(p), 6), 'n': len(a)}


repr_ = {}
# 全 53 样本内
for g in EXTRA:
    repr_[g] = {'all53': {'spearman': cor(sc_all, l2cpm[g]), 'pearson': cor(sc_all, l2cpm[g], 'pearson')},
                'DDLPS_only': {'spearman': cor(sc_ddlps, l2cpm[g]),
                               'pearson': cor(sc_ddlps, l2cpm[g], 'pearson')}}
rep['Aprime_replication'] = repr_

# ================= ⑤ 分量基因各自的 DD-vs-WD 效应 =================
comp_genes = {}
for g in GENES + ['CD36', 'CD8A']:
    if g not in l2cpm: continue
    v = l2cpm[g]
    dd = v[np.isin(lab, ['DDLPS_DD'])]; wd = v[np.isin(lab, ['DDLPS_WD'])]
    dd = dd[~np.isnan(dd)]; wd = wd[~np.isnan(wd)]
    lfc = float(np.mean(dd) - np.mean(wd))
    t, p = stats.ttest_ind(dd, wd, equal_var=False)
    comp_genes[g] = {'mean_DD': round(float(np.mean(dd)), 4), 'mean_WD': round(float(np.mean(wd)), 4),
                     'log2FC_DD_minus_WD': round(lfc, 4), 'welch_p': round(float(p), 6),
                     'auc_DD_vs_WD': auc(dd, wd)['auc']}
rep['component_gene_DD_vs_WD'] = comp_genes

# 评分分量内部一致性（在 DDLPS 标本内）
cons = {}
for i in range(len(GENES)):
    for j in range(i + 1, len(GENES)):
        gi, gj = GENES[i], GENES[j]
        if gi not in l2cpm or gj not in l2cpm: continue
        c = cor(l2cpm[gi][mD], l2cpm[gj][mD])
        if c: cons[f'{gi}~{gj}'] = c
rep['component_internal_corr_DDLPS'] = cons

json.dump(rep, open(os.path.join(OUT_OBJ, 'd11_gse221492.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)

# ---- 派生 CSV：供交叉校验脚本与人工核对 ----
import csv as _csv
_pa = b5['paired_detail']
with open(os.path.join(OUT_TAB, 'bulk_lipid_score_py_paired_delta.csv'), 'w',
          newline='', encoding='utf-8') as _f:
    _w = _csv.writer(_f)
    _w.writerow(['patient', 'DDLPS_DD', 'DDLPS_WD', 'delta'])
    for _p in sorted(_pa, key=lambda z: int(z['patient'])):
        _w.writerow(['P' + str(_p['patient']), _p['DD'], _p['WD'], _p['delta_all53']])
_pt = b5['paired_test_all53']
with open(os.path.join(OUT_TAB, 'bulk_lipid_score_py_summary.csv'), 'w',
          newline='', encoding='utf-8') as _f:
    _w = _csv.writer(_f)
    _w.writerow(['analysis', 'n_pairs', 'V_statistic', 'p_value', 'mean_delta',
                 'median_delta', 'n_dd_lt_wd'])
    _w.writerow(['Main (strict pairs)', _pt['n_pairs'], _pt['wilcoxon_W'],
                 _pt['wilcoxon_p'], _pt['mean_delta'], _pt['median_delta'], _pt['n_neg']])

# ---- 本文报告值断言：任何一位改变都必须失败，不得静默通过 ----
_c = b5['DD_vs_WD_DDLPS_specimens']['compare']
_a = b5['DD_vs_WD_DDLPS_specimens']['auc']
assert b5['n_strict_pairs'] == 8, b5['n_strict_pairs']
assert _pt['n_neg'] == 6, _pt['n_neg']
assert _pt['wilcoxon_p'] == 0.078125, _pt['wilcoxon_p']
assert round(_pt['mean_delta'], 3) == -0.587, _pt['mean_delta']
assert round(_c['mean1'], 3) == -0.305 and round(_c['mean0'], 3) == 0.164, _c
assert round(_c['delta'], 3) == -0.469 and round(_c['p'], 3) == 0.025, _c
assert round(_a['auc'], 3) == 0.249 and round(_a['p'], 3) == 0.013, _a
print('\n✓ 本文报告值断言通过（n=8 / 6-of-8 / P=0.078125 / mean Δ=-0.587 / '
      '非配对 -0.305,+0.164,Δ=-0.469,P=0.025 / AUC=0.249,P=0.013）')
print(json.dumps({k: v for k, v in rep.items()
                  if k not in ('score_all53', 'paired_detail')}, ensure_ascii=False, indent=2))
