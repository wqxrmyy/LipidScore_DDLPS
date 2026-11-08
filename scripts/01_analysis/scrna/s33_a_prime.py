# -*- coding: utf-8 -*-
"""
S33 · S32-a'  「评分 A ↔ mTORC2 轴」混杂核查（不降维版）

目的：把主交付件 §10 实测到的 r(RICTOR)=+0.5493 / r(MTOR)=+0.5585 做成**稳健**结果，
      并主动排除三类混杂：肿瘤纯度、CN 簇、脂肪含量（adipogenic 对照）。
      全部为加固步骤（非降级），任一稳健性检验不通过都如实登记而不隐藏。

口径（全部沿用项目既有实现，不另立方）：
  · 评分 = TCGA_Score.csv 唯一（D18 锁 A）
  · 纯度 = TCGA_clinical_real.csv 的 paper_purity（论文实测值，非 ESTIMATE 推测）
  · CN  = paper_CN.cluster / paper_ploidy
  · 偏相关 = 残差化 OLS（沿用 d7_purity.py 的 partial_corr 写法）
  · 基因由 ENSG 定位（s32_strata_probe.R 已核验的 ID），不凭印象硬编码
输出：_work/s33_a_prime.json / _work/s33_a_prime.txt
"""
import os, csv, json, math
import numpy as np
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}
lines = []
def L(s=''):
    print(s); lines.append(s)

# ---------------------------------------------------------------- 输入
score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
cli = {r['barcode']: r for r in cl}

def num(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except Exception:
        return None

# 基因：ENSG 已由 s32_strata_probe.R 核验；ADIPOQ/PLIN1/CPT1A 来自 d7_purity.py
ENS = {
    'RICTOR': 'ENSG00000133030', 'MTOR': 'ENSG00000198793', 'CD36': 'ENSG00000135218',
    'PPARG': 'ENSG00000132170', 'FABP4': 'ENSG00000170323',
    'ADIPOQ': 'ENSG00000181092', 'PLIN1': 'ENSG00000166819', 'CPT1A': 'ENSG00000110090',
    'MDM2': 'ENSG00000135679', 'MKI67': 'ENSG00000148773',
}
import pyreadr
df = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
df = df[list(df.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in df.index}
cols = [str(c) for c in df.columns]

missing = [g for g, e in ENS.items() if e not in idx]
L('S33 S32-a\' 「评分 A ↔ mTORC2 轴」混杂核查')
L('=' * 78)
L('基因可得性：%s' % ('全部命中 ✓' if not missing else '★缺失 %s' % missing))
L('')

E = {}
for g, e in ENS.items():
    row = idx.get(e)
    if not row:
        continue
    vals = {}
    for c in cols:
        if c in score:
            try:
                vals[c] = float(df.loc[row, c])
            except Exception:
                pass
    E[g] = vals

# ------------------------------------------------------- 偏相关工具
def partial_corr(x, y, zs):
    """控制一个或多个协变量后的偏相关（残差化）"""
    X = np.asarray(x, float); Y = np.asarray(y, float)
    Z = np.column_stack([np.ones(len(X))] + [np.asarray(z, float) for z in zs])
    k = Z.shape[1] - 1
    if len(X) < k + 4:
        return None
    def resid(a):
        beta, *_ = np.linalg.lstsq(Z, a, rcond=None)
        return a - Z @ beta
    rx, ry = resid(X), resid(Y)
    n = len(X)
    if rx.std() == 0 or ry.std() == 0:
        return {'n': n, 'partial_r': None, 'p': None, 'note': '零变异，不可估计'}
    r, _ = stats.pearsonr(rx, ry)
    d = n - 2 - k
    if d < 1 or abs(r) >= 1:
        return None
    t = r * math.sqrt(d / (1 - r * r))
    p = float(2 * stats.t.sf(abs(t), d))
    return {'n': n, 'k_cov': k, 'partial_r': round(float(r), 4), 'p': round(p, 5)}

# ------------------------------------------------- 主分析：raw + 三种校正
L('【主表】评分 A 与轴基因的关联：raw / 控制纯度 / 控制纯度+CN / 控制纯度+脂肪含量')
L('%-9s %5s %16s %20s %20s %22s' % ('gene', 'n', 'raw r (P)', 'adj purity', 'adj purity+CN', 'adj purity+ADIPOQ'))
L('-' * 100)

def build(g):
    x, y, pu, cn, ad = [], [], [], [], []
    for c, s in score.items():
        if c not in E.get(g, {}):
            continue
        r0 = cli.get(c)
        if not r0:
            continue
        p = num(r0.get('paper_purity')); k = num(r0.get('paper_CN.cluster'))
        a = E.get('ADIPOQ', {}).get(c)
        x.append(E[g][c]); y.append(s); pu.append(p); cn.append(k); ad.append(a)
    return x, y, pu, cn, ad

main = {}
for g in ['RICTOR', 'MTOR', 'CD36', 'PPARG', 'FABP4', 'CPT1A', 'PLIN1']:
    if g not in E:
        continue
    x, y, pu, cn, ad = build(g)
    m = [i for i, v in enumerate(pu) if v is not None]
    if len(m) < 10:
        L('%-9s 样本不足 (%d)' % (g, len(m))); continue
    xa = [x[i] for i in m]; ya = [y[i] for i in m]; za = [pu[i] for i in m]
    cna = [cn[i] for i in m]; ada = [ad[i] for i in m]
    raw_r, raw_p = stats.pearsonr(xa, ya)
    rec = {'n': len(m), 'raw_r': round(float(raw_r), 4), 'raw_p': round(float(raw_p), 5)}

    # (1) 控制纯度
    rec['adj_purity'] = partial_corr(xa, ya, [za])
    # (2) 控制纯度 + CN 簇（CN 有缺失则只在有值子集上做）
    m2 = [i for i in range(len(xa)) if cna[i] is not None]
    rec['adj_purity_CN'] = partial_corr([xa[i] for i in m2], [ya[i] for i in m2],
                                        [[za[i] for i in m2], [cna[i] for i in m2]]) if len(m2) >= 12 else None
    rec['n_purity_CN'] = len(m2)
    # (3) 控制纯度 + ADIPOQ（脂肪含量对照：若仍显著 → 非单纯脂肪含量）
    m3 = [i for i in range(len(xa)) if ada[i] is not None]
    rec['adj_purity_ADIPOQ'] = partial_corr([xa[i] for i in m3], [ya[i] for i in m3],
                                            [[za[i] for i in m3], [ada[i] for i in m3]]) if len(m3) >= 12 else None
    main[g] = rec

    def fmt(d):
        if not d or d.get('partial_r') is None:
            return '(不可估计)'
        return 'r=%+.4f P=%.3g' % (d['partial_r'], d['p'])
    L('%-9s %5d %16s %20s %20s %22s' % (
        g, len(m), 'r=%+.4f P=%.3g' % (raw_r, raw_p),
        fmt(rec['adj_purity']), fmt(rec['adj_purity_CN']), fmt(rec['adj_purity_ADIPOQ'])))
rep['main'] = main
L('')

# ------------------------------------------------------- 留一法（LOO）
L('【留一法】逐个剔除样本后 raw r 的稳定性（n-1）')
L('%-9s %10s %10s %12s %10s %10s  %s' % ('gene', 'r_full', 'r_min', 'r_max', 'sign_flip', 'P_max', '结论'))
L('-' * 96)
loo = {}
for g in ['RICTOR', 'MTOR', 'CD36']:
    if g not in main:
        continue
    x, y, pu, cn, ad = build(g)
    m = [i for i, v in enumerate(pu) if v is not None]
    xs = [x[i] for i in m]; ys = [y[i] for i in m]
    rs, ps = [], []
    for j in range(len(xs)):
        xx = [v for k, v in enumerate(xs) if k != j]
        yy = [v for k, v in enumerate(ys) if k != j]
        r, p = stats.pearsonr(xx, yy)
        rs.append(r); ps.append(p)
    rs = np.array(rs); ps = np.array(ps)
    flip = bool((rs.min() < 0) or (rs.max() > 0) and (rs.min() * rs.max() < 0))
    concl = '★稳健（全部同号且 P_max<0.05）' if (not flip and ps.max() < 0.05) else \
            ('方向一致但个别不显著' if not flip else '★存在翻转 → 不稳健')
    loo[g] = {'r_full': round(float(rs.mean()), 4), 'r_min': round(float(rs.min()), 4),
              'r_max': round(float(rs.max()), 4), 'sign_flip': flip,
              'p_max': round(float(ps.max()), 5), 'conclusion': concl}
    L('%-9s %10.4f %10.4f %10.4f %12s %10.3g  %s' % (
        g, main[g]['raw_r'], rs.min(), rs.max(), str(flip), ps.max(), concl))
rep['leave_one_out'] = loo
L('')

# --------------------------------------------- CN 簇分层（评分 vs 轴基因在簇内是否仍相关）
L('【CN 簇内关联】按 paper_CN.cluster 分层，簇内再测 r（簇内纯度同质，可去 CN 混杂）')
L('%-9s %6s %6s %16s' % ('gene', 'CN', 'n', 'within-cluster r (P)'))
L('-' * 46)
cn_tab = {}
for g in ['RICTOR', 'MTOR', 'CD36']:
    if g not in E:
        continue
    x, y, pu, cn, ad = build(g)
    grp = {}
    for i in range(len(x)):
        if cn[i] is None:
            continue
        grp.setdefault(int(cn[i]), []).append(i)
    cn_tab[g] = {}
    for k in sorted(grp):
        ii = grp[k]
        if len(ii) < 8:
            L('%-9s %6d %6d  (过小，跳过)' % (g, k, len(ii)))
            continue
        r, p = stats.pearsonr([x[i] for i in ii], [y[i] for i in ii])
        cn_tab[g][k] = {'n': len(ii), 'r': round(float(r), 4), 'p': round(float(p), 5)}
        L('%-9s %6d %6d  r=%+.4f  P=%.3g' % (g, k, len(ii), r, p))
rep['within_CN_cluster'] = cn_tab
L('')

# --------------------------------------- 评分自身 vs 纯度 / CN（混杂源画像）
L('【混杂源画像】评分 A 自身与纯度、CN 簇、脂肪含量的关系')
def pair(a_map, b_get, tag):
    xs, ys = [], []
    for c, s in score.items():
        b = b_get(c)
        if b is None:
            continue
        xs.append(s); ys.append(b)
    if len(xs) < 10:
        L('  %-24s 样本不足 (%d)' % (tag, len(xs))); return None
    r, p = stats.pearsonr(xs, ys)
    rs, ps = stats.spearmanr(xs, ys)
    L('  %-24s n=%d  Pearson r=%+.4f P=%.3g | Spearman rho=%+.4f P=%.3g' % (tag, len(xs), r, p, rs, ps))
    return {'n': len(xs), 'pearson': [round(float(r), 4), round(float(p), 5)],
            'spearman': [round(float(rs), 4), round(float(ps), 5)]}
prof = {}
prof['score_vs_purity'] = pair(score, lambda c: num(cli.get(c, {}).get('paper_purity')), 'score vs paper_purity')
prof['score_vs_CNcluster'] = pair(score, lambda c: num(cli.get(c, {}).get('paper_CN.cluster')), 'score vs paper_CN.cluster')
prof['score_vs_ADIPOQ'] = pair(score, lambda c: E.get('ADIPOQ', {}).get(c), 'score vs ADIPOQ')
prof['score_vs_PLIN1'] = pair(score, lambda c: E.get('PLIN1', {}).get(c), 'score vs PLIN1')
prof['score_vs_MKI67'] = pair(score, lambda c: E.get('MKI67', {}).get(c), 'score vs MKI67 (增殖对照)')
prof['score_vs_MDM2'] = pair(score, lambda c: E.get('MDM2', {}).get(c), 'score vs MDM2 (扩增对照)')
rep['confounder_profile'] = prof
L('')

# ------------------------------------------------------- 生存相关（为 6F 打底）
L('【6F 前置】评分 / RICTOR / CD36 与 OS 的单变量 Cox（n=50）')
try:
    from lifelines import CoxPHFitter
    import pandas as pd
    rows = []
    for c, s in score.items():
        r0 = cli.get(c)
        if not r0:
            continue
        t = num(r0.get('OS_time')); e = num(r0.get('OS_status'))
        if t is None or e is None or t <= 0:
            continue
        rows.append({'S': s, 'RICTOR': E.get('RICTOR', {}).get(c), 'CD36': E.get('CD36', {}).get(c),
                     'MTOR': E.get('MTOR', {}).get(c), 'T': t, 'E': int(e)})
    d = pd.DataFrame(rows).dropna()
    L('  n=%d, events=%d' % (len(d), int(d['E'].sum())))
    cph_ = {}
    for v in ['S', 'RICTOR', 'CD36', 'MTOR']:
        cc = CoxPHFitter().fit(d[[v, 'T', 'E']], 'T', 'E')
        s_ = cc.summary.loc[v]
        cph_[v] = {'HR': round(float(s_['exp(coef)']), 4), 'P': round(float(s_['p']), 5),
                   'CI95': [round(float(s_['exp(coef) lower 95%']), 4), round(float(s_['exp(coef) upper 95%']), 4)]}
        L('  %-8s HR=%.4f  95%%CI [%.3f, %.3f]  P=%.4g' % (
            v, cph_[v]['HR'], cph_[v]['CI95'][0], cph_[v]['CI95'][1], cph_[v]['P']))
    rep['cox_univariate'] = cph_
except ImportError:
    L('  lifelines 未安装（跳过；6F 用 R survival 做）')
L('')

json.dump(rep, open(os.path.join(OUT, 's33_a_prime.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
open(os.path.join(OUT, 's33_a_prime.txt'), 'w', encoding='utf-8').write('\n'.join(lines))
print('\n[写出] s33_a_prime.json / s33_a_prime.txt')
