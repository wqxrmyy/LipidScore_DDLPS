# -*- coding: utf-8 -*-
"""方向性总检：评分究竟是「脂代谢轴」还是「去分化/增殖轴的读数」？
① TCGA：评分 vs OS（Cox）；② 评分 vs FNCLCC grade / 核分裂率 / 坏死（分化—增殖代理）；
③ 卡 D（CN 簇）是否经 grade/核分裂率校正后仍成立；
④ 评分 vs 纯度 / 年龄 / 性别 的残余混杂。
"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
cli = {r['barcode']: r for r in cl}


def num(x):
    if x is None: return None
    s = str(x).strip()
    if s in ('', 'NA', 'nan', '#N/A', 'Not Reported'): return None
    try:
        v = float(s); return None if math.isnan(v) else v
    except Exception: return None


def _f(v):
    """None/'' -> nan，其余转 float，保证 numpy 可处理"""
    out = []
    for x in v:
        if x is None:
            out.append(np.nan); continue
        try:
            out.append(float(x))
        except Exception:
            out.append(np.nan)
    return np.asarray(out, float)


def pear(x, y):
    x, y = _f(x), _f(y)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5: return None
    r, p = stats.pearsonr(x, y)
    return {'r': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def spear(x, y):
    x, y = _f(x), _f(y)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5: return None
    r, p = stats.spearmanr(x, y)
    return {'rho': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def partial(x, y, *zs):
    X = _f(x); Y = _f(y)
    Z = np.column_stack([_f(z) for z in zs])
    keep = ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z).any(axis=1))
    X, Y, Z = X[keep], Y[keep], Z[keep]
    if len(X) < 6 + Z.shape[1]: return None
    A = np.column_stack([np.ones(len(Z)), Z])
    bxx, *_ = np.linalg.lstsq(A, X, rcond=None); byy, *_ = np.linalg.lstsq(A, Y, rcond=None)
    rx, ry = X - A @ bxx, Y - A @ byy
    n = len(X); r, _ = stats.pearsonr(rx, ry); d = n - 2 - Z.shape[1]
    t = r * math.sqrt(d / max(1e-9, 1 - r * r)); p = float(2 * stats.t.sf(abs(t), d))
    return {'r': round(float(r), 4), 'p': round(p, 5), 'n': n, 'n_covars': Z.shape[1]}


rep = {}
rec = []
for s, sv in score.items():
    r0 = cli.get(s)
    if not r0: continue
    rec.append({
        'sample': s, 'score': sv,
        'grade': num(r0.get('paper_FNCLCC.grade')),
        'mitotic': num(r0.get('paper_mitotic.rate')),
        'necrosis_pct': (lambda s: (1.0 if s and s.strip().startswith('<') else (2.0 if s and s.strip() else None)))(r0.get('paper_Necrosis.percentage')),
        'necrosis_present': (1.0 if (r0.get('paper_Necrosis.present') or '').strip().lower() == 'yes' else (0.0 if (r0.get('paper_Necrosis.present') or '').strip().lower() == 'no' else None)),
        'purity': num(r0.get('paper_purity')),
        'cncluster': num(r0.get('paper_CN.cluster')),
        'ploidy': num(r0.get('paper_ploidy')),
        'wgd': num(r0.get('paper_Genome.doublings')),
        'cin': num(r0.get('paper_CIN')),
        'cancdna': num(r0.get('paper_Cancer.DNA.fraction')),
        'age': num(r0.get('paper_age.at.diagnosis')),
        'os_time': num(r0.get('OS_time')) or num(r0.get('paper_OS.days')),
        'os_status': num(r0.get('OS_status')),
        'histo': (r0.get('paper_histology') or '').strip(),
        'short': (r0.get('paper_short.histo') or '').strip(),
        'gender': (r0.get('paper_gender') or '').strip(),
    })
rep['n'] = len(rec)
rep['histology_values'] = sorted({r['histo'] for r in rec})
rep['short_histo_values'] = sorted({r['short'] for r in rec})
rep['grade_dist'] = {str(k): sum(1 for r in rec if r['grade'] == k) for k in sorted({r['grade'] for r in rec if r['grade'] is not None})}

# ---------- ① 评分 vs 分化/增殖代理 ----------
S = [r['score'] for r in rec]
rep['score_vs_differentiation'] = {
    'FNCLCC_grade': {'pearson': pear([r['grade'] for r in rec], S), 'spearman': spear([r['grade'] for r in rec], S)},
    'mitotic_rate': {'pearson': pear([r['mitotic'] for r in rec], S), 'spearman': spear([r['mitotic'] for r in rec], S)},
    'necrosis_pct': {'pearson': pear([r['necrosis_pct'] for r in rec], S), 'spearman': spear([r['necrosis_pct'] for r in rec], S)},
    'necrosis_present': {'pearson': pear([r['necrosis_present'] for r in rec], S)},
    'age': {'pearson': pear([r['age'] for r in rec], S)},
}
# 按 grade 分组的评分
for g in sorted({r['grade'] for r in rec if r['grade'] is not None}):
    v = [r['score'] for r in rec if r['grade'] == g]
    if len(v) >= 2:
        rep.setdefault('score_by_grade', {})[f'grade_{int(g)}'] = {'n': len(v), 'mean': round(float(np.mean(v)), 4),
                                                                  'sd': round(float(np.std(v, ddof=1)), 4)}
# 按坏死/核分裂分组
for key, lab_ in [('necrosis_present', {1.0: 'necrosis_yes', 0.0: 'necrosis_no'})]:
    gg = {}
    for k, nm in lab_.items():
        v = [r['score'] for r in rec if r[key] == k]
        if len(v) >= 2:
            gg[nm] = {'n': len(v), 'mean': round(float(np.mean(v)), 4),
                      'sd': round(float(np.std(v, ddof=1)), 4)}
    if len(gg) == 2:
        a = [r['score'] for r in rec if r[key] == 1.0]; b = [r['score'] for r in rec if r[key] == 0.0]
        if len(a) >= 2 and len(b) >= 2:
            t, p = stats.ttest_ind(a, b, equal_var=False)
            gg['welch_p'] = round(float(p), 5)
    rep['score_by_' + key] = gg
# 核分裂率中位分组
mi = [r for r in rec if r['mitotic'] is not None]
if len(mi) >= 10:
    md = float(np.median([r['mitotic'] for r in mi]))
    hi = [r['score'] for r in mi if r['mitotic'] > md]; lo = [r['score'] for r in mi if r['mitotic'] <= md]
    t, p = stats.ttest_ind(hi, lo, equal_var=False)
    rep['score_by_mitotic_median'] = {'cut': md, 'n_high': len(hi), 'mean_high': round(float(np.mean(hi)), 4),
                                      'n_low': len(lo), 'mean_low': round(float(np.mean(lo)), 4),
                                      'welch_p': round(float(p), 5)}

# ---------- ② 评分 vs OS（Cox）----------
try:
    from lifelines import CoxPHFitter
    import pandas as pd
    df = pd.DataFrame([{'T': r['os_time'], 'E': r['os_status'], 'score': r['score'],
                        'score_z': (r['score'] - np.mean(S)) / np.std(S, ddof=1),
                        'grade': r['grade'], 'mitotic': r['mitotic'], 'purity': r['purity']}
                       for r in rec if r['os_time'] and r['os_status'] is not None])
    df['E'] = df['E'].astype(int)
    rep['OS_n'] = int(len(df)); rep['OS_events'] = int(df['E'].sum())
    cph = CoxPHFitter().fit(df[['T', 'E', 'score_z']], 'T', 'E')
    s = cph.summary.loc['score_z']
    rep['Cox_score_perSD'] = {'HR': round(float(np.exp(s['coef'])), 4),
                              'CI': [round(float(np.exp(s['coef lower 95%'])), 4),
                                     round(float(np.exp(s['coef upper 95%'])), 4)],
                              'p': round(float(s['p']), 5)}
    # 二分位
    med = float(np.median(df['score']))
    df['hi'] = (df['score'] > med).astype(int)
    c2 = CoxPHFitter().fit(df[['T', 'E', 'hi']], 'T', 'E')
    s2 = c2.summary.loc['hi']
    rep['Cox_score_highVSlow'] = {'HR': round(float(np.exp(s2['coef'])), 4),
                                  'CI': [round(float(np.exp(s2['coef lower 95%'])), 4),
                                         round(float(np.exp(s2['coef upper 95%'])), 4)],
                                  'p': round(float(s2['p']), 5),
                                  'n_high': int(df['hi'].sum()), 'n_low': int((1 - df['hi']).sum()),
                                  'median_cut': round(med, 4)}
    # 校正 grade / mitotic
    d3 = df.dropna(subset=['grade', 'mitotic'])
    if len(d3) >= 15:
        c3 = CoxPHFitter().fit(d3[['T', 'E', 'score_z', 'grade', 'mitotic']], 'T', 'E')
        s3 = c3.summary.loc['score_z']
        rep['Cox_score_adj_grade_mitotic'] = {'n': int(len(d3)), 'HR': round(float(np.exp(s3['coef'])), 4),
                                              'CI': [round(float(np.exp(s3['coef lower 95%'])), 4),
                                                     round(float(np.exp(s3['coef upper 95%'])), 4)],
                                              'p': round(float(s3['p']), 5)}
except Exception as e:
    rep['Cox_error'] = repr(e)[:300]

# ---------- ③ 卡 D 是否经分化/增殖校正后仍成立 ----------
rep['D_robustness_extra'] = {
    'score~CNcluster_ctrl_purity': partial([r['cncluster'] for r in rec], S, [r['purity'] for r in rec]),
    'score~CNcluster_ctrl_grade': partial([r['cncluster'] for r in rec], S, [r['grade'] for r in rec]),
    'score~CNcluster_ctrl_mitotic': partial([r['cncluster'] for r in rec], S, [r['mitotic'] for r in rec]),
    'score~CNcluster_ctrl_grade_mitotic_purity': partial(
        [r['cncluster'] for r in rec], S,
        [r['grade'] for r in rec], [r['mitotic'] for r in rec], [r['purity'] for r in rec]),
    'score~WGD_ctrl_grade_mitotic': partial([r['wgd'] for r in rec], S, [r['grade'] for r in rec], [r['mitotic'] for r in rec]),
    'CNcluster_vs_grade': spear([r['cncluster'] for r in rec], [r['grade'] for r in rec]),
    'CNcluster_vs_mitotic': spear([r['cncluster'] for r in rec], [r['mitotic'] for r in rec]),
    'CNcluster_vs_purity': pear([r['cncluster'] for r in rec], [r['purity'] for r in rec]),
    'WGD_vs_grade': spear([r['wgd'] for r in rec], [r['grade'] for r in rec]),
}
# 各 CN 簇的 grade / mitotic / score 画像
byc = {}
for r in rec:
    if r['cncluster'] is None: continue
    byc.setdefault(int(r['cncluster']), []).append(r)
prof = {}
for c, g in sorted(byc.items()):
    ga = [x['grade'] for x in g if x['grade'] is not None]
    mi = [x['mitotic'] for x in g if x['mitotic'] is not None]
    pu = [x['purity'] for x in g if x['purity'] is not None]
    prof[f'cluster_{c}'] = {'n': len(g), 'score_mean': round(float(np.mean([x['score'] for x in g])), 4),
                            'grade_mean': round(float(np.mean(ga)), 3) if ga else None, 'grade_n': len(ga),
                            'mitotic_mean': round(float(np.mean(mi)), 3) if mi else None,
                            'purity_mean': round(float(np.mean(pu)), 3) if pu else None}
rep['CNcluster_profile'] = prof

# ---------- ④ 12q 表达剂量 vs 分化（解释 A′ 为何在 TCGA 里"看似成立"）----------
dfg = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
dfg = dfg[list(dfg.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in dfg.index}
ENS = {'MDM2': 'ENSG00000135679', 'FRS2': 'ENSG00000166225', 'CD36': 'ENSG00000135218',
       'PPARG': 'ENSG00000132170', 'CD8A': 'ENSG00000153563', 'FASN': 'ENSG00000169710',
       'SCD': 'ENSG00000099194', 'SREBF1': 'ENSG00000072310'}
gex = {}
for g, e in ENS.items():
    row = idx.get(e)
    if not row: continue
    gex[g] = {r['sample']: float(dfg.loc[row, r['sample']]) for r in rec
              if r['sample'] in set(dfg.columns)}
t = {}
for g, d in gex.items():
    xs = [r['score'] for r in rec if r['sample'] in d]
    ys = [d[r['sample']] for r in rec if r['sample'] in d]
    t[g] = {'spearman_vs_score': spear(xs, ys)}
    # 与 grade 的关系
    xg = [r['grade'] for r in rec if r['sample'] in d and r['grade'] is not None]
    yg = [d[r['sample']] for r in rec if r['sample'] in d and r['grade'] is not None]
    t[g]['vs_grade'] = spear(xg, yg)
rep['genes_vs_score_and_grade'] = t

json.dump(rep, open(os.path.join(OUT, 'd11_direction_check.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
