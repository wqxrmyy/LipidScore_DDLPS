# -*- coding: utf-8 -*-
"""D15b · 评分改版的对照实验（修正版：全部变体在「定义现行评分的那一个矩阵」上重算）

关键修正：D15 首版用 TCGA_LPS_real_log2cpm.rds 重算 6 基因评分，
得到与 TCGA_Score.csv 仅 r = 0.8967 / 最大差 0.5262 的结果 —— 说明两个矩阵并非同一标准化体系。
本版一律改用 **TCGA_6genes_expr.rds**（即 score_formula.py 中已证实可精确复现 TCGA_Score.csv 的那个矩阵），
以保证「orig6 基线」= 手稿实际评分，变体之间的比较才是同基准的。

评价变量（响应侧，不受评分基准影响）：CN 簇、纯度、FNCLCC grade、核分裂率、OS、CD36 表达。
CD36 基因本身不在 6 基因矩阵内，故含 CD36 的变体只能跨基质构造，本版一律**不列为候选**，仅记录其存在。
"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

# ---------- 现行评分 ----------
score_csv = {r['Sample']: float(r['Score']) for r in
             csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cli = {r['barcode']: r for r in
       csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}

res = pyreadr.read_r(os.path.join(RAW, 'TCGA_6genes_expr.rds'))
g6 = res[list(res.keys())[0]]
genes6 = [str(g) for g in g6.index]
cols6 = [str(c) for c in g6.columns]
M = g6.values.astype(float)
rep['matrix_genes'] = genes6
rep['matrix_shape'] = list(g6.shape)


def num(x):
    if x is None: return None
    s = str(x).strip()
    if s in ('', 'NA', 'nan', '#N/A', 'Not Reported'): return None
    try:
        v = float(s); return None if math.isnan(v) else v
    except Exception: return None


def _f(v):
    out = []
    for x in v:
        try: out.append(float(x))
        except Exception: out.append(np.nan)
    return np.asarray(out, float)


def pear(x, y):
    x, y = _f(x), _f(y)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5 or len(set(x)) < 3: return None
    r, p = stats.pearsonr(x, y)
    return {'r': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def spear(x, y):
    x, y = _f(x), _f(y)
    m = ~(np.isnan(x) | np.isnan(y)); x, y = x[m], y[m]
    if len(x) < 5 or len(set(x)) < 3: return None
    r, p = stats.spearmanr(x, y)
    return {'rho': round(float(r), 4), 'p': round(float(p), 6), 'n': len(x)}


def partial(x, y, *zs):
    X, Y = _f(x), _f(y)
    Z = np.column_stack([_f(z) for z in zs])
    keep = ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z).any(axis=1))
    X, Y, Z = X[keep], Y[keep], Z[keep]
    if len(X) < 8 + Z.shape[1]: return None
    A = np.column_stack([np.ones(len(Z)), Z])
    bx, *_ = np.linalg.lstsq(A, X, rcond=None); by, *_ = np.linalg.lstsq(A, Y, rcond=None)
    rx, ry = X - A @ bx, Y - A @ by
    n = len(X); r, _ = stats.pearsonr(rx, ry); d = n - 2 - Z.shape[1]
    t = r * math.sqrt(d / max(1e-9, 1 - r * r)); p = float(2 * stats.t.sf(abs(t), d))
    return {'r': round(float(r), 4), 'p': round(p, 5), 'n': n, 'n_covars': Z.shape[1]}


def zrow(v):
    return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)


# ---------- 变体定义（仅限 6 基因矩阵内的基因）----------
VARIANTS = {
    'orig6':         ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN'],
    'v5_nofasn':     ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG'],
    'v4_core':       ['SCD', 'HMGCR', 'ACLY', 'PPARG'],
    'v4_nopparg':    ['SREBF1', 'SCD', 'HMGCR', 'ACLY'],
    'v3_hmgaclyscd': ['SCD', 'HMGCR', 'ACLY'],
    'v_scd1':        ['SCD'],
    'v_hmgcr1':      ['HMGCR'],
    'v_acly1':       ['ACLY'],
    'v_pparg1':      ['PPARG'],
    'v_fasn1':       ['FASN'],
    'v_srebf11':     ['SREBF1'],
}
rep['variants'] = VARIANTS

gi = {g: i for i, g in enumerate(genes6)}
S = {}
for v, gl in VARIANTS.items():
    miss = [g for g in gl if g not in gi]
    if miss:
        S[v] = None; rep.setdefault('variants_missing_genes', {})[v] = miss; continue
    S[v] = np.nanmean(np.vstack([zrow(M[gi[g]]) for g in gl]), axis=0)

# 基线校验：orig6 是否精确复现 TCGA_Score.csv
cur = np.array([score_csv.get(c, np.nan) for c in cols6], float)
ok = ~np.isnan(cur)
rr_, _ = stats.pearsonr(S['orig6'][ok], cur[ok])
rep['baseline_check'] = {'pearson_vs_TCGA_Score_csv': round(float(rr_), 6),
                         'max_abs_diff': round(float(np.max(np.abs(S['orig6'][ok] - cur[ok]))), 6),
                         'verdict': 'OK' if rr_ > 0.999 else 'MISMATCH'}

# ---------- 工作台 ----------
work = []
for i, c in enumerate(cols6):
    r0 = cli.get(c)
    if not r0: continue
    work.append({'sample': c, 'i': i,
                 'grade': num(r0.get('paper_FNCLCC.grade')),
                 'mitotic': num(r0.get('paper_mitotic.rate')),
                 'purity': num(r0.get('paper_purity')),
                 'cncluster': num(r0.get('paper_CN.cluster')),
                 'os_time': num(r0.get('OS_time')) or num(r0.get('paper_OS.days')),
                 'os_status': num(r0.get('OS_status'))})
rep['n_work'] = len(work)
rep['n_with_os'] = sum(1 for w in work if w['os_time'] and w['os_status'] is not None)

# CD36 表达（响应侧对照，另取一个基质 —— 已标注）
CD36 = {}
try:
    big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
    big = big[list(big.keys())[0]]
    bidx = {str(k).split('.')[0]: str(k) for k in big.index}
    row = bidx.get('ENSG00000135218')
    if row:
        for c in [str(x) for x in big.columns]:
            try: CD36[c] = float(big.loc[row, c])
            except Exception: pass
except Exception as e:
    rep['cd36_error'] = repr(e)[:150]
rep['cd36_n'] = len(CD36)

# ---------- 逐变体评价 ----------
out = {}
for v in VARIANTS:
    sc = S[v]
    if sc is None: continue
    vs = np.array([sc[w['i']] for w in work], float)
    d = {'genes': VARIANTS[v], 'n_genes': len(VARIANTS[v])}
    d['T1_r_CNcluster'] = pear([w['cncluster'] for w in work], vs)
    d['T2_partial_CN_purity'] = partial([w['cncluster'] for w in work], vs, [w['purity'] for w in work])
    d['T3_partial_CN_purity_grade_mitotic'] = partial(
        [w['cncluster'] for w in work], vs,
        [w['purity'] for w in work], [w['grade'] for w in work], [w['mitotic'] for w in work])
    d['T4_spearman_grade'] = spear([w['grade'] for w in work], vs)
    d['T5_spearman_mitotic'] = spear([w['mitotic'] for w in work], vs)
    if CD36:
        d['T9_pearson_CD36'] = pear([CD36.get(w['sample'], np.nan) for w in work], vs)
    out[v] = d

# Cox + Harrell C
try:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    import pandas as pd
    for v in VARIANTS:
        sc = S[v]
        if sc is None: continue
        rows = []
        for w in work:
            if not (w['os_time'] and w['os_status'] is not None): continue
            rows.append({'T': w['os_time'], 'E': int(w['os_status']), 's': float(sc[w['i']]),
                         'grade': w['grade'], 'mitotic': w['mitotic']})
        df = pd.DataFrame(rows)
        df['sz'] = (df['s'] - df['s'].mean()) / df['s'].std(ddof=1)
        try:
            s1 = CoxPHFitter().fit(df[['T', 'E', 'sz']], 'T', 'E').summary.loc['sz']
            out[v]['T6_Cox_unadj'] = {'HR': round(float(np.exp(s1['coef'])), 4),
                                      'CI': [round(float(np.exp(s1['coef lower 95%'])), 4),
                                             round(float(np.exp(s1['coef upper 95%'])), 4)],
                                      'p': round(float(s1['p']), 5)}
        except Exception as e:
            out[v]['T6_Cox_unadj'] = {'error': repr(e)[:100]}
        try:
            d3 = df.dropna(subset=['grade', 'mitotic'])
            s2 = CoxPHFitter().fit(d3[['T', 'E', 'sz', 'grade', 'mitotic']], 'T', 'E').summary.loc['sz']
            out[v]['T7_Cox_adj'] = {'n': int(len(d3)), 'HR': round(float(np.exp(s2['coef'])), 4),
                                    'CI': [round(float(np.exp(s2['coef lower 95%'])), 4),
                                           round(float(np.exp(s2['coef upper 95%'])), 4)],
                                    'p': round(float(s2['p']), 5)}
        except Exception as e:
            out[v]['T7_Cox_adj'] = {'error': repr(e)[:100]}
        try:
            out[v]['T8_harrell_C'] = round(float(concordance_index(df['T'], -df['s'], df['E'])), 4)
        except Exception:
            out[v]['T8_harrell_C'] = None
except Exception as e:
    rep['cox_error'] = repr(e)[:300]
rep['TCGA'] = out

# ---------- 合并 D15 首版的 GSE221492 结果（该部分基准正确，无需重算）----------
try:
    old = json.load(open(os.path.join(OUT, 'd15_score_revision.json'), encoding='utf-8'))
    rep['GSE221492'] = old['GSE221492']
    rep['GSE221492_note'] = '来自 d15_score_revision.py：GSE221492 全部变体在同一矩阵内计算，基准自洽。'
except Exception as e:
    rep['gse_error'] = repr(e)[:150]

# ---------- 裁决（阈值同 D15 首版，预先声明）----------
b = 'orig6'
bt3 = rep['TCGA'][b]['T3_partial_CN_purity_grade_mitotic']
bc = rep['TCGA'][b].get('T8_harrell_C')
bhr = rep['TCGA'][b].get('T6_Cox_unadj', {}).get('HR')
g = rep.get('GSE221492', {})
bg1 = (g.get(b) or {}).get('G1_AUC_DD_vs_WD')
verdict = {}
for v in VARIANTS:
    if v == b or v not in rep['TCGA']: continue
    d = rep['TCGA'][v]
    t3 = d['T3_partial_CN_purity_grade_mitotic']
    cv = d.get('T8_harrell_C')
    gv = (g.get(v) or {}).get('G1_AUC_DD_vs_WD')
    c = {'a_G1_auc_before': bg1.get('auc') if bg1 else None,
         'a_G1_auc_after': gv.get('auc') if gv else None,
         'a_fixes_direction': bool(bg1 and gv and bg1['auc'] < 0.5 and gv['auc'] > 0.5),
         'b_T3_before': bt3['r'] if bt3 else None,
         'b_T3_after': t3['r'] if t3 else None,
         'b_gain_pct': round((t3['r'] - bt3['r']) / abs(bt3['r']) * 100, 1) if (t3 and bt3) else None,
         'b_worth': bool(t3 and bt3 and (t3['r'] - bt3['r']) / abs(bt3['r']) >= 0.25),
         'c_C_before': bc, 'c_C_after': cv,
         'c_gain': round(cv - bc, 4) if (cv is not None and bc is not None) else None,
         'c_worth': bool(cv is not None and bc is not None and cv - bc >= 0.03),
         'd_HR_before': bhr, 'd_HR_after': d.get('T6_Cox_unadj', {}).get('HR')}
    c['worth_changing'] = bool(c['a_fixes_direction'] or c['b_worth'] or c['c_worth'])
    verdict[v] = c
rep['verdict'] = verdict

json.dump(rep, open(os.path.join(OUT, 'd15b_score_revision.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)

# ---------- 可读报告 ----------
L = []
L.append('===== D15b 基线校验 =====')
L.append(json.dumps(rep['baseline_check'], ensure_ascii=False))
L.append(f"矩阵基因: {genes6}  形状 {rep['matrix_shape']}  n_work={rep['n_work']}  有OS={rep['n_with_os']}")
L.append('')
L.append('--- TCGA（同一 6 基因矩阵基准，n=50）---')
L.append(f"{'variant':15} {'ng':>3} {'T1 r':>7} {'T2 ctrlP':>9} {'T3 full':>8} {'grade':>7} {'mitot':>7} {'HR/SD':>7} {'HRadj':>7} {'C':>7} {'CD36 r':>7}")
for v, d in rep['TCGA'].items():
    f = lambda x: (x or {}).get('r') if isinstance(x, dict) and 'r' in x else None
    L.append(f"{v:15} {d['n_genes']:>3} "
             f"{f(d.get('T1_r_CNcluster'))!s:>7} {f(d.get('T2_partial_CN_purity'))!s:>9} "
             f"{f(d.get('T3_partial_CN_purity_grade_mitotic'))!s:>8} "
             f"{(d.get('T4_spearman_grade') or {}).get('rho')!s:>7} "
             f"{(d.get('T5_spearman_mitotic') or {}).get('rho')!s:>7} "
             f"{d.get('T6_Cox_unadj',{}).get('HR')!s:>7} {d.get('T7_Cox_adj',{}).get('HR')!s:>7} "
             f"{d.get('T8_harrell_C')!s:>7} {f(d.get('T9_pearson_CD36'))!s:>7}")
L.append('')
L.append('--- GSE221492（同一矩阵基准，53 样本；DD=17/WD=17，严格配对 8 对）---')
L.append(f"{'variant':15} {'AUC DDvWD':>10} {'p':>9} {'paired d':>9} {'neg/8':>6} {'W p':>9}")
for v, d in (rep.get('GSE221492') or {}).items():
    if 'error' in d: L.append(f'{v:15} ERROR'); continue
    a = d.get('G1_AUC_DD_vs_WD') or {}; p2 = d.get('G2_paired') or {}
    L.append(f"{v:15} {a.get('auc')!s:>10} {a.get('p')!s:>9} {p2.get('mean_delta')!s:>9} "
             f"{str(p2.get('n_neg'))+'/'+str(p2.get('n_pairs')):>6} {p2.get('wilcoxon_p')!s:>9}")
L.append('')
L.append('--- 裁决表（阈值预先声明：a 修好方向 / b T3+25% / c C+0.03）---')
L.append(json.dumps(verdict, ensure_ascii=False, indent=1))
open(os.path.join(OUT, 'd15b_report.txt'), 'w', encoding='utf-8').write('\n'.join(L))
print('OK')
