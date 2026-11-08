# -*- coding: utf-8 -*-
"""D15 · 评分改版是否值得？—— 多变体 × 多终点的对照实验

变体（全部按「跨样本 z 分 → 取均值」的同一原式）：
  orig6         : SREBF1 SCD HMGCR ACLY PPARG FASN   （现行评分）
  v5_nofasn     : 去 FASN
  v4_core       : 去 FASN + SREBF1        （SCD HMGCR ACLY PPARG）
  v3_hmgaclyscd : SCD HMGCR ACLY
  v_scd1        : SCD 单基因
  v_cd36        : CD36 单基因              （跨队列唯一复现的关联）
  v_pparg_cd36  : PPARG + CD36             （卡 C 的靶轴）

终点：
  TCGA（DDLPS-only, n=50）
   T1 r(score, CN 簇) 原始
   T2 r(score, CN 簇) 控制纯度                     ← 卡 D 主指标
   T3 r(score, CN 簇) 控制 纯度+grade+核分裂        ← 卡 D 全控指标
   T4 评分 vs FNCLCC grade（斯皮尔曼）              ← 「增殖读数」风险
   T5 评分 vs 核分裂率
   T6 Cox OS HR / SD（未校正）
   T7 Cox OS HR / SD（校正 grade+核分裂）
   T8 Harrell C（未校正）
   T9 评分 vs CD36 表达

  GSE221492（bulk LPS atlas, 53 样本；DD=17, WD=17, 严格配对 8 对）
   G1 评分对 DDLPS_DD vs DDLPS_WD 的 AUC           ← 现行 0.2491（反向！）
   G2 配对 8 对的均值差 / 负号数 / Wilcoxon p
   G3 评分 vs MDM2（A′，现行 ρ=+0.121 NS）
   G4 评分 vs CD36                                （现行 ρ=+0.509）
   G5 评分 vs CD8A                                （现行 ρ=+0.056 NS）

裁决口径（预先声明，避免事后挑指标）：
  改版「值得」须至少满足其一 ——
    (a) G1 的 AUC 由 <0.5 翻转为 >0.5（修好跨队列方向矛盾）；
    (b) T3 相对 orig6 提升 ≥ 25%（卡 D 主指标实质性增强）；
    (c) T8 Harrell C 提升 ≥ 0.03（判别度实质性增强）。
  全部不满足 → 不改版（结论：性价比不足）。
"""
import os, gzip, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
GZ = os.path.join(OUT, 'gse221492', 'GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz')

ENS = {'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
       'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710',
       'CD36': 'ENSG00000135218', 'MDM2': 'ENSG00000135679', 'CD8A': 'ENSG00000153563'}

VARIANTS = {
    'orig6':        ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN'],
    'v5_nofasn':    ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG'],
    'v4_core':      ['SCD', 'HMGCR', 'ACLY', 'PPARG'],
    'v3_hmgaclyscd': ['SCD', 'HMGCR', 'ACLY'],
    'v_scd1':       ['SCD'],
    'v_cd36':       ['CD36'],
    'v_pparg_cd36': ['PPARG', 'CD36'],
}
rep = {'variants': {k: v for k, v in VARIANTS.items()}}


# ============================ 工具 ============================
def _f(v):
    out = []
    for x in v:
        try:
            out.append(float(x))
        except Exception:
            out.append(np.nan)
    return np.asarray(out, float)


def num(x):
    if x is None: return None
    s = str(x).strip()
    if s in ('', 'NA', 'nan', '#N/A', 'Not Reported'): return None
    try:
        v = float(s); return None if math.isnan(v) else v
    except Exception:
        return None


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


def auc(x1, x0):
    x1 = _f(x1); x0 = _f(x0)
    x1 = x1[~np.isnan(x1)]; x0 = x0[~np.isnan(x0)]
    if len(x1) < 2 or len(x0) < 2: return None
    U = 0.0
    for a in x1:
        U += (x0 < a).sum() + 0.5 * (x0 == a).sum()
    a_ = U / (len(x1) * len(x0))
    try:
        _, p = stats.mannwhitneyu(x1, x0, alternative='two-sided')
        p = round(float(p), 6)
    except Exception:
        p = None
    return {'auc': round(float(a_), 4), 'n1': len(x1), 'n0': len(x0), 'p': p,
            'mean1': round(float(np.mean(x1)), 4), 'mean0': round(float(np.mean(x0)), 4)}


# ============================ TCGA ============================
score_csv = {r['Sample']: float(r['Score']) for r in
             csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cli = {r['barcode']: r for r in
       csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}

# --- 表达矩阵 ---
g6 = pyreadr.read_r(os.path.join(RAW, 'TCGA_6genes_expr.rds'))
g6 = g6[list(g6.keys())[0]]
cols6 = [str(c) for c in g6.columns]

big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
bidx = {str(i).split('.')[0]: str(i) for i in big.index}

T = {}   # gene -> {sample: expr}
for g, e in ENS.items():
    row = bidx.get(e)
    if row is None: continue
    d = {}
    for c in [str(x) for x in big.columns]:
        try:
            d[c] = float(big.loc[row, c])
        except Exception:
            pass
    T[g] = d
rep['tcga_genes_loaded'] = sorted(T.keys())

# 现行评分（TCGA_Score.csv）与由大矩阵重算的 orig6 是否一致
samples = [c for c in cols6 if c in score_csv and c in cli]
cur = _f([score_csv[s] for s in samples])


def zsc(vals):
    v = _f(vals); m = np.nanmean(v); sd = np.nanstd(v, ddof=1)
    return (v - m) / sd


def tcga_score(genes, ss):
    zs = [zsc([T[g].get(s, np.nan) for s in ss]) for g in genes if g in T]
    return np.nanmean(np.vstack(zs), axis=0)


recomputed_orig = tcga_score(VARIANTS['orig6'], samples)
mm = ~(np.isnan(recomputed_orig) | np.isnan(cur))
rr, _ = stats.pearsonr(recomputed_orig[mm], cur[mm])
rep['tcga_recompute_check'] = {
    'pearson_vs_TCGA_Score_csv': round(float(rr), 6),
    'max_abs_diff': round(float(np.max(np.abs(recomputed_orig[mm] - cur[mm]))), 6),
    'verdict': 'OK' if rr > 0.999 else 'MISMATCH'}

# --- TCGA 临床工作台 ---
work = []
for i, s in enumerate(samples):
    r0 = cli[s]
    work.append({
        'sample': s, 'cur': cur[i], 'recomp': recomputed_orig[i],
        'grade': num(r0.get('paper_FNCLCC.grade')),
        'mitotic': num(r0.get('paper_mitotic.rate')),
        'purity': num(r0.get('paper_purity')),
        'cncluster': num(r0.get('paper_CN.cluster')),
        'os_time': num(r0.get('OS_time')) or num(r0.get('paper_OS.days')),
        'os_status': num(r0.get('OS_status')),
    })
rep['tcga_n'] = len(work)

TCGA = {}
for vname, genes in VARIANTS.items():
    sc = tcga_score(genes, samples)
    d = {'genes': genes, 'n_genes': len(genes)}
    d['T1_r_CNcluster'] = pear([w['cncluster'] for w in work], sc)
    d['T2_partial_CN_purity'] = partial([w['cncluster'] for w in work], sc, [w['purity'] for w in work])
    d['T3_partial_CN_purity_grade_mitotic'] = partial(
        [w['cncluster'] for w in work], sc,
        [w['purity'] for w in work], [w['grade'] for w in work], [w['mitotic'] for w in work])
    d['T4_spearman_grade'] = spear([w['grade'] for w in work], sc)
    d['T5_spearman_mitotic'] = spear([w['mitotic'] for w in work], sc)
    d['T9_pearson_CD36'] = pear([T['CD36'].get(w['sample'], np.nan) for w in work], sc)
    TCGA[vname] = d
    d['_score'] = sc

# --- Cox + Harrell C ---
try:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    import pandas as pd
    base = [w for w in work if w['os_time'] and w['os_status'] is not None]
    for vname in VARIANTS:
        sc = TCGA[vname]['_score']
        idx = {w['sample']: i for i, w in enumerate(work)}
        rows = [{'T': w['os_time'], 'E': int(w['os_status']), 's': sc[idx[w['sample']]],
                 'grade': w['grade'], 'mitotic': w['mitotic']}
                for w in base if not np.isnan(sc[idx[w['sample']]])]
        df = pd.DataFrame(rows)
        df['sz'] = (df['s'] - df['s'].mean()) / df['s'].std(ddof=1)
        try:
            c1 = CoxPHFitter().fit(df[['T', 'E', 'sz']], 'T', 'E')
            s1 = c1.summary.loc['sz']
            TCGA[vname]['T6_Cox_unadj'] = {
                'HR': round(float(np.exp(s1['coef'])), 4),
                'CI': [round(float(np.exp(s1['coef lower 95%'])), 4),
                       round(float(np.exp(s1['coef upper 95%'])), 4)],
                'p': round(float(s1['p']), 5)}
        except Exception as e:
            TCGA[vname]['T6_Cox_unadj'] = {'error': repr(e)[:120]}
        try:
            d3 = df.dropna(subset=['grade', 'mitotic'])
            c2 = CoxPHFitter().fit(d3[['T', 'E', 'sz', 'grade', 'mitotic']], 'T', 'E')
            s2 = c2.summary.loc['sz']
            TCGA[vname]['T7_Cox_adj'] = {
                'n': int(len(d3)), 'HR': round(float(np.exp(s2['coef'])), 4),
                'CI': [round(float(np.exp(s2['coef lower 95%'])), 4),
                       round(float(np.exp(s2['coef upper 95%'])), 4)],
                'p': round(float(s2['p']), 5)}
        except Exception as e:
            TCGA[vname]['T7_Cox_adj'] = {'error': repr(e)[:120]}
        try:
            TCGA[vname]['T8_harrell_C'] = round(float(concordance_index(
                df['T'], -df['s'], df['E'])), 4)
        except Exception as e:
            TCGA[vname]['T8_harrell_C'] = None
    rep['tcga_OS_n'] = int(len(base))
    rep['tcga_OS_events'] = int(sum(w['os_status'] for w in base))
except Exception as e:
    rep['cox_error'] = repr(e)[:300]

# 分基因在 TCGA 内与 grade / 核分裂 / CN 簇的关系（诊断用）
gene_diag = {}
for g, d in T.items():
    xs = [w['sample'] for w in work if w['sample'] in d]
    if len(xs) < 10: continue
    gene_diag[g] = {
        'vs_CNcluster': spear([cli[s]['paper_CN.cluster'] for s in xs], [d[s] for s in xs]),
        'vs_grade': spear([num(cli[s].get('paper_FNCLCC.grade')) for s in xs], [d[s] for s in xs]),
        'vs_mitotic': spear([num(cli[s].get('paper_mitotic.rate')) for s in xs], [d[s] for s in xs]),
    }
rep['tcga_gene_diagnostics'] = gene_diag

for vname in TCGA: TCGA[vname].pop('_score', None)
rep['TCGA'] = TCGA

# ============================ GSE221492 ============================
samples2, mat = None, {}
WANT = sorted({g for gs in VARIANTS.values() for g in gs} | {'MDM2', 'CD8A'})
with gzip.open(GZ, 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples2 = hdr[2:]
    keep = set(WANT)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 3: continue
        sym = p[1].strip()
        if sym in keep and sym not in mat:
            mat[sym] = np.array([float(x) if x not in ('', 'NA') else np.nan for x in p[2:]])

colsum = np.zeros(len(samples2))
with gzip.open(GZ, 'rt', encoding='utf-8', errors='ignore') as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 3: continue
        try:
            colsum += np.array([float(x) if x not in ('', 'NA') else 0.0 for x in p[2:]])
        except Exception:
            pass
l2 = {}
for g in mat:
    l2[g] = np.log2(mat[g] / np.where(colsum > 0, colsum, np.nan) * 1e6 + 1)


def parse(s):
    b = s.split('_')
    return b[0].lstrip('Pp'), '_'.join(b[1:])


meta = [{'sample': s, 'patient': parse(s)[0], 'label': parse(s)[1]} for s in samples2]
lab = np.array([m['label'] for m in meta])
mA = lab == 'DDLPS_DD'; mB = lab == 'DDLPS_WD'
pairs = {}
for m in meta:
    pairs.setdefault(m['patient'], {})[m['label']] = m['sample']
strict = [(p, d['DDLPS_WD'], d['DDLPS_DD']) for p, d in sorted(pairs.items())
          if 'DDLPS_WD' in d and 'DDLPS_DD' in d]
pos = {m['sample']: i for i, m in enumerate(meta)}

G = {}
for vname, genes in VARIANTS.items():
    gs = [g for g in genes if g in l2]
    if not gs:
        G[vname] = {'error': 'no genes present'}; continue
    zs = np.vstack([zsc(l2[g]) for g in gs])
    sc = np.nanmean(zs, axis=0)
    d = {'genes_used': gs}
    d['G1_AUC_DD_vs_WD'] = auc(sc[mA], sc[mB])
    deltas = [float(sc[pos[dd]] - sc[pos[w]]) for _, w, dd in strict]
    if len(deltas) >= 3:
        _, wp = stats.wilcoxon(deltas)
        d['G2_paired'] = {'n_pairs': len(deltas), 'mean_delta': round(float(np.mean(deltas)), 4),
                          'n_neg': sum(1 for x in deltas if x < 0), 'n_pos': sum(1 for x in deltas if x > 0),
                          'wilcoxon_p': round(float(wp), 6)}
        d['G2_paired']['detail'] = [round(x, 4) for x in deltas]
    d['G3_spearman_MDM2'] = spear(sc, l2.get('MDM2', np.full(len(samples2), np.nan)))
    d['G4_spearman_CD36'] = spear(sc, l2.get('CD36', np.full(len(samples2), np.nan)))
    d['G5_spearman_CD8A'] = spear(sc, l2.get('CD8A', np.full(len(samples2), np.nan)))
    G[vname] = d
rep['GSE221492'] = G
rep['GSE221492_label_counts'] = {k: int((lab == k).sum()) for k in sorted(set(lab))}
rep['GSE221492_strict_pairs'] = len(strict)

# ============================ 裁决 ============================
base = 'orig6'
verdict = {}
b_T3 = TCGA[base]['T3_partial_CN_purity_grade_mitotic']
b_G1 = G[base]['G1_AUC_DD_vs_WD']
b_C = TCGA[base].get('T8_harrell_C')
for v in VARIANTS:
    if v == base or 'error' in G.get(v, {}): continue
    c = {}
    g1 = G[v]['G1_AUC_DD_vs_WD']
    c['a_G1_auc_before'] = b_G1['auc'] if b_G1 else None
    c['a_G1_auc_after'] = g1['auc'] if g1 else None
    c['a_fixes_direction'] = bool(b_G1 and g1 and b_G1['auc'] < 0.5 and g1['auc'] > 0.5)
    t3 = TCGA[v]['T3_partial_CN_purity_grade_mitotic']
    c['b_T3_before'] = b_T3['r'] if b_T3 else None
    c['b_T3_after'] = t3['r'] if t3 else None
    c['b_gain_pct'] = (round((t3['r'] - b_T3['r']) / abs(b_T3['r']) * 100, 1)
                       if (t3 and b_T3) else None)
    c['b_worth'] = bool(t3 and b_T3 and (t3['r'] - b_T3['r']) / abs(b_T3['r']) >= 0.25)
    cv = TCGA[v].get('T8_harrell_C')
    c['c_C_before'] = b_C; c['c_C_after'] = cv
    c['c_worth'] = bool(cv is not None and b_C is not None and cv - b_C >= 0.03)
    c['worth_changing'] = bool(c['a_fixes_direction'] or c['b_worth'] or c['c_worth'])
    verdict[v] = c
rep['verdict'] = verdict

json.dump(rep, open(os.path.join(OUT, 'd15_score_revision.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps({'tcga_recompute_check': rep['tcga_recompute_check'],
                  'tcga_n': rep['tcga_n'], 'tcga_OS': [rep.get('tcga_OS_n'), rep.get('tcga_OS_events')],
                  'TCGA': rep['TCGA'], 'GSE221492': rep['GSE221492'],
                  'verdict': rep['verdict']}, ensure_ascii=False, indent=2))
