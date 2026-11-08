# -*- coding: utf-8 -*-
"""D12 · P4 统计复核：卡 D（CN 复杂度 ↔ 脂代谢）的稳健性总检
① 复现主效应；② 留一法；③ 置换检验；④ 簇6 拆分检验；
⑤ 升级为连续 CIN 剂量（paper_CIN / ploidy / genome doublings 等此前未用的列）；
⑥ 纯度偏相关；⑦ 与 12q 扩增子的独立性再确认。
"""
import os, csv, json, math
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rng = np.random.default_rng(20260917)

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
cli = {r['barcode']: r for r in cl}


def num(x):
    if x is None: return None
    try:
        v = float(str(x).strip())
        return None if math.isnan(v) else v
    except Exception:
        return None


def pear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    r, p = stats.pearsonr(x, y)
    return round(float(r), 4), float(p), len(x)


def spear(x, y):
    r, p = stats.spearmanr(x, y)
    return round(float(r), 4), float(p)


def partial(x, y, z):
    X = np.asarray(x, float); Y = np.asarray(y, float); Z = np.asarray(z, float)
    keep = ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z))
    X, Y, Z = X[keep], Y[keep], Z[keep]
    if len(X) < 6: return None, None, len(X)
    A = np.column_stack([np.ones_like(Z), Z])
    bx, *_ = np.linalg.lstsq(A, X, rcond=None); by, *_ = np.linalg.lstsq(A, Y, rcond=None)
    rx, ry = X - A @ bx, Y - A @ by
    n = len(X); r, _ = stats.pearsonr(rx, ry); d = n - 3
    t = r * math.sqrt(d / (1 - r * r)); p = float(2 * stats.t.sf(abs(t), d))
    return round(float(r), 4), round(p, 5), n


rep = {}

# ---------- 0. 建立工作表 ----------
VARS = ['paper_CN.cluster', 'paper_CIN', 'paper_ploidy', 'paper_Genome.doublings',
        'paper_Cancer.DNA.fraction', 'paper_Subclona.genome.fraction',
        'paper_Mutation.Load', 'paper_avg.CN.12q13.15', 'paper_purity',
        'paper_leukocyte.fraction', 'paper_mitotic.rate', 'paper_Necrosis.percentage',
        'paper_FNCLCC.grade']
rec = []
for s, sv in score.items():
    r0 = cli.get(s)
    if not r0: continue
    d = {'sample': s, 'score': sv}
    for v in VARS:
        d[v] = num(r0.get(v))
    rec.append(d)
rep['n_samples'] = len(rec)

# 描述统计 + 与评分的相关（未校正 / 纯度校正）
desc = {}
for v in VARS:
    xs = [r[v] for r in rec if r[v] is not None]
    srs = [(r[v], r['score']) for r in rec if r[v] is not None]
    if len(xs) < 8: 
        desc[v] = {'n_obs': len(xs)}; continue
    d = {'n_obs': len(xs), 'n_unique': len(set(xs)),
         'min': round(min(xs), 4), 'max': round(max(xs), 4),
         'mean': round(float(np.mean(xs)), 4), 'sd': round(float(np.std(xs, ddof=1)), 4),
         'iqr': round(float(np.percentile(xs, 75) - np.percentile(xs, 25)), 4)}
    if len(set(xs)) > 2:
        r_, p_, n_ = pear([a for a, _ in srs], [b for _, b in srs])
        rs, ps = spear([a for a, _ in srs], [b for _, b in srs])
        d['pearson'] = {'r': r_, 'p': round(p_, 6), 'n': n_}
        d['spearman'] = {'rho': rs, 'p': round(ps, 6)}
    # 纯度校正
    tri = [(r[v], r['score'], r['paper_purity']) for r in rec
           if r[v] is not None and r['paper_purity'] is not None]
    if len(tri) >= 10 and len(set(t[0] for t in tri)) > 2:
        pr, pp, pn = partial([t[0] for t in tri], [t[1] for t in tri], [t[2] for t in tri])
        if pr is not None:
            d['partial_r_ctrl_purity'] = {'r': pr, 'p': pp, 'n': pn}
    desc[v] = d
rep['descriptives_and_correlations'] = desc

# ---------- 1. 卡 D 主效应复现 + 留一法 ----------
d_cn = [r for r in rec if r['paper_CN.cluster'] is not None]
x = [r['paper_CN.cluster'] for r in d_cn]; y = [r['score'] for r in d_cn]
r_main, p_main, n_main = pear(x, y)
rep['D_main_CNcluster'] = {'pearson': {'r': r_main, 'p': round(p_main, 8), 'n': n_main},
                           'spearman': dict(zip(['rho', 'p'], spear(x, y)))}

loo = []
for i in range(len(d_cn)):
    xs = [x[j] for j in range(len(x)) if j != i]; ys = [y[j] for j in range(len(y)) if j != i]
    rr, pp, _ = pear(xs, ys)
    loo.append({'dropped': d_cn[i]['sample'], 'dropped_cluster': int(d_cn[i]['paper_CN.cluster']),
                'r': rr, 'p': round(pp, 6)})
rs_loo = sorted(l['r'] for l in loo)
rep['D_leave_one_out'] = {'n_tests': len(loo), 'r_min': rs_loo[0], 'r_max': rs_loo[-1],
                          'r_median': round(float(np.median(rs_loo)), 4),
                          'n_still_p_lt_05': sum(1 for l in loo if l['p'] < 0.05),
                          'worst5': sorted(loo, key=lambda l: l['r'])[:5],
                          'best3': sorted(loo, key=lambda l: -l['r'])[:3]}

# 整簇剔除
byc = {}
for r in d_cn: byc.setdefault(int(r['paper_CN.cluster']), []).append(r)
rep['D_drop_whole_cluster'] = {}
for c in sorted(byc):
    xs = [r['paper_CN.cluster'] for r in d_cn if int(r['paper_CN.cluster']) != c]
    ys = [r['score'] for r in d_cn if int(r['paper_CN.cluster']) != c]
    rr, pp, nn = pear(xs, ys)
    rep['D_drop_whole_cluster'][f'drop_cluster_{c}'] = {'n_pairs_left': len(set(xs)), 'n': nn,
                                                        'r': rr, 'p': round(pp, 6)}

# ---------- 2. 置换检验 ----------
def perm_test(x, y, nperm=20000):
    x = np.asarray(x, float); y = np.asarray(y, float)
    obs, _, _ = pear(x, y)
    cnt = 0
    for _ in range(nperm):
        rr, _, _ = pear(rng.permutation(x), y)
        if abs(rr) >= abs(obs): cnt += 1
    return {'r_obs': obs, 'n_perm': nperm, 'p_perm': round((cnt + 1) / (nperm + 1), 6)}

rep['D_permutation_CNcluster'] = perm_test(x, y)

# Kruskal-Wallis 置换
groups = [v for v in byc.values() if len(v) >= 2]
H, pk = stats.kruskal(*[[r['score'] for r in g] for g in groups])
cnt = 0
allsc = [r['score'] for r in d_cn]; sizes = [len(g) for g in groups]; NP = 20000
for _ in range(NP):
    perm = rng.permutation(allsc)
    idx = 0; gg = []
    for s in sizes:
        gg.append(perm[idx:idx + s]); idx += s
    HH, _ = stats.kruskal(*gg)
    if HH >= H: cnt += 1
rep['D_kruskal_permutation'] = {'H': round(float(H), 4), 'p_asym': round(float(pk), 6),
                                'n_perm': NP, 'p_perm': round((cnt + 1) / (NP + 1), 6),
                                'group_sizes': {k: len(v) for k, v in sorted(byc.items())}}

# 簇6 vs 其余（Mann-Whitney）
g6 = [r['score'] for r in d_cn if int(r['paper_CN.cluster']) == 6]
go = [r['score'] for r in d_cn if int(r['paper_CN.cluster']) != 6]
if len(g6) >= 2:
    U, pu = stats.mannwhitneyu(g6, go, alternative='two-sided')
    rep['D_cluster6_vs_rest'] = {'n6': len(g6), 'mean6': round(float(np.mean(g6)), 4),
                                 'n_rest': len(go), 'mean_rest': round(float(np.mean(go)), 4),
                                 'U': round(float(U), 2), 'p': round(float(pu), 6)}
    # 簇6 内 bootstrap（4 例，只能给区间）
    bs = [float(np.mean(rng.choice(g6, size=len(g6), replace=True))) for _ in range(20000)]
    rep['D_cluster6_bootstrap_mean'] = {'ci95': [round(float(np.percentile(bs, 2.5)), 4),
                                                 round(float(np.percentile(bs, 97.5)), 4)]}

# ---------- 3. 升级为连续 CIN 剂量 ----------
rep['D_continuous_upgrade'] = {}
for v in ['paper_CIN', 'paper_ploidy', 'paper_Genome.doublings',
          'paper_Cancer.DNA.fraction', 'paper_Subclona.genome.fraction', 'paper_Mutation.Load']:
    tri = [(r[v], r['score'], r['paper_purity']) for r in rec if r[v] is not None]
    if len(tri) < 10 or len(set(t[0] for t in tri)) < 3: continue
    rr, pp, nn = pear([t[0] for t in tri], [t[1] for t in tri])
    rs, ps = spear([t[0] for t in tri], [t[1] for t in tri])
    pr, ppp, pn = partial([t[0] for t in tri], [t[1] for t in tri], [t[2] for t in tri])
    rep['D_continuous_upgrade'][v] = {'n': nn, 'pearson': {'r': rr, 'p': round(pp, 6)},
                                      'spearman': {'rho': rs, 'p': round(ps, 6)},
                                      'partial_ctrl_purity': {'r': pr, 'p': ppp, 'n': pn},
                                      'n_with_purity': sum(1 for t in tri if t[2] is not None)}

# ---------- 4. 这些 CN 变量之间的相关（是否彼此独立）----------
mm = {}
VL = ['paper_CN.cluster', 'paper_CIN', 'paper_ploidy', 'paper_Genome.doublings',
      'paper_Cancer.DNA.fraction', 'paper_Subclona.genome.fraction', 'paper_purity',
      'paper_avg.CN.12q13.15']
for i in range(len(VL)):
    for j in range(i + 1, len(VL)):
        pr = [(r[VL[i]], r[VL[j]]) for r in rec if r[VL[i]] is not None and r[VL[j]] is not None]
        if len(pr) < 10: continue
        rr, pp, nn = pear([a for a, _ in pr], [b for _, b in pr])
        mm[f'{VL[i]} ~ {VL[j]}'] = {'r': rr, 'p': round(pp, 5), 'n': nn}
rep['CN_var_matrix'] = mm

# ---------- 5. 与 12q 扩增子表达的独立性（沿用既有列）+ MDM2 表达 ----------
df = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
df = df[list(df.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in df.index}
ENS = {'MDM2': 'ENSG00000135679', 'FRS2': 'ENSG00000166225', 'CDK4': 'ENSG00000135446',
       'HMGA2': 'ENSG00000149948', 'CD36': 'ENSG00000135218', 'PPARG': 'ENSG00000132170',
       'CD8A': 'ENSG00000153563', 'FASN': 'ENSG00000169710', 'SCD': 'ENSG00000099194',
       'SREBF1': 'ENSG00000072310', 'ACLY': 'ENSG00000131473', 'HMGCR': 'ENSG00000113161'}
gex = {}
for g, e in ENS.items():
    row = idx.get(e)
    if not row: continue
    gex[g] = {}
    for r in rec:
        try:
            gex[g][r['sample']] = float(df.loc[row, r['sample']])
        except Exception:
            pass
rep['genes_loaded'] = sorted(gex.keys())

it = {}
for v in ['paper_CIN', 'paper_ploidy', 'paper_Genome.doublings', 'paper_CN.cluster']:
    sub = {}
    for g in ['MDM2', 'FRS2', 'CDK4', 'HMGA2']:
        pr = [(r[v], gex[g][r['sample']]) for r in rec
              if r[v] is not None and r['sample'] in gex.get(g, {})]
        if len(pr) < 10: continue
        rr, pp, nn = pear([a for a, _ in pr], [b for _, b in pr])
        sub[g] = {'r': rr, 'p': round(pp, 5), 'n': nn}
    it[v] = sub
rep['CN_var_vs_amplicon_expression'] = it

# ---------- 6. 评分分量循环性检查：评分 vs 其 6 个分量 ----------
comp = {}
for g in ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN', 'CD36', 'CD8A']:
    pr = [(r['score'], gex[g][r['sample']]) for r in rec if r['sample'] in gex.get(g, {})]
    if len(pr) < 10: continue
    rr, pp, nn = pear([a for a, _ in pr], [b for _, b in pr])
    comp[g] = {'r': rr, 'p': round(pp, 6), 'n': nn}
rep['score_vs_genes'] = comp

json.dump(rep, open(os.path.join(OUT, 'd12_p4_robustness.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps({k: rep[k] for k in ['n_samples', 'D_main_CNcluster', 'D_leave_one_out',
                                      'D_drop_whole_cluster', 'D_permutation_CNcluster',
                                      'D_kruskal_permutation', 'D_cluster6_vs_rest',
                                      'D_continuous_upgrade', 'score_vs_genes']},
                 ensure_ascii=False, indent=2))
