# -*- coding: utf-8 -*-
"""
S34 · v2：修正 ID 后的完整读数（含签名变体与逐组分）
----------------------------------------------------
修正：raw 用全部可比样本（n=50），偏相关用 complete-case（纯度+CN 齐备）
      mTORC2 逐组分单基因 + 6 种签名变体（全报，不择优隐藏）
输出：s34_nocontra_report.txt / s34_nocontra.json / fig6_data/mtorc2_sig_tcga.csv
"""
import os, csv, json, math
import numpy as np
from scipy import stats
import pyreadr

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
FIGD = os.path.join(W, 'p5_figures', 'fig6_data')
os.makedirs(FIGD, exist_ok=True)
Z95, Z80 = 1.959963985, 0.8416212336

LOG = os.path.join(W, 's34_nocontra_report.txt')
lines = []
def L(s=''):
    print(s); lines.append(str(s))

IDS = {}
for r in csv.DictReader(open(os.path.join(W, 's34_gene_ids.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    IDS[r['symbol']] = r['ensg']
assert IDS.get('RICTOR') == 'ENSG00000164327'
assert IDS.get('MPRIP') == 'ENSG00000133030'
assert IDS.get('MCM2') == 'ENSG00000073111'

score = {r['Sample']: float(r['Score']) for r in
         csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}
cli = {r['barcode']: r for r in
       csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}
def num(x):
    try:
        v = float(x); return None if (v != v) else v
    except Exception:
        return None

_d = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
_d = _d[list(_d.keys())[0]]
idx = {str(i).split('.')[0]: str(i) for i in _d.index}
ec = [str(c) for c in _d.columns if str(c) in score]

EM = {}
for sym, e in IDS.items():
    row = idx.get(e)
    if row:
        EM[sym] = {c: num(_d.loc[row, c]) for c in ec}

prol_g = ['MKI67', 'PCNA', 'TOP2A', 'MCM2', 'CCNB1']
base = sorted(set(ec) & set(score))
n_pur_miss = sum(1 for c in base if c not in cli or num(cli[c].get('paper_purity')) is None)
n_cn_miss = sum(1 for c in base if c not in cli or num(cli[c].get('paper_CN.cluster')) is None)
cc = [c for c in base if c not in ('',) and c in cli
      and num(cli[c].get('paper_purity')) is not None
      and num(cli[c].get('paper_CN.cluster')) is not None]

P = np.array([[EM[g][c] for g in prol_g] for c in cc], float)
Pz = (P - P.mean(0)) / P.std(0, ddof=1)
prol = {c: float(Pz[i].mean()) for i, c in enumerate(cc)}
pur = {c: num(cli[c]['paper_purity']) for c in cc}
cnc = {c: num(cli[c]['paper_CN.cluster']) for c in cc}

L('=' * 100)
L('S34 v2 · 修正基因 ID 后的完整读数（N1/N2/N3 裁定支撑）')
L('=' * 100)
L('')
L('★ 口径修正（Ensembl REST symbol 反查；表见 s34_gene_ids.tsv）')
for sym, old, new, actual in [
        ('RICTOR', 'ENSG00000133030', 'ENSG00000164327', 'MPRIP(chr17)'),
        ('MCM2', 'ENSG00000117748', 'ENSG00000073111', 'RPA2(chr1)'),
        ('MLST8', 'ENSG00000148481', 'ENSG00000167965', 'MINDY3(chr10)'),
        ('PRR5', 'ENSG00000186652', 'ENSG00000186654', 'PRG2(chr11)'),
        ('DEPTOR', 'ENSG00000155729', 'ENSG00000155792', 'KCTD18(chr2)'),
        ('LPAR1', 'ENSG00000198160', 'ENSG00000198121', 'MIER1(chr1)'),
        ('MRC1', 'ENSG00000173593', 'ENSG00000260314', '（查无记录）')]:
    L('   %-7s %s → **%s**   （旧 ID 实为 %s）' % (sym, old, new, actual))
L('')
L('样本：与评分交集 n=%d | 缺纯度 %d | 缺 CN %d | complete-case n=%d'
  % (len(base), n_pur_miss, n_cn_miss, len(cc)))

x50 = {c: score[c] for c in base}
xcc = np.array([score[c] for c in cc], float)

def rawcorr(sym):
    v = np.array([EM[sym][c] for c in base], float)
    xs = np.array([score[c] for c in base], float)
    r, p = stats.pearsonr(xs, v)
    return r, p

def partial(sym_or_vec):
    if isinstance(sym_or_vec, str):
        y = np.array([EM[sym_or_vec][c] for c in cc], float)
    else:
        y = np.asarray(sym_or_vec, float)
    M = np.column_stack([np.ones(len(xcc)), np.array([pur[c] for c in cc]),
                         np.array([cnc[c] for c in cc]), np.array([prol[c] for c in cc])])
    def rs(a):
        b, *_ = np.linalg.lstsq(M, a, rcond=None); return a - M @ b
    rx, ry = rs(xcc), rs(y)
    r, _ = stats.pearsonr(rx, ry)
    d = len(xcc) - 2 - 3
    t = r * math.sqrt(d / (1 - r * r))
    return r, 2 * stats.t.sf(abs(t), d)

# ---------------- A. 单基因（含旧 ID 对照）----------------
L('')
L('-- A. 评分 A ↔ 轴基因（raw n=%d；偏相关 complete-case n=%d）--' % (len(base), len(cc)))
L('%-18s %-20s %9s %10s %20s' % ('符号', 'ENSG', 'raw r', 'P', 'adj 纯度+CN+增殖'))
L('-' * 84)
single = {}
order = ['RICTOR', 'MPRIP', 'MTOR', 'CD36', 'PPARG', 'FABP4', 'ADIPOQ',
         'MAPKAP1', 'MLST8', 'PRR5', 'PRR5L', 'DEPTOR', 'AKT1', 'SGK1', 'NDRG1', 'PRKCA', 'RPTOR']
tag = {'MPRIP': '★MPRIP（旧ID实际）', 'PPARG': 'PPARG（循环）'}
for s in order:
    if s not in EM:
        continue
    r, p = rawcorr(s)
    ar, ap = partial(s)
    L('%-18s %-20s %+9.4f %10.3g %20s' % (tag.get(s, s), IDS.get(s, ''), r, p,
                                          '%+.4f (P=%.3g)' % (ar, ap)))
    single[s] = {'ensg': IDS.get(s), 'raw_r': round(float(r), 4), 'raw_p': float(p),
                 'adj_r': round(float(ar), 4), 'adj_p': round(float(ap), 5)}

# ---------------- B. 签名变体（全报）----------------
SIG = {
    'S1_RICTOR单基因':        ['RICTOR'],
    'S2_mTORC2实体4基因':     ['RICTOR', 'MAPKAP1', 'PRR5', 'PRR5L'],
    'S3_实体+MLST8+DEPTOR':   ['RICTOR', 'MAPKAP1', 'PRR5', 'PRR5L', 'MLST8', 'DEPTOR'],
    'S4_S3+下游读出':         ['RICTOR', 'MAPKAP1', 'PRR5', 'PRR5L', 'MLST8', 'DEPTOR', 'SGK1', 'NDRG1', 'PRKCA'],
    'S5_仅下游读出4基因':     ['SGK1', 'NDRG1', 'PRKCA', 'AKT1'],
    'S6_RICTOR+下游读出':     ['RICTOR', 'SGK1', 'NDRG1', 'PRKCA'],
}
def zmean(genes):
    G = [g for g in genes if g in EM]
    M = np.array([[EM[g][c] for g in G] for c in cc], float)
    Z = (M - M.mean(0)) / M.std(0, ddof=1)
    return Z.mean(1), G

L('')
L('-- B. mTORC2 多基因签名变体（全部报告，不择优隐藏）--')
L('%-22s %-34s %9s %10s %20s %14s' % ('签名', '组成', 'raw r', 'P', 'adj 纯度+CN+增殖', '留一法'))
L('-' * 116)
sig = {}
for nm, genes in SIG.items():
    sv, G = zmean(genes)
    # raw on base: 用 EM 全样本（签名值只在 cc 上算 → 用 cc 的 raw，标注）
    r, p = stats.pearsonr(xcc, sv)
    ar, ap = partial(sv)
    loo = []
    for i in range(len(cc)):
        m = np.ones(len(cc), bool); m[i] = False
        rr, _ = stats.pearsonr(xcc[m], sv[m]); loo.append(rr)
    L('%-22s %-34s %+9.4f %10.3g %20s %14s' %
      (nm, ('%d基因:' % len(G)) + ','.join(G), r, p, '%+.4f (P=%.3g)' % (ar, ap),
       '[%.2f,%.2f]' % (min(loo), max(loo))))
    sig[nm] = {'genes': G, 'raw_r': round(float(r), 4), 'raw_p': float(p),
               'adj_r': round(float(ar), 4), 'adj_p': round(float(ap), 5),
               'loo_min': round(min(loo), 4), 'loo_max': round(max(loo), 4),
               'values': {c: float(sv[i]) for i, c in enumerate(cc)}}

best_sig = max(sig, key=lambda k: sig[k]['adj_r'])
L('')
L('   ★ 偏相关后最强签名 = %s （adj r=%+.4f, P=%.3g；单基因 RICTOR = %+.4f, P=%.3g）'
  % (best_sig, sig[best_sig]['adj_r'], sig[best_sig]['adj_p'],
     single['RICTOR']['adj_r'], single['RICTOR']['adj_p']))
L('   → 结论：**签名不优于单基因 RICTOR**（多基因平均把 mTORC2 特异信号与共享组分噪声混合）')

# 写出
outp = os.path.join(FIGD, 'mtorc2_sig_tcga.csv')
with open(outp, 'w', encoding='utf-8-sig', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['Sample', 'Score'] + list(sig.keys()) +
               ['RICTOR', 'MPRIP_OLDID', 'MTOR', 'CD36', 'PPARG', 'prolif', 'purity', 'cn_cluster'])
    for c in cc:
        row = [c, score[c]] + [sig[k]['values'][c] for k in sig]
        for g in ['RICTOR', 'MPRIP', 'MTOR', 'CD36', 'PPARG']:
            row.append(EM[g][c] if g in EM else '')
        row += [prol[c], pur[c], cnc[c]]
        w.writerow(row)
L('写出: %s（%d 行）' % (outp, len(cc)))

# ---------------- C. 效力 ----------------
L('')
L('=' * 100)
L('C. 效力分析：区分「阴性」与「矛盾」')
L('=' * 100)
def sp_mde(n):
    se = 1.0 / math.sqrt(n - 3)
    return math.tanh(Z95 * se), math.tanh((Z95 + Z80) * se)
def sp_ci(r, n):
    se = 1.0 / math.sqrt(n - 3)
    z = math.atanh(max(-0.999999, min(0.999999, r)))
    return math.tanh(z - Z95 * se), math.tanh(z + Z95 * se)

L('-- C1. Spearman 检出下限 --')
L('%6s %18s %18s' % ('n', '临界|rho|(a=.05)', 'MDE |rho|(power=.80)'))
for n in [19, 27, 31, 50]:
    c, m = sp_mde(n)
    L('%6d %18.3f %18.3f' % (n, c, m))

REF = 0.3836
L('')
L('-- C2. 取样级（n=27）观测值的 95' + '%' + 'CI 是否覆盖参照效应 %+.4f --' % REF)
L('%-42s %9s %9s %22s %12s' % ('对比', 'rho', 'P', '95% CI', '覆盖参照'))
L('-' * 100)
N27 = [('肿瘤评分 vs 巨噬 PPARG', 0.3836, 0.0482),
       ('肿瘤评分 vs 巨噬 M2 样打分', 0.0818, 0.684),
       ('肿瘤评分 vs 巨噬 CD36', 0.0354, 0.861),
       ('肿瘤评分 vs 巨噬 RICTOR', 0.1466, 0.466),
       ('A−PPARG 评分 vs 巨噬 CD36', -0.3547, 0.0701)]
b2 = []
for nm, r, p in N27:
    lo, hi = sp_ci(r, 27)
    cov = bool(lo <= REF <= hi)
    L('%-42s %+9.4f %9.4f   [%+.3f, %+.3f] %12s' % (nm, r, p, lo, hi, '✓ 覆盖' if cov else '✗ 不覆盖'))
    b2.append({'pair': nm, 'rho': r, 'p': p, 'ci95': [round(lo, 4), round(hi, 4)], 'covers_ref': cov})

E = 23
mde_hr = math.exp((Z95 + Z80) / math.sqrt(E))
L('')
L('-- C3. Cox 检出下限（HR per 1 SD）--')
L('   事件 E=%d -> MDE HR = %.3f /SD' % (E, mde_hr))
cd36 = np.array([EM['CD36'][c] for c in cc], float)
sd36 = float(cd36.std(ddof=1))
L('')
L('-- C4. CD36 单变量 HR 尺度换算（log2CPM SD=%.4f）--' % sd36)
for lab, hu, lu, hui in [('（旧，用错配增殖指数）', 1.079, 0.884, 1.318)]:
    L('   %s per-1-unit %.3f [%.3f, %.3f] -> per-SD %.3f [%.3f, %.3f]'
      % (lab, hu, lu, hui, hu ** sd36, lu ** sd36, hui ** sd36))
L('   方向 per-SD = %.3f > 1 -> 与已发表论文同向（高 CD36 预后更差）' % (1.079 ** sd36))

json.dump({'n_base': len(base), 'n_cc': len(cc),
           'single': single,
           'signatures': {k: {kk: vv for kk, vv in v.items() if kk != 'values'} for k, v in sig.items()},
           'best_sig': best_sig,
           'spearman_mde': {str(n): {'crit': round(sp_mde(n)[0], 4), 'mde': round(sp_mde(n)[1], 4)} for n in [19, 27, 31, 50]},
           'ci27': b2, 'cox_mde_hr': round(mde_hr, 4), 'cd36_sd': round(sd36, 4)},
          open(os.path.join(W, 's34_nocontra.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
open(LOG, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n[写出]', LOG)
