# -*- coding: utf-8 -*-
"""D7 可行性评估：TCGA-SARC MDM2 拷贝数状态是否可用于检验假说 A。
一手数据：TCGA_clinical_real.csv（TCGA-SARC 论文自带 CNA 列）+ TCGA_Score.csv
"""
import os, csv, json, math
import numpy as np

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

# ---------- 1. 读 Score ----------
sp = os.path.join(RAW, 'TCGA_Score.csv')
rows = list(csv.DictReader(open(sp, encoding='utf-8-sig', errors='ignore')))
rep['score_file'] = sp
rep['score_nrow'] = len(rows)
rep['score_cols'] = list(rows[0].keys())
score_col = None
for c in rows[0].keys():
    if 'score' in c.lower() or 'lipid' in c.lower():
        score_col = c
rep['score_col'] = score_col
scores = {}
for r in rows:
    b = r.get('Sample') or r.get('sample') or r.get('barcode') or r.get('Sample_short3')
    if b is None:
        b = list(r.values())[0]
    try:
        scores[str(b)] = float(r[score_col])
    except Exception:
        pass
rep['score_parsed'] = len(scores)
rep['score_head'] = list(scores.items())[:5]

# ---------- 2. 读临床 / CNA ----------
cl = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig', errors='ignore')))
cn_cols = [c for c in cl[0].keys() if c.startswith('paper_AMP.') or c.startswith('paper_avg.CN') or c == 'paper_CN.cluster' or c == 'paper_purity' or c == 'paper_ploidy']
rep['cn_cols'] = cn_cols
rep['clinical_nrow'] = len(cl)

def num(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except Exception:
        return None

# ---------- 3. MDM2 CN 分布（范围限制检验） ----------
dist = {}
for c in cn_cols:
    vals = [num(r.get(c)) for r in cl]
    vals = [v for v in vals if v is not None]
    if not vals:
        raw = [str(r.get(c)).strip() for r in cl if str(r.get(c)).strip() not in ('', 'NA', 'nan')]
        cats = {}
        for x in raw:
            cats[x] = cats.get(x, 0) + 1
        dist[c] = {'type': 'categorical', 'n_obs': len(raw), 'categories': dict(sorted(cats.items(), key=lambda kv: -kv[1])[:12])}
        continue
    a = np.array(vals, dtype=float)
    dist[c] = {
        'type': 'numeric', 'n_obs': len(a),
        'min': round(float(a.min()), 4), 'q25': round(float(np.percentile(a, 25)), 4),
        'median': round(float(np.median(a)), 4), 'q75': round(float(np.percentile(a, 75)), 4),
        'max': round(float(a.max()), 4), 'sd': round(float(a.std(ddof=1)), 4),
        'iqr': round(float(np.percentile(a, 75) - np.percentile(a, 25)), 4),
        # 常见扩增阈值（GISTIC log2 ratio 惯例：>0.2 / >0.3 / >1 分别对应 broad gain / gain / amplification）
        'n_gt_0.2': int((a > 0.2).sum()), 'n_gt_0.3': int((a > 0.3).sum()),
        'n_gt_0.5': int((a > 0.5).sum()), 'n_gt_1.0': int((a > 1.0).sum()),
        'n_le_0.2': int((a <= 0.2).sum()),
    }
rep['cn_distribution'] = dist

# ---------- 4. 关联检验：Score ~ CN ----------
def pearson(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 4: return None
    r = float(np.corrcoef(x, y)[0, 1])
    n = len(x)
    t = r * math.sqrt((n - 2) / max(1e-12, (1 - r * r)))
    from scipy import stats
    p = float(2 * stats.t.sf(abs(t), n - 2))
    return round(r, 4), round(p, 5), n

def spearman(x, y):
    from scipy import stats
    if len(x) < 4: return None
    rs, p = stats.spearmanr(x, y)
    return round(float(rs), 4), round(float(p), 5), len(x)

# 样本键匹配
def key_of(b):
    b = str(b)
    if '-' in b:
        parts = b.split('-')
        return '-'.join(parts[:3])
    return b

cli_by_key = {}
for r in cl:
    k = key_of(r.get('barcode', ''))
    cli_by_key[k] = r

assoc = {}
pairs = {}
for c in cn_cols:
    xs, ys = [], []
    for b, s in scores.items():
        r = cli_by_key.get(key_of(b))
        if not r: continue
        v = num(r.get(c))
        if v is None: continue
        xs.append(v); ys.append(s)
    if len(xs) >= 6:
        assoc[c] = {'pearson': pearson(xs, ys), 'spearman': spearman(xs, ys)}
        pairs[c] = (xs, ys)
rep['score_vs_cn'] = assoc

# ---------- 5. 分层：仅在 DDLPS 内看 MDM2 剂量 ----------
histo_col = 'paper_histology'
histo_counts = {}
for r in cl:
    histo_counts[str(r.get(histo_col)).strip()] = histo_counts.get(str(r.get(histo_col)).strip(), 0) + 1
rep['histology_counts'] = histo_counts

mdm2 = 'paper_AMP.12q15.MDM2'
strat = {}
if mdm2 in pairs:
    xs, ys = pairs[mdm2]
    # 按组织学分层
    for lab in ['DDLPS', 'WDLS', 'WDLPS', 'DDLPS ']:
        xs2, ys2 = [], []
        for b, s in scores.items():
            r = cli_by_key.get(key_of(b))
            if not r: continue
            if str(r.get(histo_col)).strip() != lab.strip(): continue
            v = num(r.get(mdm2))
            if v is None: continue
            xs2.append(v); ys2.append(s)
        if len(xs2) >= 6:
            strat[lab.strip()] = {'n': len(xs2), 'pearson': pearson(xs2, ys2), 'cn_range': [round(min(xs2), 3), round(max(xs2), 3)]}
rep['mdm2_within_histology'] = strat

# ---------- 6. 判定：是否具备「分层/交互」检验所需的对照分布 ----------
if mdm2 in dist and dist[mdm2].get('type') == 'numeric':
    d = dist[mdm2]
    rep['verdict_inputs'] = {
        'n_total': d['n_obs'],
        'n_amplified_gt0.2': d['n_gt_0.2'],
        'n_neutral_le0.2': d['n_le_0.2'],
        'sd': d['sd'], 'iqr': d['iqr'],
        'pct_amplified_gt0.2': round(100 * d['n_gt_0.2'] / d['n_obs'], 1),
        'pct_amplified_gt0.5': round(100 * d['n_gt_0.5'] / d['n_obs'], 1),
    }

# ---------- 7. 其他 12q13-15 成员是否与评分相关（用于新假说挖掘） ----------
rep['amplicon_scan'] = {c: assoc.get(c) for c in cn_cols if c.startswith('paper_AMP.')}

os.makedirs(OUT, exist_ok=True)
json.dump(rep, open(os.path.join(OUT, 'd7_mdm2_feasibility.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2)[:6000])
