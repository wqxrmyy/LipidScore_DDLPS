# -*- coding: utf-8 -*-
"""D7 第二步：核查 MDM2 状态的三种可能来源，并定位手稿 R=-0.102/P=0.48 的真实出处。
   来源1: TCGA CNA 列（已验证为近常量）
   来源2: TCGA 表达矩阵中的 MDM2 表达（本脚本重点）
   来源3: 12q13-15 共扩增伙伴基因
"""
import os, json, math
import numpy as np
import pyreadr

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

def pear(x, y):
    from scipy import stats
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5: return None
    r, p = stats.pearsonr(x, y)
    return [round(float(r), 4), round(float(p), 5), int(len(x))]

def spear(x, y):
    from scipy import stats
    if len(x) < 5: return None
    rs, p = stats.spearmanr(x, y)
    return [round(float(rs), 4), round(float(p), 5), int(len(x))]

# ---- score ----
import csv
srows = list(csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig')))
score = {r['Sample']: float(r['Score']) for r in srows}

def k3(b):
    return '-'.join(str(b).split('-')[:3])

# ---- 表达矩阵 ----
for fn in ['TCGA_LPS_real_log2cpm.rds', 'TCGA_LPS_real_counts_filtered.rds', 'TCGA_6genes_expr.rds']:
    p = os.path.join(RAW, fn)
    if not os.path.exists(p): continue
    try:
        res = pyreadr.read_r(p)
    except Exception as e:
        rep.setdefault('load_errors', {})[fn] = str(e); continue
    name = list(res.keys())[0]
    df = res[name]
    info = {'key': name, 'shape': list(df.shape), 'index_head': [str(x) for x in list(df.index[:4])],
            'cols_head': [str(c) for c in list(df.columns[:4])]}
    genes = [str(g) for g in df.index]
    cols = [str(c) for c in df.columns]
    # 判断方向：基因在 index 还是 column
    orient = 'genes_as_index' if any('MDM2' == g or g.upper().startswith('MDM2') for g in genes) else 'genes_as_columns'
    info['orient'] = orient
    info['n_genes'] = df.shape[0] if orient == 'genes_as_index' else df.shape[1]
    info['has_MDM2'] = ('MDM2' in genes) or ('MDM2' in cols)
    rep.setdefault('matrices', {})[fn] = info

    if orient == 'genes_as_index':
        sub = df.loc[[g for g in genes if g in ('MDM2', 'CDK4', 'HMGA2', 'DDIT3', 'FRS2', 'PPARG', 'CD36', 'FASN', 'SCD', 'SREBF1', 'PLIN1', 'ADIPOQ')]]
        sample_map = {}
        for c in cols:
            smp = score.get(c)
            if smp is None:
                smp = score.get(k3(c))
            if smp is None:
                for s in score:
                    if c.startswith(s) or s.startswith(c):
                        smp = score[s]; break
            if smp is not None:
                sample_map[c] = smp
        info['n_samples_matched'] = len(sample_map)
        out = {}
        for g in list(sub.index):
            x = []
            y = []
            for c, s in sample_map.items():
                try:
                    v = float(sub.loc[g, c])
                except Exception:
                    continue
                if math.isnan(v): continue
                x.append(v); y.append(s)
            out[g] = {'pearson': pear(x, y), 'spearman': spear(x, y),
                      'range': [round(min(x), 3), round(max(x), 3)] if x else None}
        rep.setdefault('score_vs_expression', {})[fn] = out

json.dump(rep, open(os.path.join(OUT, 'd7_expression.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
