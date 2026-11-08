# -*- coding: utf-8 -*-
"""手稿评分 C 的分量一致性：FASN/SREBF1 的「负贡献」问题在 C 上是否也存在？
同时给出 C 与 6 个分量基因、以及脂代谢相关基因的相关，用于判断 C 的构造偏向。
"""
import os, csv, json
import numpy as np, pyreadr
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

C = {}
for r in csv.DictReader(open(os.path.join(RES, '02_feature_scores.csv'), encoding='utf-8-sig')):
    k = r.get('') or r.get('Sample')
    try: C[str(k).strip()] = float(r['Lipid_Metabolism_Score'])
    except Exception: pass
A = {r['Sample']: float(r['Score']) for r in
     csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}

big = pyreadr.read_r(os.path.join(RAW, 'TCGA_LPS_real_log2cpm.rds'))
big = big[list(big.keys())[0]]
bidx = {str(k).split('.')[0]: str(k) for k in big.index}
ENS = {'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
       'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710',
       'CD36': 'ENSG00000135218', 'FABP4': 'ENSG00000170323', 'LPL': 'ENSG00000175445',
       'ADIPOQ': 'ENSG00000181092', 'PLIN1': 'ENSG00000166819', 'CPT1A': 'ENSG00000110090',
       'CD8A': 'ENSG00000153563'}

samples = sorted(set(A) & set(C) & set(str(c) for c in big.columns))
rep['n'] = len(samples)
vals = {}
for g, e in ENS.items():
    row = bidx.get(e)
    if not row: continue
    v = []
    for s in samples:
        try: v.append(float(big.loc[row, s]))
        except Exception: v.append(np.nan)
    vals[g] = np.array(v, float)

a = np.array([A[s] for s in samples]); c = np.array([C[s] for s in samples])
out = {}
for g, v in vals.items():
    ra, pa = stats.pearsonr(v, a); rc, pc = stats.pearsonr(v, c)
    out[g] = {'vs_A': {'r': round(float(ra), 4), 'p': round(float(pa), 6)},
              'vs_C': {'r': round(float(rc), 4), 'p': round(float(pc), 6)}}
rep['component_vs_score'] = out

# C 是否为「上调和/下调」方向：C 与 A 的符号一致率
rep['sign_agreement_A_C'] = round(float(np.mean(np.sign(a - a.mean()) == np.sign(c - c.mean()))), 4)
rep['corr_A_C'] = {'pearson': round(float(stats.pearsonr(a, c)[0]), 4),
                   'spearman': round(float(stats.spearmanr(a, c)[0]), 4)}

json.dump(rep, open(os.path.join(OUT, 'd15g_C_components.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=1))
