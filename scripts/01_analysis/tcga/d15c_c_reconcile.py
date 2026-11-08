# -*- coding: utf-8 -*-
"""核对：S4 用的 Lipid_Metabolism_Score 是否等于 TCGA_Score.csv 的 Score"""
import os, csv, json
import numpy as np
from scipy import stats

RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")
RES = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
rep = {}

fs = os.path.join(RES, '02_feature_scores.csv')
rep['feature_scores_exists'] = os.path.exists(fs)
rows = list(csv.DictReader(open(fs, encoding='utf-8-sig')))
rep['feature_scores_nrow'] = len(rows)
rep['feature_scores_cols'] = list(rows[0].keys()) if rows else []
rep['feature_scores_index_col'] = list(rows[0].keys())[0] if rows else None
key = list(rows[0].keys())[0]

sc = {r['Sample']: float(r['Score']) for r in
      csv.DictReader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig'))}

pairs = []
for r in rows:
    k = r[key]
    if k in sc and r.get('Lipid_Metabolism_Score') not in (None, '', 'NA'):
        try:
            pairs.append((k, float(r['Lipid_Metabolism_Score']), sc[k]))
        except Exception:
            pass
rep['n_matched'] = len(pairs)
if len(pairs) >= 5:
    a = np.array([p[1] for p in pairs]); b = np.array([p[2] for p in pairs])
    r_, p_ = stats.pearsonr(a, b)
    rep['compare_LMS_vs_TCGAScore'] = {
        'pearson': round(float(r_), 6), 'p': float(p_),
        'max_abs_diff': round(float(np.max(np.abs(a - b))), 6),
        'mean_LMS': round(float(a.mean()), 4), 'mean_Score': round(float(b.mean()), 4),
        'sd_LMS': round(float(a.std(ddof=1)), 4), 'sd_Score': round(float(b.std(ddof=1)), 4),
        'n_identical': int(np.sum(np.abs(a - b) < 1e-9))}
    # 是否只是线性变换
    slope, inter = np.polyfit(b, a, 1)
    rep['linear_fit_LMS_on_Score'] = {'slope': round(float(slope), 5), 'intercept': round(float(inter), 5),
                                      'r2': round(float(r_ ** 2), 6)}

# 用 LMS 列算 Harrell C，看是否复现 0.6175
cl = {r['barcode']: r for r in
      csv.DictReader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig'))}


def num(x):
    try:
        v = float(str(x).strip()); return None if np.isnan(v) else v
    except Exception:
        return None


try:
    from lifelines.utils import concordance_index
    def cval(getter):
        T, E_, S_ = [], [], []
        for k in sc:
            r0 = cl.get(k)
            if not r0: continue
            t = num(r0.get('OS_time')) or num(r0.get('paper_OS.days'))
            e = num(r0.get('OS_status'))
            v = getter(k)
            if t and e is not None and v is not None:
                T.append(t); E_.append(int(e)); S_.append(v)
        if len(T) < 10: return None
        return {'n': len(T), 'events': int(sum(E_)),
                'harrell_C': round(float(concordance_index(T, [-x for x in S_], E_)), 4)}
    rep['C_using_TCGA_Score_csv'] = cval(lambda k: sc.get(k))
    lms = {p[0]: p[1] for p in pairs}
    rep['C_using_Lipid_Metabolism_Score'] = cval(lambda k: lms.get(k))
except Exception as e:
    rep['c_error'] = repr(e)[:200]

json.dump(rep, open(os.path.join(OUT, 'd15c_c_reconcile.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
print(json.dumps(rep, ensure_ascii=False, indent=2))
