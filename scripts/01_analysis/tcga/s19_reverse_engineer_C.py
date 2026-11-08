import os
# -*- coding: utf-8 -*-
"""
S19 · 反推并冻结 C 评分（02_feature_scores.csv 的 Lipid_Metabolism_Score）的官方定义
D16 已裁定：C 生成脚本取不回，授权 P4 以现有文件反推冻结。
判据：反推出的公式必须在 253 行上精确复现 C（Pearson r = 1.000000 或残差≈0）。
"""
import pyreadr, csv, math, statistics, itertools, json

base = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s19_reverse_engineer_C.json")

# ---- 读 C ----
c_map = {}
with open(base + 'results/02_feature_scores.csv', encoding='utf-8-sig') as f:
    rows = list(csv.reader(f))
for r in rows[1:]:
    c_map[r[0].strip()] = float(r[1])

# ---- 读三个队列 6 基因表达（log2 尺度） ----
expr = {}
for q in ['TCGA', 'GSE21122', 'GSE159659']:
    df = list(pyreadr.read_r(base + 'raw/' + q + '_6genes_expr.rds').values())[0]
    for sid in df.columns:
        sidc = sid[:-7] if sid.endswith('.CEL.gz') else sid
        expr[sidc] = {g: float(df.at[g, sid]) for g in df.index}

genes = ['SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN']
common = sorted(set(c_map) & set(expr))
print('样本数 expr=%d  C=%d  匹配=%d' % (len(expr), len(c_map), len(common)))

results = {}

def score(vals):
    n = len(common)
    ys = [c_map[s] for s in common]
    mx = sum(vals) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in vals)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(vals, ys))
    r = sxy / math.sqrt(sxx * syy)
    slope = sxy / sxx
    icpt = my - slope * mx
    resid = [y - (slope * x + icpt) for x, y in zip(vals, ys)]
    return {'r': r, 'slope': slope, 'icpt': icpt, 'resid_max': max(abs(v) for v in resid)}

def try_candidate(name, vals):
    results[name] = score(vals)

# ---- 候选公式（全部在 6 基因 log2 表达上构造） ----
# 1. 6 基因 log2 均值
try_candidate('mean_log2_6gene', [sum(expr[s][g] for g in genes)/6 for s in common])
# 2. 6 基因 log2 中位数
try_candidate('median_log2_6gene', [statistics.median([expr[s][g] for g in genes]) for s in common])
# 3. 去 SREBF1 的 5 基因均值
try_candidate('mean_log2_5gene_noSREBF1', [sum(expr[s][g] for g in genes if g != 'SREBF1')/5 for s in common])
# 4. PPARG 单独
try_candidate('PPARG_only', [expr[s]['PPARG'] for s in common])
# 5. 每队列内 6 基因 z 分后均值（=A）
def zscore_mean():
    # 每队列分别 z 分
    per_q = {}
    for q in ['TCGA', 'GSE21122', 'GSE159659']:
        df = list(pyreadr.read_r(base + 'raw/' + q + '_6genes_expr.rds').values())[0]
        # 逐基因跨样本 z 分
        import numpy as np
        zdf = df.copy()
        for g in df.index:
            v = np.array([float(x) for x in df.loc[g]])
            zdf.loc[g] = (v - v.mean()) / v.std()
        for sid in df.columns:
            sidc = sid[:-7] if sid.endswith('.CEL.gz') else sid
            per_q[sidc] = sum(zdf.at[g, sid] for g in df.index) / 6
    return [per_q[s] for s in common]
try_candidate('A_zscore_mean_6gene', zscore_mean())

# 6. 反 log2 后均值（线性空间）
try_candidate('mean_linear_6gene', [sum(2 ** expr[s][g] for g in genes)/6 for s in common])

# 7. ssGSEA 简化：秩和（每个样本内 6 基因按表达排序，取某种统计）
#    （6 基因太少，ssGSEA 意义有限，跳过）

# ---- 打印结果，按 r 排序 ----
print('\n=== 候选公式与 C 的拟合 ===')
for name, sc in sorted(results.items(), key=lambda kv: -abs(kv[1]['r'])):
    print('%-28s  r=%.8f  resid_max=%.6f  slope=%.4f  icpt=%.4f' % (
        name, sc['r'], sc['resid_max'], sc['slope'], sc['icpt']))

json.dump(results, open(OUT, 'w'), indent=2, ensure_ascii=False)
print('\n已写入', OUT)
