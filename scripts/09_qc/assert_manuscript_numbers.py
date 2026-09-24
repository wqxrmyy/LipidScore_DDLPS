#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GSE221492 主结果交叉校验 —— 两套并存实现 + 本文报告值 三重断言。

背景
----
`scripts/07_bulk_validation/` 下有**两套并存实现**：
  · 01_lipid_score_gse221492.py  (Python) —— 本文正文所报数字的原始来源
  · 01_lipid_score_bulk.R        (R)      —— 独立重实现
两套必须是同一口径。历史上 R 侧曾有两处口径错误，会把本文报告的第 3 位小数算错：
  [FIX-1] library size 算在「合并重复 hgnc_symbol 之后」-> 分母偏小 0.47%
  [FIX-2] DDLPS_sWD 被并入 DDLPS_WD -> 配对集从 8 对变 10 对

本脚本检查三层
--------------
  L1 两实现一致：逐患者 Δ、n、DD<WD 计数、P 值、mean/median Δ
  L2 本文报告值：n=8 / 6-of-8 / P=0.078 / mean Δ=-0.587 /
                 非配对 DD -0.305、WD +0.164、Δ -0.469、Welch P 0.025 / AUC 0.249、P 0.013
  L3 （可选 --manuscript）稿件正文确实印着这些数字 —— 防止「代码对、正文错」

用法
----
  python scripts/09_qc/assert_manuscript_numbers.py
  python scripts/09_qc/assert_manuscript_numbers.py --manuscript /path/to/Manuscript.md
  python scripts/09_qc/assert_manuscript_numbers.py --root /path/to/repo

退出码：0 = 全部通过；1 = 有失败项（逐条打印）。
"""
import argparse
import csv
import json
import os
import re
import sys

_ap = argparse.ArgumentParser(description='GSE221492 主结果交叉校验')
_ap.add_argument('--root', default=None, help='仓库根（默认由本脚本位置推断）')
_ap.add_argument('--manuscript', default=None, help='稿件正文（可选，做 L3 文字核对）')
_ap.add_argument('--tol-delta', type=float, default=1e-4,
                 help='两实现 Δ 的绝对容差（Python 侧 JSON 只保留 4 位小数）')
_a = _ap.parse_args()

ROOT = _a.root or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
TAB = os.path.join(ROOT, 'results', 'tables')
OBJ = os.path.join(ROOT, 'results', 'objects')

fails = []
warns = []


def check(cond, msg):
    print(('  ✅ ' if cond else '  ❌ ') + msg)
    if not cond:
        fails.append(msg)


def warn(cond, msg):
    if not cond:
        print('  ⚠️  ' + msg)
        warns.append(msg)


def read_csv(p):
    with open(p, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


# ---------- 载入两实现 ----------
py_sum = os.path.join(TAB, 'bulk_lipid_score_py_summary.csv')
py_del = os.path.join(TAB, 'bulk_lipid_score_py_paired_delta.csv')
r_sum = os.path.join(TAB, 'bulk_lipid_score_summary.csv')
r_del = os.path.join(TAB, 'bulk_lipid_score_paired_delta.csv')
py_json = os.path.join(OBJ, 'd11_gse221492.json')

for p in (py_sum, py_del, r_sum, r_del, py_json):
    if not os.path.exists(p):
        print('❌ 缺少输入: %s\n   请先跑两套实现：\n'
              '     python scripts/07_bulk_validation/01_lipid_score_gse221492.py\n'
              '     Rscript scripts/07_bulk_validation/01_lipid_score_bulk.R' % p)
        sys.exit(1)

J = json.load(open(py_json, encoding='utf-8'))
B5 = J['B5_selfcheck']
PT = B5['paired_test_all53']
CMP = B5['DD_vs_WD_DDLPS_specimens']['compare']
AUC = B5['DD_vs_WD_DDLPS_specimens']['auc']

py_s = {r['analysis']: r for r in read_csv(py_sum)}
r_s = {r['analysis']: r for r in read_csv(r_sum)}
py_d = {re.sub(r'^P', '', r['patient']): float(r['delta']) for r in read_csv(py_del)}
r_d = {re.sub(r'^P', '', r['patient']): float(r['delta']) for r in read_csv(r_del)}

PY_MAIN = 'Main (strict pairs)'
R_MAIN = next((k for k in r_s if 'strict' in k.lower()), None)
if R_MAIN is None:
    print('❌ R summary 表里找不到主分析行，实际行: %s' % list(r_s))
    sys.exit(1)

print('=' * 78)
print('L1 · 两实现一致性')
print('=' * 78)
check(set(py_d) == set(r_d),
      '配对患者集合一致: %d 人  %s' % (len(py_d), sorted(py_d, key=int)))
if set(py_d) == set(r_d):
    diffs = {k: abs(py_d[k] - r_d[k]) for k in py_d}
    mx = max(diffs.values())
    check(mx <= _a.tol_delta,
          '逐患者 Δ 一致（最大绝对差 %.2e <= 容差 %.0e；差异来自 Python 侧 JSON 仅存 4 位小数）'
          % (mx, _a.tol_delta))
    if diffs:
        worst = max(diffs, key=diffs.get)
        print('      最大差落在患者 P%s: py=%+.6f  R=%+.6f' % (worst, py_d[worst], r_d[worst]))

check(int(r_s[R_MAIN]['n_pairs']) == int(py_s[PY_MAIN]['n_pairs']),
      'n_pairs 一致: R=%s  Python=%s' % (r_s[R_MAIN]['n_pairs'], py_s[PY_MAIN]['n_pairs']))
check(int(r_s[R_MAIN]['n_dd_lt_wd']) == int(py_s[PY_MAIN]['n_dd_lt_wd']),
      'DD<WD 计数一致: R=%s  Python=%s' % (r_s[R_MAIN]['n_dd_lt_wd'],
                                           py_s[PY_MAIN]['n_dd_lt_wd']))
check(abs(float(r_s[R_MAIN]['p_value']) - float(py_s[PY_MAIN]['p_value'])) < 1e-9,
      'Wilcoxon P 一致: R=%s  Python=%s' % (r_s[R_MAIN]['p_value'], py_s[PY_MAIN]['p_value']))
check(abs(float(r_s[R_MAIN]['mean_delta']) - float(py_s[PY_MAIN]['mean_delta'])) <= _a.tol_delta,
      'mean Δ 一致: R=%.6f  Python=%.6f' % (float(r_s[R_MAIN]['mean_delta']),
                                            float(py_s[PY_MAIN]['mean_delta'])))
check(abs(float(r_s[R_MAIN]['median_delta']) - float(py_s[PY_MAIN]['median_delta'])) <= _a.tol_delta,
      'median Δ 一致: R=%.6f  Python=%.6f' % (float(r_s[R_MAIN]['median_delta']),
                                              float(py_s[PY_MAIN]['median_delta'])))

# 防呆：mean 与 median 必须被区分（历史上正是把 median 当成了报告值）
md = abs(float(py_s[PY_MAIN]['mean_delta']) - float(py_s[PY_MAIN]['median_delta']))
warn(md > 1e-6, 'mean Δ 与 median Δ 相差 %.4f —— 务必确认引用的是 mean（正文口径）' % md)
mcol = [k for k in py_s[PY_MAIN] if k.startswith('mean')]
warn(bool(mcol), '输出表缺少 mean_delta 列 —— 只给 median 会让「报告哪个值」重新变成隐患')

print()
print('=' * 78)
print('L2 · 本文报告值断言（3 位小数，两实现都查）')
print('=' * 78)
for tag, r in (('Python', py_s[PY_MAIN]), ('R', r_s[R_MAIN])):
    n = int(r['n_pairs'])
    lt = int(r['n_dd_lt_wd'])
    p = float(r['p_value'])
    mean = float(r['mean_delta'])
    check(n == 8, '%s: 严格配对数 n == 8（实得 %d）' % (tag, n))
    check(lt == 6, '%s: DD<WD == 6（正文 “six of eight”，实得 %d）' % (tag, lt))
    check(abs(p - 0.078125) < 1e-9, '%s: Wilcoxon P == 0.078125（实得 %.6f）' % (tag, p))
    check(round(mean, 3) == -0.587, '%s: mean Δ -> -0.587（实得 %.6f -> %+.3f）'
          % (tag, mean, round(mean, 3)))

check(round(CMP['mean1'], 3) == -0.305, '非配对 DD 均值 -> -0.305（实得 %.4f）' % CMP['mean1'])
check(round(CMP['mean0'], 3) == 0.164, '非配对 WD 均值 -> +0.164（实得 %.4f）' % CMP['mean0'])
check(round(CMP['delta'], 3) == -0.469, '非配对 Δ -> -0.469（实得 %.4f）' % CMP['delta'])
check(round(CMP['p'], 3) == 0.025, '非配对 Welch P -> 0.025（实得 %.6f）' % CMP['p'])
check(round(AUC['auc'], 3) == 0.249, 'AUC -> 0.249（实得 %.4f）' % AUC['auc'])
check(round(AUC['p'], 3) == 0.013, 'AUC P -> 0.013（实得 %.6f）' % AUC['p'])
check(B5['n_strict_pairs'] == 8, 'JSON n_strict_pairs == 8（实得 %d）' % B5['n_strict_pairs'])

# ---------- L3：稿件正文核对（可选） ----------
if _a.manuscript:
    print()
    print('=' * 78)
    print('L3 · 稿件正文核对: %s' % os.path.basename(_a.manuscript))
    print('=' * 78)
    txt = open(_a.manuscript, encoding='utf-8').read()
    for pat, desc in [
        (r'eight patients with strictly paired components', '正文写明 strict paired / 8 例'),
        (r'six showed a lower score', '正文写明 6 of 8 更低'),
        (r'0\.587', '正文含 Δ = -0.587'),
        (r'P\s*=\s*0\.078', '正文含 P = 0.078'),
        (r'0\.305', '正文含 DD 均值 -0.305'),
        (r'0\.164', '正文含 WD 均值 +0.164'),
        (r'0\.469', '正文含非配对 Δ -0.469'),
        (r'0\.249', '正文含 AUC 0.249'),
    ]:
        check(bool(re.search(pat, txt)), '正文: %s' % desc)
else:
    print()
    print('L3 · 稿件正文核对 …… 跳过（未提供 --manuscript）')

print()
print('=' * 78)
print('交叉校验: %s' % ('✅ 全部通过' if not fails else '❌ %d 项失败' % len(fails)))
for f in fails:
    print('   -', f)
if warns:
    print('提示 %d 条:' % len(warns))
    for w in warns:
        print('   ·', w)
sys.exit(0 if not fails else 1)
