# -*- coding: utf-8 -*-
"""细胞谱系构成梳理：L1/L2/L3 全谱 + 组织分布 + 7 对配对的 WD vs DD 组成对比。
用途：为 Step 5（细胞类型定义）与 Step 7（Kitagawa 分解选层）提供事实基础。"""
import os, gzip, json
from collections import defaultdict, Counter

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
G = os.path.join(W, 'gse221493')
ANN = os.path.join(G, 'GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz')
OUT = os.path.join(W, 's39_s3c_lineage.txt')

rows = []
with gzip.open(ANN, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 10:
            continue
        rows.append((p[0], p[1], p[2], p[4], p[5], p[6], p[7]))   # bc, pat, tissue, grade, L1, L2, L3

log = []
def A(s):
    log.append(s); print(s, flush=True)

n = len(rows)
A('细胞总数 %d' % n)

ORDER = ['nl', 'lipoma', 'WDLPS', 'DDLPS-WD', 'DDLPS-sWD', 'DDLPS-DD', 'DDLPS-DDK']
tissues = [t for t in ORDER if any(r[2] == t for r in rows)]
tissues += [t for t in sorted({r[2] for r in rows}) if t not in tissues]

# ---------- 1. L1 全谱 × 组织 ----------
A('')
A('=' * 100)
A('1. L1 谱系大类 × 组织（细胞数 / 该组织内行百分比 %）')
A('=' * 100)
l1s = sorted({r[4] for r in rows}, key=lambda x: -sum(1 for r in rows if r[4] == x))
hdr = '%-24s %8s' % ('L1 lineage', 'total')
for t in tissues:
    hdr += ' %12s' % t
A(hdr)
for l1 in l1s:
    tot = sum(1 for r in rows if r[4] == l1)
    line = '%-24s %8d' % (l1, tot)
    for t in tissues:
        c = sum(1 for r in rows if r[4] == l1 and r[2] == t)
        nt = sum(1 for r in rows if r[2] == t)
        line += ' %7d(%4.1f)' % (c, 100.0 * c / nt if nt else 0)
    A(line)
line = '%-24s %8d' % ('-- total --', n)
for t in tissues:
    line += ' %12d' % sum(1 for r in rows if r[2] == t)
A(line)

# ---------- 2. L2 类别 ----------
A('')
A('=' * 100)
A('2. L2 中间类（%d 类）' % len({r[5] for r in rows}))
A('=' * 100)
l2c = Counter(r[5] for r in rows)
for k, v in l2c.most_common():
    l1 = next(r[4] for r in rows if r[5] == k)
    A('  %-34s %7d   [%s]' % (k, v, l1))

# ---------- 3. L3 类别 ----------
A('')
A('=' * 100)
A('3. L3 精细类（%d 类）' % len({r[6] for r in rows}))
A('=' * 100)
l3c = Counter(r[6] for r in rows)
for k, v in l3c.most_common():
    l2 = next(r[5] for r in rows if r[6] == k)
    l1 = next(r[4] for r in rows if r[6] == k)
    A('  %-34s %7d   [%s > %s]' % (k, v, l1, l2))

# ---------- 4. 7 对配对：L1 组成 WD vs DD ----------
PAIRS = ['P1', 'P4', 'P7', 'P10', 'P14', 'P15', 'P16']
A('')
A('=' * 100)
A('4. ★ 严格配对 7 例：L1 组成 WD vs DD（组织内 %）')
A('=' * 100)
A('%-24s' % 'L1 lineage' + ''.join(' %14s' % p for p in PAIRS))
comp = {}
for l1 in l1s:
    line = '%-24s' % l1
    for p in PAIRS:
        nd = sum(1 for r in rows if r[1] == p and r[2] == 'DDLPS-DD')
        nw = sum(1 for r in rows if r[1] == p and r[2] == 'DDLPS-WD')
        cd = sum(1 for r in rows if r[1] == p and r[2] == 'DDLPS-DD' and r[4] == l1)
        cw = sum(1 for r in rows if r[1] == p and r[2] == 'DDLPS-WD' and r[4] == l1)
        fd = 100.0 * cd / nd if nd else float('nan')
        fw = 100.0 * cw / nw if nw else float('nan')
        comp[(p, l1)] = (fw, fd)
        line += ' %6.1f>%5.1f' % (fw, fd)
    A(line)

A('')
A('  （格式 WD>DD；粗看哪些谱系在去分化时扩张/收缩）')
A('')
A('  逐谱系方向汇总（WD→DD）：')
for l1 in l1s:
    ds = [comp[(p, l1)][1] - comp[(p, l1)][0] for p in PAIRS]
    up = sum(1 for d in ds if d > 0)
    A('    %-24s DD 高于 WD 的患者数 %d/7   Δ中位 %+6.2f 个百分点' % (l1, up, sorted(ds)[3]))

# ---------- 5. 脂肪谱系相关类的去向 ----------
A('')
A('=' * 100)
A('5. 脂肪谱系 / 肿瘤类相关的 L2、L3 分布')
A('=' * 100)
KEY_L2 = [k for k in l2c if any(s in k for s in ('Adipo', 'Tumor', 'Stem', 'Invasion', 'ECM', 'Hypox', 'Angio', 'Cycling'))]
for k in sorted(KEY_L2, key=lambda x: -l2c[x]):
    A('')
    A('  [L2] %s  (n=%d)' % (k, l2c[k]))
    A('       %s' % ''.join(' %14s' % t for t in tissues))
    line = '       %14s' % 'n'
    for t in tissues:
        c = sum(1 for r in rows if r[5] == k and r[2] == t)
        nt = sum(1 for r in rows if r[2] == t)
        line += ' %7d(%4.1f)' % (c, 100.0 * c / nt if nt else 0)
    A(line)

txt = '\n'.join(log)
open(OUT, 'w', encoding='utf-8').write(txt)
print('\nwritten ->', OUT)
