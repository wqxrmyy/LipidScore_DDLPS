# -*- coding: utf-8 -*-
"""
S34 · 基因 ID 权威表（symbol -> ENSG），Ensembl REST symbol 反查生成
------------------------------------------------------------------
背景：S34 首次核验发现项目口径锁中 6 个 ENSG 与 Ensembl 不符，其中
      **RICTOR = ENSG00000133030 实为 MPRIP（chr17）**，RICTOR 正确为 ENSG00000164327（chr5）。
      性质同当年 FASN 探针错配（S14 撤回）。
做法：一律走 symbol -> ENSG 反查（权威方向），不再用「拟用 ID 去 lookup」。
输出：_work/s34_gene_ids.tsv（单一口径源，供后续脚本读取，杜绝二次硬编码）
"""
import json, os, io, urllib.request, time

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")

SYMS = [
    # 评分六基因（D18 锁 A）
    'SREBF1', 'SCD', 'HMGCR', 'ACLY', 'PPARG', 'FASN',
    # mTORC2 轴
    'RICTOR', 'MTOR', 'MAPKAP1', 'MLST8', 'PRR5', 'PRR5L', 'DEPTOR',
    'AKT1', 'SGK1', 'NDRG1', 'PRKCA', 'RPTOR',
    # 已发表论文相关
    'CD36', 'FABP4', 'MRC1', 'IL10', 'CD163', 'LPAR1', 'TGFB1', 'TGFBR1', 'TGFBR2',
    # 混杂对照
    'MKI67', 'PCNA', 'TOP2A', 'MCM2', 'CCNB1', 'ADIPOQ', 'PLIN1', 'MDM2',
    # 误配对照（本轮查出的错配目标，用于显式展示「此前实际测的是谁」）
    'MPRIP', 'RPA2', 'MINDY3', 'PRG2', 'KCTD18', 'MIER1',
]

# 此前项目口径锁 / 我第一版拟用（用于对照）
PREV = {
    'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
    'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710',
    'RICTOR': 'ENSG00000133030', 'MTOR': 'ENSG00000198793', 'RPTOR': 'ENSG00000141564',
    'CD36': 'ENSG00000135218', 'FABP4': 'ENSG00000170323', 'AKT1': 'ENSG00000142208',
    'MAPKAP1': 'ENSG00000119487', 'MLST8': 'ENSG00000148481', 'PRR5': 'ENSG00000186652',
    'PRR5L': 'ENSG00000135362', 'DEPTOR': 'ENSG00000155729',
    'SGK1': 'ENSG00000118515', 'NDRG1': 'ENSG00000104419', 'PRKCA': 'ENSG00000154229',
    'MRC1': 'ENSG00000173593', 'IL10': 'ENSG00000136634', 'CD163': 'ENSG00000177575',
    'MKI67': 'ENSG00000148773', 'PCNA': 'ENSG00000132646', 'TOP2A': 'ENSG00000131747',
    'MCM2': 'ENSG00000117748', 'CCNB1': 'ENSG00000134057', 'ADIPOQ': 'ENSG00000181092',
    'LPAR1': 'ENSG00000198160', 'TGFB1': 'ENSG00000105329',
}

U = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s?content-type=application/json'
rows = []
lines = []
def L(s=''):
    print(s); lines.append(str(s))

L('S34 · 基因 ID 权威表（Ensembl REST symbol 反查）')
L('=' * 92)
L('%-9s %-20s %-20s %-9s %-6s %-7s %s' %
  ('symbol', '权威 ENSG(chr)', '此前拟用 ENSG', '一致', 'chr', 'biotype', '此前实际测的是'))
L('-' * 92)
n_chg = 0
for s in SYMS:
    try:
        with urllib.request.urlopen(U % s, timeout=45) as r:
            o = json.load(r)
        ensg = o.get('id'); chro = str(o.get('seq_region_name')); bio = o.get('biotype')
    except Exception as e:
        L('%-9s 查询失败: %s' % (s, e)); continue
    prev = PREV.get(s)
    same = ('—' if prev is None else ('✓' if prev == ensg else '✗ 需改'))
    if prev is not None and prev != ensg:
        n_chg += 1
    # 此前拟用 ID 实际是什么（仅对不一致者查）
    actual = ''
    if prev and prev != ensg:
        try:
            with urllib.request.urlopen(
                    'https://rest.ensembl.org/lookup/id/%s?content-type=application/json' % prev,
                    timeout=45) as r:
                oo = json.load(r)
            actual = '%s (chr%s)' % (oo.get('display_name'), oo.get('seq_region_name'))
        except Exception:
            actual = '?'
    L('%-9s %-20s %-20s %-9s %-6s %-7s %s' %
      (s, '%s' % ensg, prev or '-', same, chro, (bio or '')[:7], actual))
    rows.append((s, ensg, prev or '', 'same' if (prev is None or prev == ensg) else 'changed',
                 chro, bio or '', actual))
    time.sleep(0.12)

L('')
L('★ 需改 %d 项（此前拟用 ID 实际指向别的基因）' % n_chg)

# 与 TCGA 矩阵比对
import pyreadr
_d = pyreadr.read_r(os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw", "TCGA_LPS_real_log2cpm.rds"))
_d = _d[list(_d.keys())[0]]
have = set(str(i).split('.')[0] for i in _d.index)
L('')
L('-- TCGA-SARC 矩阵命中（60660 基因）--')
miss = [r for r in rows if r[1] not in have]
L('   权威 ID 命中 %d / %d %s' % (len(rows) - len(miss), len(rows),
                                 '' if not miss else ('| ★缺: ' + ', '.join('%s(%s)' % (m[0], m[1]) for m in miss))))

with io.open(os.path.join(W, 's34_gene_ids.tsv'), 'w', encoding='utf-8-sig', newline='') as fh:
    fh.write('symbol\tensg\tsymbol_prev\tstatus\tchr\tbiotype\tprev_actual\n')
    for r in rows:
        fh.write('\t'.join(str(x) for x in r) + '\n')
io.open(os.path.join(W, 's34_gene_ids_log.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n[写出] s34_gene_ids.tsv / s34_gene_ids_log.txt')
