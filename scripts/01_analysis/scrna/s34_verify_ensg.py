# -*- coding: utf-8 -*-
"""
S34 · Ensembl REST 核验 Fig6 所需全部 ENSG ID（沿用项目 s18_ensembl_ids.json 的同一方法）

铁律：不凭印象硬编码基因 ID（本项目曾因 FASN 探针错配整条撤回 S14）。
做法：POST https://rest.ensembl.org/lookup/id -> 校验 display_name == 期望 symbol。
输出：_work/s34_ensg_verify.json / .txt（供 s34_nocontra.py 读取，避免二次硬编码）
"""
import json, os, urllib.request

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
EXPECT = {
    # 已锁（D18 六基因 + 既有）
    'SREBF1': 'ENSG00000072310', 'SCD': 'ENSG00000099194', 'HMGCR': 'ENSG00000113161',
    'ACLY': 'ENSG00000131473', 'PPARG': 'ENSG00000132170', 'FASN': 'ENSG00000169710',
    # 轴基因（既有）
    'RICTOR': 'ENSG00000133030', 'MTOR': 'ENSG00000198793', 'RPTOR': 'ENSG00000141564',
    'CD36': 'ENSG00000135218', 'FABP4': 'ENSG00000170323', 'AKT1': 'ENSG00000142208',
    # ★ 本轮新增（须核验）
    'MAPKAP1': 'ENSG00000119487', 'MLST8': 'ENSG00000148481', 'PRR5': 'ENSG00000186652',
    'PRR5L': 'ENSG00000135362', 'DEPTOR': 'ENSG00000155729',
    'SGK1': 'ENSG00000118515', 'NDRG1': 'ENSG00000104419', 'PRKCA': 'ENSG00000154229',
    # 对照 / 分层用
    'MRC1': 'ENSG00000173593', 'IL10': 'ENSG00000136634', 'CD163': 'ENSG00000177575',
    'MKI67': 'ENSG00000148773', 'PCNA': 'ENSG00000132646', 'TOP2A': 'ENSG00000131747',
    'MCM2': 'ENSG00000117748', 'CCNB1': 'ENSG00000134057', 'ADIPOQ': 'ENSG00000181092',
    'LPAR1': 'ENSG00000198160', 'TGFB1': 'ENSG00000105329',
}

ids = sorted(set(EXPECT.values()))
url = 'https://rest.ensembl.org/lookup/id'
body = json.dumps({'ids': ids}).encode()
req = urllib.request.Request(url, data=body,
                            headers={'Content-Type': 'application/json',
                                     'Accept': 'application/json',
                                     'User-Agent': 'wb-fig6-idcheck/1.0'})
out = {}
try:
    with urllib.request.urlopen(req, timeout=60) as r:
        out = json.load(r)
    err = None
except Exception as e:
    err = '%s: %s' % (type(e).__name__, e)

lines = []
def L(s=''):
    print(s); lines.append(str(s))

L('S34 · Ensembl REST ID 核验')
L('=' * 78)
L('端点: %s | 请求 %d 个 ID' % (url, len(ids)))
if err:
    L('★ 请求失败: %s' % err)
else:
    L('返回 %d 条' % len(out))
L('')
L('%-10s %-20s %-16s %-10s %s' % ('symbol', '本脚本拟用 ENSG', 'Ensembl display', '一致', 'biotype'))
L('-' * 84)
ok = bad = miss = 0
res = {}
for sym, e in sorted(EXPECT.items()):
    rec = out.get(e) if isinstance(out, dict) else None
    if not rec:
        L('%-10s %-20s %-16s %-10s %s' % (sym, e, '(未返回)', '★缺', '-'))
        miss += 1
        res[sym] = {'ensg': e, 'display': None, 'agree': False}
        continue
    disp = rec.get('display_name')
    bio = rec.get('biotype')
    ag = (disp == sym)
    L('%-10s %-20s %-16s %-10s %s' % (sym, e, disp or '-', '✓' if ag else '✗ 冲突', bio or '-'))
    if ag:
        ok += 1
    else:
        bad += 1
    res[sym] = {'ensg': e, 'display': disp, 'biotype': bio, 'agree': bool(ag), 'raw': rec}

L('')
L('一致 %d | 冲突 %d | 未返回 %d' % (ok, bad, miss))
L('')
if bad:
    L('★ 冲突项（须改口径锁后再用）:')
    for s, v in res.items():
        if not v['agree'] and v['display']:
            L('   %s: 拟用 %s 实为 %s' % (s, v['ensg'], v['display']))
if miss and not err:
    L('★ 未返回项: %s（Ensembl 对已退休 ID 会返回空）' %
      ', '.join(s for s, v in res.items() if v['display'] is None))

json.dump({'endpoint': url, 'error': err, 'result': res},
          open(os.path.join(W, 's34_ensg_verify.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
open(os.path.join(W, 's34_ensg_verify.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n[写出] s34_ensg_verify.json / .txt')
