# -*- coding: utf-8 -*-
"""S18：核对审计链所用全部 Ensembl ID 与官方符号是否一致"""
import json, urllib.request, os

CAND = {
    'MDM2':   ['ENSG00000135679'],
    'FRS2':   ['ENSG00000166225'],
    'CD36':   ['ENSG00000135218'],
    'PPARG':  ['ENSG00000132170'],
    'SCD':    ['ENSG00000099194'],
    'FASN':   ['ENSG00000169710', 'ENSG00000071082'],
    'HMGCR':  ['ENSG00000113161'],
    'ACLY':   ['ENSG00000131473'],
    'ADIPOQ': ['ENSG00000181092'],
    'PLIN1':  ['ENSG00000166819'],
    'CPT1A':  ['ENSG00000110090'],
    'SREBF1': ['ENSG00000072310'],
    'CDK4':   ['ENSG00000135446'],
    'HMGA2':  ['ENSG00000149948'],
    'FABP4':  ['ENSG00000170323'],
    'LPL':    ['ENSG00000175445'],
    'CD8A':   ['ENSG00000153563', 'ENSG00000198851'],
    'PDCD1':  ['ENSG00000188389'],
    'LIPE':   ['ENSG00000079435'],
    'DDIT3':  ['ENSG00000175197'],
}

allids = sorted({i for v in CAND.values() for i in v})
body = json.dumps({'ids': allids}).encode()
req = urllib.request.Request(
    'https://rest.ensembl.org/lookup/id',
    data=body,
    headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
with urllib.request.urlopen(req, timeout=60) as r:
    res = json.loads(r.read().decode())

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s18_ensembl_ids.json")
json.dump(res, open(OUT, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

print('=== Ensembl 官方判定 ===')
for label, ids in CAND.items():
    for i in ids:
        d = res.get(i) or {}
        sym = d.get('display_name', 'NOT_FOUND')
        chr_ = d.get('seq_region_name', '?')
        desc = (d.get('description') or '')[:55]
        ok = 'OK  ' if sym == label else 'MISMATCH'
        print(f'{ok} label={label:<8} id={i} -> 官方符号 {sym:<10} chr{chr_:<3} {desc}')
print()
print('saved ->', OUT)

# 汇总不一致项
bad = [(l, i, (res.get(i) or {}).get('display_name')) for l, ids in CAND.items() for i in ids
       if (res.get(i) or {}).get('display_name') != l]
print()
print('=== 不一致清单 ===')
for l, i, s in bad:
    print(f'  {l} 被写成 {i}，实为 {s}')
