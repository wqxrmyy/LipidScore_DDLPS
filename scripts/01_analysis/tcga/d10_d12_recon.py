# -*- coding: utf-8 -*-
"""D10-D12 勘察：_work 清单 + GEO 元数据 + TCGA 临床全列 + 评分列 + CN 簇分布"""
import os, csv, json, collections

W = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
RAW = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "data", "raw")

print('=== _work files ===')
for f in sorted(os.listdir(W)):
    p = os.path.join(W, f)
    if os.path.isfile(p):
        print(f'  {os.path.getsize(p):>10} {f}')

print('\n=== json inventories ===')
for f in sorted(os.listdir(W)):
    if f.endswith('.json'):
        try:
            j = json.load(open(os.path.join(W, f), encoding='utf-8'))
            print(f'  {f}: keys={list(j.keys())[:12]}' if isinstance(j, dict) else f'  {f}: list len={len(j)}')
        except Exception as e:
            print(f'  {f}: ERR {e}')

print('\n=== TCGA_Score.csv header/rows ===')
sc = list(csv.reader(open(os.path.join(RAW, 'TCGA_Score.csv'), encoding='utf-8-sig')))
print('nrows', len(sc), 'HDR', sc[0], 'row1', sc[1] if len(sc) > 1 else None)

print('\n=== TCGA_clinical_real.csv ALL columns ===')
cl = list(csv.reader(open(os.path.join(RAW, 'TCGA_clinical_real.csv'), encoding='utf-8-sig')))
print('nrows', len(cl))
for i, c in enumerate(cl[0]):
    print(f'  [{i}] {c}')
