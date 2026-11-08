# -*- coding: utf-8 -*-
import json, os
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
j = json.load(open(os.path.join(OUT, 'd15_score_revision.json'), encoding='utf-8'))
d11 = json.load(open(os.path.join(OUT, 'd11_direction_check.json'), encoding='utf-8'))
L = []
L.append('===== D15 重算基线核对 =====')
L.append('recompute_check: ' + json.dumps(j['tcga_recompute_check'], ensure_ascii=False))
L.append(f"tcga_n={j['tcga_n']}  OS n={j.get('tcga_OS_n')} events={j.get('tcga_OS_events')}")
L.append('')
L.append('--- orig6 (D15 由 log2cpm 重算) ---')
L.append(json.dumps(j['TCGA']['orig6'], ensure_ascii=False, indent=1))
L.append('')
L.append('--- d11 既有值（用 TCGA_Score.csv） ---')
L.append('D_robustness_extra: ' + json.dumps(d11['D_robustness_extra'], ensure_ascii=False, indent=1))
L.append(f"OS n={d11.get('OS_n')} events={d11.get('OS_events')}")
L.append('Cox_score_perSD: ' + json.dumps(d11.get('Cox_score_perSD'), ensure_ascii=False))
L.append('Cox_adj: ' + json.dumps(d11.get('Cox_score_adj_grade_mitotic'), ensure_ascii=False))
L.append('')
L.append('--- 全变体 TCGA 关键指标 ---')
hdr = f"{'variant':16} {'T1 r':>7} {'T2 p.purity':>12} {'T3 full':>8} {'grade rho':>10} {'mitot rho':>10} {'HR/SD':>7} {'HRadj':>7} {'C':>7}"
L.append(hdr)
for v, d in j['TCGA'].items():
    t1 = d['T1_r_CNcluster']['r'] if d.get('T1_r_CNcluster') else None
    t2 = d['T2_partial_CN_purity']['r'] if d.get('T2_partial_CN_purity') else None
    t3 = d['T3_partial_CN_purity_grade_mitotic']['r'] if d.get('T3_partial_CN_purity_grade_mitotic') else None
    g4 = d['T4_spearman_grade']['rho'] if d.get('T4_spearman_grade') else None
    m5 = d['T5_spearman_mitotic']['rho'] if d.get('T5_spearman_mitotic') else None
    hr = d.get('T6_Cox_unadj', {}).get('HR')
    ha = d.get('T7_Cox_adj', {}).get('HR')
    c = d.get('T8_harrell_C')
    L.append(f'{v:16} {t1!s:>7} {t2!s:>12} {t3!s:>8} {g4!s:>10} {m5!s:>10} {hr!s:>7} {ha!s:>7} {c!s:>7}')
L.append('')
L.append('--- 全变体 GSE221492 关键指标 ---')
L.append(f"{'variant':16} {'AUC DDvWD':>10} {'p':>9} {'paired d':>9} {'neg/8':>6} {'W p':>8} {'MDM2 rho':>9} {'CD36 rho':>9} {'CD8A rho':>9}")
for v, d in j['GSE221492'].items():
    if 'error' in d:
        L.append(f'{v:16} ERROR {d["error"]}'); continue
    a = d.get('G1_AUC_DD_vs_WD') or {}
    p2 = d.get('G2_paired') or {}
    L.append(f"{v:16} {a.get('auc')!s:>10} {a.get('p')!s:>9} {p2.get('mean_delta')!s:>9} "
             f"{str(p2.get('n_neg'))+'/'+str(p2.get('n_pairs')):>6} {p2.get('wilcoxon_p')!s:>8} "
             f"{(d.get('G3_spearman_MDM2') or {}).get('rho')!s:>9} {(d.get('G4_spearman_CD36') or {}).get('rho')!s:>9} "
             f"{(d.get('G5_spearman_CD8A') or {}).get('rho')!s:>9}")
L.append('')
L.append('--- TCGA 单基因诊断（各分量自身的性质） ---')
L.append(json.dumps(j['tcga_gene_diagnostics'], ensure_ascii=False, indent=1))
L.append('')
L.append('--- 裁决表 ---')
L.append(json.dumps(j['verdict'], ensure_ascii=False, indent=1))
open(os.path.join(OUT, 'd15_report.txt'), 'w', encoding='utf-8').write('\n'.join(L))
print('written', len('\n'.join(L)))
