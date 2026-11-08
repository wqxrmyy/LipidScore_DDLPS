# -*- coding: utf-8 -*-
"""Build reference_tiers.csv for the 39 manuscript references.
Uses Crossref-verified metadata (refs_crossref.json) + P1 tier assignment.
"""
import json, os, csv, re

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
DEST = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "重构手稿", "v2", "docs")
os.makedirs(DEST, exist_ok=True)

cr = {r['ref']: r for r in json.load(open(os.path.join(OUT, 'refs_crossref.json'), encoding='utf-8'))}

# ref -> (Tier, 用途(新叙事中的角色), 证据等级说明, 处置建议)
T = {
 1:  ('Tier1','DDLPS 临床病理与预后基线','高','保留：DDLPS 定义与临床背景首选引文'),
 2:  ('Tier1','DDLPS 免疫微环境（CD4/M2 浸润与预后）','高','保留：TME 免疫表型在 DDLPS 的既有证据'),
 3:  ('Tier2','DDLPS 化疗（Cabazitaxel）II 期','中','保留但合引：治疗背景一段带过'),
 4:  ('Tier2','DDLPS Selinexor III 期','中','保留但合引：治疗背景一段带过'),
 5:  ('Tier1','DDLPS 临床病程与预后因子','高','保留：预后因子缺口论证的依据'),
 6:  ('Tier1','软组织肉瘤性别依赖预后','高','保留但改用途：仅作“性别差异已有报告”，不得再用于支撑男性特异结论（见 A1/D2）'),
 7:  ('Tier1','WD/DD 腹膜后 LPS 脂质组+转录组预后标志物','高','保留：与本课题最接近的既有脂代谢研究'),
 8:  ('Tier1','脂肪肉瘤-肌肉代谢串扰维持肿瘤生长','高','保留：脂质竞争生态位概念的直接来源'),
 9:  ('Tier2','脂代谢调控铁死亡','中','保留：机制背景，可与 20/27 合引'),
10:  ('Tier1','mTORC2–PPARγ/CD36–TAM 促 LPS 进展','高','保留：PPARG/CD36 轴唯一直接实验证据，核心引文'),
11:  ('Tier1','DDLPS 免疫图谱与新抗原','高','保留：免疫排斥表型的支撑'),
12:  ('Tier2','DDLPS 未满足临床需求（Delphi）','中','保留：引言临床动机，可与 3/4 合引'),
13:  ('Tier2','PPARγ 激动剂+曲贝替定诱导脂肪分化（黏脂肉瘤）','低','⚠️ 问题项：当前引用的是 Correction（勘误）条目，须替换为原始研究文章'),
14:  ('Tier2','年龄与亚型对可手术 LPS 结局的影响','中','保留：亚型/年龄分层背景'),
15:  ('Tier3','(原用于) 性别-肿瘤代谢差异','不可用','❌ 删除：Enzinger 1962 为组织学分类文献，被错用于支撑“sex variation in cancer metabolism”（台账 A4）；且无 DOI。替换为 38/39 或专门的性别代谢文献'),
16:  ('Tier1','脂肪肉瘤免疫微环境与临床意义','高','保留：TME 免疫框架'),
17:  ('Tier1','CD36 摄取氧化脂质致 CD8+T 细胞功能障碍','高','保留：假说 C 的机制核心引文'),
18:  ('Tier2','谷氨酰胺阻断逆转免疫逃逸','中','保留：代谢竞争旁证，可与 28 合引'),
19:  ('Tier1','PPARγ 为代谢/分化/生长的核调控因子','高','保留：PPARG 生物学基础'),
20:  ('Tier1','CD36-自噬-脂代谢互作与肿瘤进展','高','保留：CD36 综述，支撑 CD36 双向性'),
21:  ('Tier1','CD36 在肿瘤进展中的作用与靶向价值','高','保留：CD36 综述（与 20 分工：机制 vs 靶向）'),
22:  ('Tier2','软组织肉瘤细胞周期三相表型','中','保留：细胞周期激活表型的背景，可与 23 合引'),
23:  ('Tier2','CENPF 作为 LPS 预后与免疫微环境标志物','中','保留：LPS 细胞周期/免疫交叉的既有报告'),
24:  ('Tier1','表观调控代谢程序介导 LPS 细胞可塑性','高','保留：假说 B（去分化-代谢耦合）核心引文'),
25:  ('Tier2','六基因评分预测 LPS 生存与免疫浸润','中','⚠️ 谨慎：旧范式同类（六基因评分）研究。可作为“既有评分类工作”对照引用，但不得用于支撑本评分的新颖性'),
26:  ('Tier1','DDLPS 脂肪生成基因增强子表观破坏','高','保留：去分化程序的关键证据'),
27:  ('Tier2','肿瘤代谢重编程综述（糖/脂/氨基酸）','中','保留：代谢重编程总论，可压缩'),
 28:  ('Tier1','脂代谢重编程塑造肿瘤免疫微环境','高','保留：脂代谢→免疫的核心枢纽引文'),
 29:  ('Tier1','中介分析（bootstrap）方法学','高','保留：假说 C 中介分析的方法学依据，须与方法段一致'),
 30:  ('Tier1','DDLPS 免疫检查点治疗（MD Anderson 经验）','高','保留：免疫检查点临床相关性'),
 31:  ('Tier1','CDK4/6+PD-1 序贯 II 期（DDLPS）','高','保留：联合治疗依据'),
 32:  ('Tier2','肉瘤靶向治疗与检查点抑制剂','中','保留：治疗综述，可与 30/31 合引'),
 33:  ('Tier1','CDK4/6 抑制对 CDK4 扩增 LPS 的抗增殖作用','高','保留：CDK4/6 在 LPS 的机制基础（与 MDM2/CDK4 共扩增相关）'),
 34:  ('Tier2','FASN 抑制在 HCC 临床前模型','中','保留：FASN 靶向跨界证据，可压缩'),
 35:  ('Tier1','FASN 与肿瘤脂肪生成表型','高','保留：FASN/脂肪生成基础引文'),
 36:  ('Tier2','CRISPR 敲除 PPARG 对前脂肪细胞分化影响','中','保留：PPARG 功能方法学旁证'),
 37:  ('Tier1','新辅助免疫检查点阻断随机 II 期（DDLPS/UPS）','高','保留：ICI 在 DDLPS 的最高级别临床证据'),
 38:  ('Tier1','免疫应答的性别差异','高','保留但改用途：性别×免疫背景；不得再用于支撑男性特异预后结论'),
 39:  ('Tier1','雄激素受体驱动 CD8+T 干性与性别差异抗肿瘤免疫','高','保留：性别×免疫机制（替换 15 的主引文）'),
}

rows = []
for n in sorted(T):
    tier, role, ev, action = T[n]
    r = cr.get(n, {})
    rows.append({
        'Ref_No': n,
        'DOI': r.get('doi') or 'NO_DOI',
        'Crossref_status': r.get('cr_status', ''),
        'Crossref_title': (r.get('cr_title') or '')[:160],
        'Crossref_journal': r.get('cr_journal', ''),
        'Crossref_year': r.get('cr_year', ''),
        'Tier': tier,
        'Role_in_new_narrative': role,
        'Evidence_strength': ev,
        'Action': action,
    })

path = os.path.join(DEST, 'reference_tiers.csv')
with open(path, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print('wrote', path)
print('total', len(rows))
from collections import Counter
print(Counter(r['Tier'] for r in rows))
print(Counter(r['Crossref_status'] for r in rows))
