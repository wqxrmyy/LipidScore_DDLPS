# 分析计划（预注册，执行前锁定）

## 主要终点
1. 脂质代谢评分在 WD vs DD 组织中的差异（GSE221492 bulk）
2. 评分在单细胞亚群中的分布（GSE221493 scRNA）
3. PPARG2 靶基因程序与评分的解耦联

## 次要终点
1. 细胞组成差异是否驱动 bulk 评分差异
2. TGF-β 通路在 DD 微环境中的活性
3. 轨迹分析中 PPARG2 靶基因的变化方向

## 统计方法
- bulk: paired Wilcoxon（WD vs DD 配对）
- scRNA: AUCell 评分 + 混合效应模型
- 多重比较: BH 校正

## 停止规则
- 如果 Harmony 校正后 WD 和 DD 样本仍完全分离 → 不强行合并，改为分层分析
- 如果某基因在 scRNA 中 dropout > 80% → 触发限制记录，不静默移除

## 版本
- version: v1
- locked_at: "YYYY-MM-DD"
