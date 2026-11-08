# 判定注释文件 col5 的真实语义（交叉表纯度）
import gzip, os, collections

G = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "gse221493")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s39_s4_col5_ident.txt")
f = os.path.join(G, 'GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz')

rows = []
with gzip.open(f, 'rt', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        rows.append(line.rstrip('\n').split('\t'))

w = open(OUT, 'w', encoding='utf-8')
def p(s=''):
    print(s); w.write(s + '\n')

p('== col5 语义判定（n=%d 细胞）==' % len(rows))
p('')
col5 = [r[4] for r in rows]
for ci, name in [(5, 'col6/L1'), (6, 'col7/L2'), (7, 'col8/L3')]:
    lab = [r[ci] for r in rows]
    # 每个 col5 值 -> 标签分布
    m = collections.defaultdict(collections.Counter)
    for v, l in zip(col5, lab):
        m[v][l] += 1
    # 纯度：各 col5 组内主标签占比（加权）
    tot = 0; hit = 0; impure = 0; multi = 0
    for v, c in m.items():
        s = sum(c.values()); top = c.most_common(1)[0][1]
        tot += s; hit += top
        if top < s: impure += 1
        if len(c) > 1: multi += 1
    # 反向：每个标签 -> col5 值个数
    rev = collections.defaultdict(set)
    for v, l in zip(col5, lab):
        rev[l].add(v)
    nmulti = sum(1 for l, s in rev.items() if len(s) > 1)
    p('col5 x %-9s : 加权纯度 %.6f | col5 组内混入他类的组数 %d/%d | 一个标签对应多个 col5 值的标签数 %d/%d'
      % (name, hit / tot, multi, len(m), nmulti, len(rev)))

p('')
p('--- col5 全部取值（按细胞数降序）---')
cnt = collections.Counter(col5)
p('uniq = %d' % len(cnt))
for k, v in cnt.most_common():
    p('  %-10s %7d' % (k, v))

p('')
p('--- 抽样核对：col5 各值的代表性 L3 标签 ---')
lab3 = [r[7] for r in rows]
m = collections.defaultdict(collections.Counter)
for v, l in zip(col5, lab3):
    m[v][l] += 1
for k, _ in cnt.most_common(20):
    c = m[k].most_common(3)
    p('  col5=%-10s -> %s' % (k, '; '.join('%s(%d)' % (a, b) for a, b in c)))

p('')
p('--- QC 前 L1/L2/L3 类数 ---')
for ci, name in [(5, 'L1'), (6, 'L2'), (7, 'L3')]:
    p('  %s uniq = %d' % (name, len(set(r[ci] for r in rows))))

w.close()
print('\nwrote ->', OUT)
