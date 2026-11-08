# 核对注释 barcodes 文件的实际列定义（grade 列异常专项）
import gzip, os, collections

G = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "gse221493")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s39_s4_barcode_cols.txt")
f = os.path.join(G, 'GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz')

L = []
with gzip.open(f, 'rt', encoding='utf-8', errors='replace') as fh:
    for i, line in enumerate(fh):
        if i >= 300000:
            break
        L.append(line.rstrip('\n').split('\t'))

ncol = len(L[0])
w = open(OUT, 'w', encoding='utf-8')
def p(s=''):
    print(s); w.write(s + '\n')

p('== 注释 barcodes 列定义核对 ==')
p('文件: %s' % os.path.basename(f))
p('行数(前 300000 取样): %d   列数: %d' % (len(L), ncol))
p('')
p('--- 前 3 行原文 ---')
for i in range(3):
    p('L%d  %s' % (i, ' | '.join(L[i])))
p('')
p('--- 逐列取值统计 ---')
for c in range(ncol):
    vals = [r[c] for r in L if len(r) > c]
    cnt = collections.Counter(vals)
    top = cnt.most_common(6)
    p('col%-3d uniq=%-7d top: %s' % (c + 1, len(cnt),
      ', '.join('%s(%d)' % (k, v) for k, v in top)))
p('')
p('--- 与 Step 3 metadata 列对照 ---')
p('Step 3 metadata 里的列名顺序: orig.ident, nCount_RNA, nFeature_RNA, patient, tissue,')
p('                            site, grade, L1, L2, L3, umap1, umap2, percent.mt, percent.hb')
w.close()
print('\nwrote ->', OUT)
