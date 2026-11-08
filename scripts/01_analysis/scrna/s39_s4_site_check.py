# site（解剖部位）在配对患者内部是否为常量 + Adipo diff 分布预备
import gzip, os, collections

G = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "gse221493")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s39_s4_site_check.txt")
f = os.path.join(G, 'GSE221493_sc_LPS_atlas_normalized_counts_barcodes.tsv.gz')

PAT, TIS, SITE, L2 = 1, 2, 3, 6
mn = collections.defaultdict(collections.Counter)   # (pat, tis) -> site
ad = collections.Counter()                          # (pat, tis) -> Adipo diff n
tot = collections.Counter()
with gzip.open(f, 'rt', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        r = line.rstrip('\n').split('\t')
        k = (r[PAT], r[TIS])
        mn[k][r[SITE]] += 1
        tot[k] += 1
        if r[L2] == 'Adipo diff':
            ad[k] += 1

w = open(OUT, 'w', encoding='utf-8')
def p(s=''):
    print(s); w.write(s + '\n')

PAIRS = ['P1', 'P4', 'P7', 'P10', 'P14', 'P15', 'P16']
p('== site 在配对患者内部的一致性 ==')
p('%-6s %-12s %-18s %8s | %-12s %-18s %8s | site一致?' % (
    'pat', 'WD tissue', 'WD site', 'nWD', 'DD tissue', 'DD site', 'nDD'))
ok = 0
for pa in PAIRS:
    wd = mn.get((pa, 'DDLPS-WD'), collections.Counter())
    dd = mn.get((pa, 'DDLPS-DD'), collections.Counter())
    ws = ','.join(wd); ds = ','.join(dd)
    same = (set(wd) == set(dd))
    ok += same
    p('%-6s %-12s %-18s %8d | %-12s %-18s %8d | %s' % (
        pa, 'DDLPS-WD', ws, tot.get((pa, 'DDLPS-WD'), 0),
        'DDLPS-DD', ds, tot.get((pa, 'DDLPS-DD'), 0), 'YES' if same else 'NO'))
p('')
p('配对内 site 一致的患者数: %d / %d' % (ok, len(PAIRS)))
p('')
p('== 各样本 site 一览 ==')
for (pa, ti), c in sorted(mn.items()):
    p('  %-5s %-12s %-18s %7d' % (pa, ti, ','.join(c), tot[(pa, ti)]))
p('')
p('== Adipo diff 在 7 个配对患者中的分布（col7/L2 == "Adipo diff"）==')
p('%-6s %8s %8s %8s %8s' % ('pat', 'nWD', 'nDD', 'adiWD', 'adiDD'))
for pa in PAIRS:
    p('%-6s %8d %8d %8d %8d' % (
        pa, tot.get((pa, 'DDLPS-WD'), 0), tot.get((pa, 'DDLPS-DD'), 0),
        ad.get((pa, 'DDLPS-WD'), 0), ad.get((pa, 'DDLPS-DD'), 0)))
w.close()
print('\nwrote ->', OUT)
