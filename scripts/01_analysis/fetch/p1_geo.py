# -*- coding: utf-8 -*-
import json, urllib.request, urllib.parse, os
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
UA = {'User-Agent': 'Mozilla/5.0 (research-metadata-check; mailto:refcheck@example.org)'}

def fetch(url, timeout=30):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')

res = {}
# 1) NCBI eutils GEO series summary
for acc, uid in [('GSE221492', '200221492'), ('GSE221493', '200221493')]:
    for tgt in ['self', 'gds']:
        try:
            u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
                 '?db=gds&id=%s&retmode=json' % uid)
            res['eutils_%s_%s' % (acc, tgt)] = fetch(u)[:6000]
            break
        except Exception as e:
            res['eutils_%s_%s' % (acc, tgt)] = 'ERR %s: %s' % (type(e).__name__, e)

# 2) GEO sample list via eutils esearch on gsm
for acc in ['GSE221492', 'GSE221493']:
    try:
        u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
             '?db=gds&term=%s[ACCN]+AND+gse[ETYP]&retmode=json' % acc)
        res['esearch_%s' % acc] = fetch(u)[:2000]
    except Exception as e:
        res['esearch_%s' % acc] = 'ERR %s: %s' % (type(e).__name__, e)

# 3) EBI BioStudies
for acc in ['GSE221492', 'GSE221493']:
    try:
        u = 'https://www.ebi.ac.uk/biostudies/api/v1/search?query=%s&pageSize=5' % acc
        res['biostudies_%s' % acc] = fetch(u)[:4000]
    except Exception as e:
        res['biostudies_%s' % acc] = 'ERR %s: %s' % (type(e).__name__, e)

# 4) OmicsDI
for acc in ['GSE221492', 'GSE221493']:
    try:
        u = ('https://www.omicsdi.org/ws/dataset/search?query=%s&size=5' % acc)
        res['omicsdi_%s' % acc] = fetch(u)[:4000]
    except Exception as e:
        res['omicsdi_%s' % acc] = 'ERR %s: %s' % (type(e).__name__, e)

with open(os.path.join(OUT, 'geo_raw.json'), 'w', encoding='utf-8') as f:
    json.dump(res, f, ensure_ascii=False, indent=1)
for k, v in res.items():
    print('=====', k, '=====')
    print(str(v)[:900])
    print()
