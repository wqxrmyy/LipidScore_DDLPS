# -*- coding: utf-8 -*-
import json, urllib.request, urllib.parse, os, time, re
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
UA = {'User-Agent': 'Mozilla/5.0 (geo-metadata-check; mailto:refcheck@example.org)'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def fetch(url, timeout=40):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')

out = {}
for acc in ['GSE221492', 'GSE221493']:
    uid = '200' + acc[3:]
    j = json.loads(fetch(EUT + 'esummary.fcgi?db=gds&id=%s&retmode=json' % uid))
    rec = j['result'][uid]
    keep = {k: rec.get(k) for k in ('accession', 'title', 'summary', 'taxon', 'n_samples',
                                    'gpl', 'gdstype', 'pdat', 'pubmedids', 'gse', 'bioproject',
                                    'samples', 'contact', 'suppfile')}
    out[acc] = keep
    print('===', acc, '===')
    print('title      :', rec.get('title'))
    print('taxon      :', rec.get('taxon'))
    print('n_samples  :', rec.get('n_samples'))
    print('gpl        :', rec.get('gpl'))
    print('gdstype    :', rec.get('gdstype'))
    print('pdat       :', rec.get('pdat'))
    print('pubmedids  :', rec.get('pubmedids'))
    print('bioproject :', rec.get('bioproject'))
    print('contact    :', str(rec.get('contact'))[:200])
    print('samples(3) :', str(rec.get('samples'))[:300])
    print('suppfile   :', str(rec.get('suppfile'))[:300])
    print('summary    :', rec.get('summary'))
    print()

    # sample list
    es = json.loads(fetch(EUT + 'esearch.fcgi?db=gds&term=%s[ACCN]+AND+gsm[ETYP]&retmax=200&retmode=json' % acc))
    ids = es['esearchresult']['idlist']
    print(acc, 'n GSM ids =', len(ids))
    titles = []
    for i in range(0, len(ids), 100):
        chunk = ','.join(ids[i:i + 100])
        sj = json.loads(fetch(EUT + 'esummary.fcgi?db=gds&id=%s&retmode=json' % chunk))
        for k, v in sj['result'].items():
            if k in ('uids',):
                continue
            titles.append({'gsm': v.get('accession'), 'title': v.get('title'),
                           'gpl': v.get('gpl'), 'taxon': v.get('taxon'),
                           'n': v.get('n_samples'), 'ftype': v.get('ftplink')})
        time.sleep(0.4)
    out[acc + '_samples'] = titles
    print(acc, 'n sample titles =', len(titles))
    for t in titles:
        print('   ', t['gsm'], '|', t['title'])
    print()

with open(os.path.join(OUT, 'geo_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print('saved geo_summary.json')
