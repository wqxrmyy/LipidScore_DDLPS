# -*- coding: utf-8 -*-
import json, urllib.request, urllib.parse, os, gzip, io, re
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
UA = {'User-Agent': 'Mozilla/5.0 (geo-metadata-check; mailto:refcheck@example.org)'}

def fetch(url, timeout=60, binary=False):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        b = r.read()
    return b if binary else b.decode('utf-8', 'replace')

# 1) list FTP series dir for GSE221492 / GSE221493
for acc, pre in [('GSE221492', 'GSE221'), ('GSE221493', 'GSE221')]:
    try:
        u = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/' % (pre, acc)
        html = fetch(u)
        files = re.findall(r'href="([^"]+)"', html)
        print('=== FTP', acc, 'matrix dir ===')
        for f in files:
            print('   ', f)
        # try to read series matrix header
        for f in files:
            if f.endswith('_series_matrix.txt.gz'):
                b = fetch(u + f, binary=True)
                txt = gzip.decompress(b).decode('utf-8', 'replace')
                lines = txt.splitlines()
                print('--- header of', f, '(%d lines) ---' % len(lines))
                for ln in lines[:80]:
                    if ln.startswith('!Sample_') or ln.startswith('!Series_'):
                        print(ln[:400])
                break
    except Exception as e:
        print('ERR', acc, type(e).__name__, e)
    print()

# 2) Crossref for PMID-associated DOI: search by title
def cr_query(q):
    u = 'https://api.crossref.org/works?rows=2&query.bibliographic=' + urllib.parse.quote(q) + '&mailto=refcheck@example.org'
    req = urllib.request.Request(u, headers={'User-Agent': 'RefCheck/1.0 (mailto:refcheck@example.org)'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode('utf-8'))

try:
    j = cr_query('Cellular origin and clonal evolution of human dedifferentiated liposarcoma')
    for it in j['message']['items']:
        print('CR:', it.get('DOI'), '|', (it.get('container-title') or [''])[0], '|',
              (it.get('title') or [''])[0][:90], '|', it.get('issued', {}).get('date-parts'))
except Exception as e:
    print('ERR crossref', e)
