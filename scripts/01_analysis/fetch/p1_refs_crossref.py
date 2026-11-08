# -*- coding: utf-8 -*-
import os, re, json, time, urllib.request, urllib.parse

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
MS_REFS = os.path.join(OUT, 'ms_refs.txt')

# 1. re-parse references from manuscript paragraphs
from docx import Document
MS = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "Cancer Letters投稿", "Manuscript.docx")
d = Document(MS)
paras = [p.text.strip() for p in d.paragraphs]

refs = {}
started = False
for i, t in enumerate(paras):
    if not t:
        continue
    if re.match(r'^\s*References\s*$', t):
        started = True
        continue
    if started:
        m = re.match(r'^\[(\d+)\]\s*(.+)$', t.replace('\t', ' '))
        if m:
            refs[int(m.group(1))] = m.group(2).strip()
        if re.match(r'^\s*(Figure Legends|Tables?|Supplementary)\b', t):
            break
print('N refs =', len(refs), 'keys', sorted(refs)[:5], '...', sorted(refs)[-3:])

def get_doi(s):
    m = re.search(r'(10\.\d{4,9}/[^\s]+)', s)
    if not m:
        return None
    doi = m.group(1).rstrip('.').rstrip(')')
    return doi

rows = []
for n in sorted(refs):
    s = refs[n]
    doi = get_doi(s)
    # journal & year heuristic
    jm = re.search(r'\.\s*([A-Z][A-Za-z\(\)\.\s&\-]{2,60}?)\s+(\d{1,4})\s*\((\d{4})\)', s)
    year = None
    ym = re.findall(r'\((19|20)\d{2}\)', s)
    if ym:
        year = re.findall(r'\(((?:19|20)\d{2})\)', s)[-1]
    rows.append({'ref': n, 'doi': doi, 'year': year, 'raw': s})

with open(MS_REFS, 'w', encoding='utf-8') as f:
    for r in rows:
        f.write('[%d]\tDOI=%s\tYEAR=%s\n' % (r['ref'], r['doi'], r['year']))
print('refs without DOI:', [r['ref'] for r in rows if not r['doi']])
print('refs with DOI:', sum(1 for r in rows if r['doi']))

# 2. Crossref verification
MAILTO = 'refcheck@example.org'
results = []
for r in rows:
    rec = dict(r)
    if not r['doi']:
        rec.update(cr_status='NO_DOI_IN_REF', cr_title='', cr_journal='', cr_year='')
        results.append(rec)
        continue
    url = 'https://api.crossref.org/works/' + urllib.parse.quote(r['doi']) + '?mailto=' + MAILTO
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'RefCheck/1.0 (mailto:%s)' % MAILTO})
        with urllib.request.urlopen(req, timeout=25) as resp:
            j = json.loads(resp.read().decode('utf-8'))
        msg = j.get('message', {})
        title = (msg.get('title') or [''])[0]
        cont = (msg.get('container-title') or [''])
        journal = cont[0] if cont else ''
        yr = ''
        for k in ('published-print', 'published-online', 'issued', 'created'):
            if msg.get(k, {}).get('date-parts'):
                yr = str(msg[k]['date-parts'][0][0])
                break
        rec.update(cr_status='FOUND', cr_title=title, cr_journal=journal, cr_year=yr,
                   cr_type=msg.get('type', ''), cr_doi=msg.get('DOI', ''))
    except urllib.error.HTTPError as e:
        rec.update(cr_status='HTTP_%s' % e.code, cr_title='', cr_journal='', cr_year='')
    except Exception as e:
        rec.update(cr_status='ERR_%s' % type(e).__name__, cr_title='', cr_journal='', cr_year='')
    results.append(rec)
    time.sleep(0.4)

with open(os.path.join(OUT, 'refs_crossref.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=1)

print('\n=== CROSSREF SUMMARY ===')
ok = sum(1 for r in results if r['cr_status'] == 'FOUND')
print('FOUND %d / %d' % (ok, len(results)))
for r in results:
    print('[%s] %-12s %s | %s | %s' % (r['ref'], r['cr_status'], (r['cr_year'] or '?'), (r['cr_journal'] or '-')[:38], (r['cr_title'] or '')[:70]))
