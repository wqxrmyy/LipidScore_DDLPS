# -*- coding: utf-8 -*-
import os, re, json, glob
from docx import Document

MS = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "Cancer Letters投稿", "Manuscript.docx")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
os.makedirs(OUT, exist_ok=True)

d = Document(MS)
paras = [p.text.strip() for p in d.paragraphs]
with open(os.path.join(OUT, 'ms_paras.txt'), 'w', encoding='utf-8') as f:
    for i, t in enumerate(paras):
        if t:
            f.write('%d\t%s\n' % (i, t))

# locate reference section
refs = []
started = False
ref_start = None
for i, t in enumerate(paras):
    if not t:
        continue
    if re.match(r'^\s*(References|REFERENCES|参考文献)\s*$', t):
        started = True
        ref_start = i
        continue
    if started:
        if re.match(r'^\s*\d+[\.\)]\s*\S', t):
            refs.append((i, t))
print('REF SECTION header para =', ref_start)
print('N ref elements =', len(refs))
with open(os.path.join(OUT, 'ms_refs.txt'), 'w', encoding='utf-8') as f:
    for i, t in refs:
        f.write('%d\t%s\n' % (i, t))

print('N docx tables =', len(d.tables))
for i, t in enumerate(paras):
    if re.match(r'^\s*(Abstract|ABSTRACT)\s*$', t, re.I):
        print('ABSTRACT header at', i)
    if re.match(r'^\s*(Keywords|KEYWORDS)\b', t, re.I):
        print('KEYWORDS at', i, '->', t[:150])
print('TOTAL paras =', len(paras))

# section headings guess
print('--- headings-like ---')
for i, t in enumerate(paras):
    if re.match(r'^\s*\d+(\.\d+)*\s+[A-Z]', t) and len(t) < 120:
        print(i, t)
