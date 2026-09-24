#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_hygiene.py — assert that this repository tracks nothing it should not.

Why this exists
---------------
The repository is released as **code + governance locks only**. Three classes of
artefact must never be committed, because committing them is irreversible once
pushed (GitHub keeps the object forever, and public corpora harvest it):

  1. Manuscript text  (.docx/.doc/.html/.htm/.tex)  -> plagiarism-corpus risk
  2. Figures          (.pdf/.tif/.png/.jpg/.eps)    -> figure-library risk
  3. Internal working snapshots (_work/, reproducibility/logs/, _backup_*/)
     -> contain internal process artefacts

`.gitignore` alone is a *convention*: `git add -f`, a renamed pattern, or a
sandbox that regenerates `.gitignore` can silently defeat it. This script checks
the **actual index** (`git ls-files`), not the ignore rules, so it catches what
`.gitignore` misses.

Usage
-----
  python scripts/00_repo_hygiene/check_hygiene.py            # check HEAD/index
  python scripts/00_repo_hygiene/check_hygiene.py --root DIR # explicit repo

Exit code 0 = clean, 1 = violations found. Wire it into a pre-push hook or CI.
"""

import argparse
import os
import re
import subprocess
import sys

# --- forbidden tracked paths -------------------------------------------------
# Match against the POSIX-style path as reported by `git ls-files`.
FORBIDDEN_EXT = {
    '.docx', '.doc', '.odt', '.rtf', '.tex', '.html', '.htm',
    '.pdf', '.tif', '.tiff', '.png', '.jpg', '.jpeg', '.eps',
}
FORBIDDEN_DIR = (
    'manuscript/',
    'reproducibility/logs/',
    '_work/',
    '_internal/',
    'data/',
)
FORBIDDEN_ANYWHERE = (
    r'(^|/)\.Rproj\.user/',
    r'(^|/)\.RData$',
    r'(^|/)\.Rhistory$',
    r'(^|/)\.Ruserdata$',
    r'(^|/)_backup_[^/]*/',
)

SIZE_LIMIT_MB = 5.0


def git(root, *args):
    p = subprocess.run(['git'] + list(args), cwd=root,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError('git %s failed: %s' % (' '.join(args),
                                                  p.stderr.decode('utf-8', 'replace')))
    return p.stdout.decode('utf-8', 'replace')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default=None,
                    help='repository root (default: inferred from this file)')
    ap.add_argument('--size-limit-mb', type=float, default=SIZE_LIMIT_MB)
    a = ap.parse_args()

    if a.root:
        root = os.path.abspath(a.root)
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        root = os.path.abspath(os.path.join(here, '..', '..'))

    print('=' * 78)
    print('仓库卫生检查 / repository hygiene check')
    print('  root: %s' % root)
    print('=' * 78)

    tracked = [ln for ln in git(root, 'ls-files').splitlines() if ln.strip()]
    print('\n被跟踪文件数 / tracked files: %d' % len(tracked))

    problems = []

    # 1) forbidden extensions
    bad_ext = [f for f in tracked
               if os.path.splitext(f.lower())[1] in FORBIDDEN_EXT]
    for f in bad_ext:
        problems.append(('禁用扩展名 / forbidden extension', f))

    # 2) forbidden directories (as a path prefix anywhere in the tree)
    bad_dir = []
    for f in tracked:
        posix = f.replace('\\', '/')
        for d in FORBIDDEN_DIR:
            if posix.startswith(d) or ('/' + d) in ('/' + posix):
                bad_dir.append((f, d))
                break
    for f, d in bad_dir:
        problems.append(('禁用目录 / forbidden directory (%s)' % d.rstrip('/'), f))

    # 3) forbidden patterns anywhere
    for f in tracked:
        posix = f.replace('\\', '/')
        for pat in FORBIDDEN_ANYWHERE:
            if re.search(pat, posix):
                problems.append(('禁用模式 / forbidden pattern', f))
                break

    # 4) oversized tracked files
    big = []
    for f in tracked:
        p = os.path.join(root, f.replace('/', os.sep))
        if os.path.isfile(p):
            mb = os.path.getsize(p) / 1024.0 / 1024.0
            if mb > a.size_limit_mb:
                big.append((f, mb))
    for f, mb in big:
        problems.append(('超出大小上限 / oversized (%.1f MB)' % mb, f))

    # --- report --------------------------------------------------------------
    print('\n' + '-' * 78)
    if problems:
        print('❌ 发现 %d 项违规 / %d violation(s):' % (len(problems), len(problems)))
        seen = set()
        for why, f in problems:
            key = (why, f)
            if key in seen:
                continue
            seen.add(key)
            print('   [%s]  %s' % (why, f))
    else:
        print('✅ 未被跟踪任何稿件正文 / 图片 / 内部快照 / 第三方数据')
        print('✅ no manuscript text, figures, internal snapshots or data tracked')

    # --- positive evidence: what IS tracked ---------------------------------
    print('\n' + '-' * 78)
    print('已跟踪内容概览 / what IS tracked:')
    tops = {}
    for f in tracked:
        top = f.split('/')[0] if '/' in f else '(root)'
        tops[top] = tops.get(top, 0) + 1
    for k in sorted(tops):
        print('   %-28s %4d files' % (k, tops[k]))
    total_mb = sum(os.path.getsize(os.path.join(root, f.replace('/', os.sep)))
                   for f in tracked
                   if os.path.isfile(os.path.join(root, f.replace('/', os.sep)))) / 1024.0 / 1024.0
    print('   总计 / total: %d files, %.2f MB' % (len(tracked), total_mb))

    # --- commit count --------------------------------------------------------
    try:
        n = git(root, 'rev-list', '--all', '--count').strip()
        print('   提交数 / commits: %s' % n)
    except RuntimeError:
        pass

    print('=' * 78)
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
