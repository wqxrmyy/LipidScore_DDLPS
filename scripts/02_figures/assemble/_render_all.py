# -*- coding: utf-8 -*-
"""按最终印刷尺寸批量重渲染全部面板（v3 合规整改）。

用法：
    python _render_all.py            # 全部
    python _render_all.py fig1 figS5 # 只跑指定组
"""
import os, sys, subprocess, glob, time

ROOT = os.path.join(os.environ.get("LPS_FIG_ROOT", "."))
RSCRIPT = os.path.join(os.environ.get("R_BIN_DIR", "."), "Rscript.exe")

GROUPS = ["fig1", "fig2", "fig3", "fig4", "figS1", "figS2", "figS3", "figS4", "figS5"]

want = sys.argv[1:] or GROUPS
groups = [g for g in GROUPS if g in want]
# 按名字排序保证 fig1 在 figS5 之前（组内脚本名也排序）
skip = {"prep_data.R", "prep_gse221492.R", "themes.R", "colors.R", "assemble_lib.R"}

ok, fail = [], []
t0 = time.time()
for g in groups:
    d = os.path.join(ROOT, g)
    if not os.path.isdir(d):
        print(f"[skip] {g}: 目录不存在")
        continue
    scripts = sorted(f for f in os.listdir(d) if f.endswith(".R") and f not in skip)
    for s in scripts:
        p = os.path.join(d, s)
        r = subprocess.run([RSCRIPT, p], capture_output=True, text=True,
                           encoding="utf-8", errors="replace", cwd=d)
        tag = "ok  " if r.returncode == 0 else "FAIL"
        out = (r.stdout or "").strip().splitlines()
        msg = out[-1] if out else ""
        print(f"[{tag}] {g}/{s}  {msg}")
        if r.returncode != 0:
            fail.append((f"{g}/{s}", (r.stderr or "")[-800:]))
        else:
            ok.append(f"{g}/{s}")

print(f"\n=== {len(ok)+len(fail)} scripts, {len(fail)} failed, {time.time()-t0:.1f}s ===")
for name, err in fail:
    print(f"\n--- {name} ---\n{err}")
