# -*- coding: utf-8 -*-
# 串行重渲染 fig6b/6c/6d/6f（R 里 plan("multisession") 会崩 -> 一律串行）
import subprocess, os, io, sys, time
R = os.path.join(os.environ.get("R_BIN_DIR", "."), "Rscript.exe")
FIG = os.path.join(os.environ.get("LPS_FIG_ROOT", "."), "fig6")
LOG = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work", "s34_f6_render5.log")
stems = ["fig6b", "fig6c", "fig6d", "fig6f"]
out = []
for s in stems:
    out.append("===== %s =====" % s)
    t0 = time.time()
    p = subprocess.run([R, os.path.join(FIG, s + ".R")], capture_output=True)
    so = p.stdout.decode("utf-8", "replace"); se = p.stderr.decode("utf-8", "replace")
    out.append("EXIT=%d  %.1fs" % (p.returncode, time.time() - t0))
    for ln in (so + se).splitlines():
        if "Setting LC_" in ln or ln.strip().startswith(("1:", "2:", "3:", "4:")):
            continue
        out.append(ln)
txt = "\n".join(out)
io.open(LOG, "w", encoding="utf-8", newline="\n").write(txt + "\n")
print(txt)
