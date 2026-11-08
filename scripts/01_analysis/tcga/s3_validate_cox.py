# -*- coding: utf-8 -*-
"""校验自实现 Cox（Breslow+Newton）与 lifelines 的一致性。"""
import sys, os
sys.path.insert(0, os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work"))
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from s3_mediation import cox_fit, load_data

df = load_data()
t = df["time"].values.astype(float)
d = df["event"].values.astype(float)
Xz = ((df["X"] - df["X"].mean()) / df["X"].std(ddof=1)).values
PP = ((df["PPARG"] - df["PPARG"].mean()) / df["PPARG"].std(ddof=1)).values

configs = {
    "单变量 X": Xz[:, None],
    "单变量 PPARG": PP[:, None],
    "X + PPARG": np.column_stack([Xz, PP]),
}

print("%-14s %-28s %-28s" % ("模型", "自实现 beta (se)", "lifelines beta (se)"))
for name, X in configs.items():
    b, se, ll, ok = cox_fit(X, t, d)
    lf = pd.DataFrame(X, columns=["v%d" % i for i in range(X.shape[1])])
    lf["T"], lf["E"] = t, d
    cph = CoxPHFitter().fit(lf, "T", "E")
    print("-" * 74)
    for j, col in enumerate(["v%d" % i for i in range(X.shape[1])]):
        mine = "%+.6f (%.6f)" % (b[j], se[j])
        theirs = "%+.6f (%.6f)" % (cph.params_[col], cph.standard_errors_[col])
        flag = "OK" if abs(b[j] - cph.params_[col]) < 1e-4 else "**DIFF**"
        print("%-14s %-28s %-28s %s" % (name if j == 0 else "", mine, theirs, flag))
    print("%-14s loglik mine=%.6f  lifelines=%.6f" % ("", ll, cph.log_likelihood_))
