# -*- coding: utf-8 -*-
"""
S3 · 规格稳健性网格
--------------------------------------------------------------------------
手稿未提供中介分析代码，故其"未检出中介"的结论无法复现。
本脚本在同一套原始数据上枚举合理规格，检验结论对规格的敏感性——
这是判断阴性断言是否可用的必要条件。

枚举维度：
  X 形式     : 连续(z)  /  中位数二分 High-Low
  M 标准化   : 是 / 否
  协变量     : 无 / 年龄 / 年龄+性别
  间接效应   : 系数乘积 a*b  /  系数差 c-c'
  CI 方法    : 百分位法 / BCa
"""
import sys
import os
import json
import itertools
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3_mediation import load_data, ols_fit, cox_fit, pct_ci

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
SEED = 20260917
NBOOT = 1000


def zs(v):
    v = np.asarray(v, float)
    s = v.std(ddof=1)
    return (v - v.mean()) / s if s > 0 else v - v.mean()


def fit_mediation(X, M, t, d, covs=None, pt="prod"):
    """返回 (a, b, c, cprime, nie)。X/M 均已按规格预处理。"""
    C = np.zeros((len(X), 0)) if covs is None else np.asarray(covs, float)
    if C.ndim == 1:
        C = C[:, None]
    # a 路径
    ab, _ = ols_fit(M, np.column_stack([X, C]))
    a = ab[1]
    # b 路径（X + M + 协变量）
    bb, _, _, _ = cox_fit(np.column_stack([X, M, C]), t, d)
    b, cprime = bb[1], bb[0]
    # 总效应（X + 协变量）
    cb, _, _, _ = cox_fit(np.column_stack([X, C]), t, d)
    c = cb[0]
    nie = a * b if pt == "prod" else (c - cprime)
    return a, b, c, cprime, nie


def bca_ci(theta_hat, boot_vals, jack_vals, alpha=0.05):
    """BCa 区间。"""
    b = np.asarray(boot_vals, float)
    b = b[np.isfinite(b)]
    j = np.asarray(jack_vals, float)
    j = j[np.isfinite(j)]
    if b.size < 20 or j.size < 3:
        return pct_ci(b, 100 * alpha / 2, 100 * (1 - alpha / 2))
    from scipy import stats as st
    prop = np.mean(b < theta_hat)
    prop = min(max(prop, 1.0 / (b.size + 1)), 1 - 1.0 / (b.size + 1))
    z0 = st.norm.ppf(prop)
    jm = j.mean()
    num = np.sum((jm - j) ** 3)
    den = 6.0 * (np.sum((jm - j) ** 2)) ** 1.5
    a_hat = num / den if den > 0 else 0.0
    zl, zu = st.norm.ppf(alpha / 2), st.norm.ppf(1 - alpha / 2)

    def adj(z):
        return st.norm.cdf(z0 + (z0 + z) / (1 - a_hat * (z0 + z)))
    pl, pu = adj(zl), adj(zu)
    return float(np.percentile(b, 100 * pl)), float(np.percentile(b, 100 * pu))


def main():
    df = load_data()
    t = df["time"].values.astype(float)
    d = df["event"].values.astype(float)
    n = len(df)
    # load_data 未取 age/sex，这里直接从临床表补
    clin = pd.read_csv(os.path.join(
        os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1"), "data", "results",
        "03_clinical_outcomes.csv"), index_col=0)
    clin = clin.loc[df.index]
    age = pd.to_numeric(clin["Age"], errors="coerce").values
    gender = (clin["Gender"].astype(str).str.lower() == "male").astype(float).values
    rng = np.random.default_rng(SEED)

    rows = []
    specs = list(itertools.product(
        ["continuous", "median_split"],   # X 形式
        ["M_std", "M_raw"],               # M 标准化
        ["none", "age", "age+sex"],       # 协变量
    ))
    for mediator in ("PPARG", "CD36"):
        Mraw = df[mediator].values
        print("=" * 78)
        print("介质 = %s" % mediator)
        print("%-14s %-8s %-10s %-10s %-22s %-7s %s"
              % ("X形式", "M", "协变量", "效应定义", "NIE 95%CI", "PM", "结论"))
        print("-" * 78)
        for xmode, mmode, covmode in specs:
            Xb = zs(df["X"].values) if xmode == "continuous" else \
                (df["X"].values > np.median(df["X"].values)).astype(float)
            M = zs(Mraw) if mmode == "M_std" else Mraw - Mraw.mean()
            covs = {"none": None,
                    "age": age - np.nanmean(age),
                    "age+sex": np.column_stack([age - np.nanmean(age), gender])}[covmode]
            # 点估计 + jackknife + bootstrap
            a, b, c, cp, nie = fit_mediation(Xb, M, t, d, covs)
            boot = []
            for _ in range(NBOOT):
                idx = rng.integers(0, n, n)
                try:
                    _, _, _, _, v = fit_mediation(Xb[idx], M[idx], t[idx], d[idx],
                                                  None if covs is None else covs[idx])
                    if np.isfinite(v):
                        boot.append(v)
                except Exception:
                    continue
            jack = []
            for i in range(n):
                keep = np.ones(n, bool)
                keep[i] = False
                try:
                    _, _, _, _, v = fit_mediation(Xb[keep], M[keep], t[keep], d[keep],
                                                  None if covs is None else covs[keep])
                    jack.append(v)
                except Exception:
                    continue
            lo_p, hi_p = pct_ci(np.asarray(boot))
            lo_b, hi_b = bca_ci(nie, boot, jack)
            sig = (lo_b > 0 or hi_b < 0)
            pm = nie / c if c != 0 else np.nan
            print("%-14s %-8s %-10s %-10s [%+.3f, %+.3f]  %+7.0f%%  %s"
                  % (xmode, mmode, covmode, "a*b",
                     lo_p, hi_p, 100 * pm,
                     "显著" if sig else "不显著（含0）"))
            rows.append({"mediator": mediator, "x_mode": xmode, "m_mode": mmode,
                         "covariates": covmode, "effect_def": "a*b",
                         "a": a, "b": b, "c": c, "cprime": cp, "nie": nie, "pm": pm,
                         "ci_pct_lo": lo_p, "ci_pct_hi": hi_p,
                         "ci_bca_lo": lo_b, "ci_bca_hi": hi_b,
                         "significant_bca": bool(sig),
                         "significant_pct": bool(lo_p > 0 or hi_p < 0)})

        # 效应定义：系数差
        Xb = zs(df["X"].values)
        M = zs(Mraw)
        a, b, c, cp, nie = fit_mediation(Xb, M, t, d, None, pt="diff")
        boot = []
        for _ in range(NBOOT):
            idx = rng.integers(0, n, n)
            try:
                _, _, _, _, v = fit_mediation(Xb[idx], M[idx], t[idx], d[idx], None,
                                              pt="diff")
                if np.isfinite(v):
                    boot.append(v)
            except Exception:
                continue
        lo_p, hi_p = pct_ci(np.asarray(boot))
        pm = nie / c if c != 0 else np.nan
        print("%-14s %-8s %-10s %-10s [%+.3f, %+.3f]  %+7.0f%%  %s"
              % ("continuous", "M_std", "none", "c-c'", lo_p, hi_p, 100 * pm,
                 "显著" if (lo_p > 0 or hi_p < 0) else "不显著（含0）"))
        rows.append({"mediator": mediator, "x_mode": "continuous", "m_mode": "M_std",
                     "covariates": "none", "effect_def": "c-c'",
                     "nie": nie, "pm": pm, "ci_pct_lo": lo_p, "ci_pct_hi": hi_p,
                     "significant_pct": bool(lo_p > 0 or hi_p < 0)})
        print()

    # ---- 附带复核：PPARG "独立预后" 断言 P=0.019 是否可达 ----
    print("=" * 78)
    print("附带复核：手稿 3.6 称 \"high expression of PPARG was independently "
          "associated with poor survival (P=0.019)\"")
    print("-" * 78)
    pp = df["PPARG"].values
    for cut_name, grp in [
        ("中位数二分", (pp > np.median(pp)).astype(float)),
        ("三分位上 1/3", (pp > np.percentile(pp, 66.667)).astype(float)),
        ("连续(z)", zs(pp)),
    ]:
        for adj_name, C in [("未调整", np.zeros((n, 0))),
                            ("调整评分", zs(df["X"].values)[:, None]),
                            ("调整评分+年龄", np.column_stack([zs(df["X"].values),
                                                        age - np.nanmean(age)]))]:
            Xd = grp[:, None] if cut_name != "连续(z)" else grp[:, None]
            bb, se, _, _ = cox_fit(np.column_stack([Xd, C]) if C.shape[1] else Xd, t, d)
            from scipy import stats as st
            p = 2 * (1 - st.norm.cdf(abs(bb[0] / se[0])))
            print("  %-12s / %-14s HR = %.3f (%.3f-%.3f)  P = %.4f"
                  % (cut_name, adj_name, np.exp(bb[0]),
                     np.exp(bb[0] - 1.96 * se[0]), np.exp(bb[0] + 1.96 * se[0]), p))
    print()
    from s3_mediation import logrank
    print("  中位数二分 log-rank P = %.4f  （数据自带的 "
          "PPARG_CD36_analysis_stats.csv 记录为 chi-sq=2.87, P=0.0901）"
          % logrank(t, d, (pp > np.median(pp)).astype(int)))
    print("  结论：三种切点 × 三种调整方式共 9 种组合中，最小 P 值见上；"
          "手稿的 P=0.019 需另行说明其口径")

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "s3_robustness.csv"),
                              index=False, encoding="utf-8-sig")
    print("\n-> 网格结果已写出 s3_robustness.csv")


if __name__ == "__main__":
    main()
