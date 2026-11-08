# -*- coding: utf-8 -*-
"""
S3 · 中介分析重建（P4 统计与效力复核官）
--------------------------------------------------------------------------
目的：项目内不存在任何中介分析代码（已全量扫描全部 .R/.py 确认）。
      本脚本按手稿 Methods 2.13 所述方法（bootstrap 1000 次重抽样）
      在 TCGA-SARC n=50（23 例死亡事件）上重建该分析，
      并输出可直接引用的效应量、置信区间与可复现代码。

方法：
  a 路径 : OLS   M ~ X              （X 已中心化标准化）
  b 路径 : Cox   Y ~ X + M          （M 系数）
  c 总效应: Cox  Y ~ X
  NIE    : a * b                    （自然间接效应，log-hazard 尺度）
  PM     : a*b / c                  （中介比例）
  区间   : 1000 次个案重抽样的百分位法

数据源：
  data/results/02_feature_scores.csv     （Lipid_Metabolism_Score）
  data/results/03_clinical_outcomes.csv  （OS_time / OS_status）
  data/raw/TCGA_LPS_real_log2cpm.rds     （PPARG / CD36 表达，log2CPM）
"""
import os
import json
import csv
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats

BASE = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1")
RAW = os.path.join(BASE, "data", "raw")
RES = os.path.join(BASE, "data", "results")
OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
os.makedirs(OUT, exist_ok=True)

SEED = 20260917
NBOOT = 1000

ENS = {"PPARG": "ENSG00000132170", "CD36": "ENSG00000135218",
       "FASN": "ENSG00000169710", "ACLY": "ENSG00000131473",
       "SCD": "ENSG00000099194", "SREBF1": "ENSG00000072310",
       "HMGCR": "ENSG00000113161"}


# ---------------------------------------------------------------- 数据准备
def load_data():
    scores = pd.read_csv(os.path.join(RES, "02_feature_scores.csv"), index_col=0)
    clinical = pd.read_csv(os.path.join(RES, "03_clinical_outcomes.csv"), index_col=0)
    l2 = pyreadr.read_r(os.path.join(RAW, "TCGA_LPS_real_log2cpm.rds"))[None]

    # ENSEMBL -> SYMBOL（每个基因在本矩阵中唯一）
    sym = {}
    for g, eid in ENS.items():
        hit = [i for i in l2.index if i.split(".")[0] == eid]
        assert len(hit) == 1, (g, hit)
        sym[g] = hit[0]
    expr = pd.DataFrame({g: l2.loc[r].astype(float) for g, r in sym.items()})

    common = expr.index.intersection(scores.index).intersection(clinical.index)
    expr = expr.loc[common]
    scores = scores.loc[common]
    clinical = clinical.loc[common]

    df = pd.DataFrame({
        "X": scores["Lipid_Metabolism_Score"].astype(float),
        "PPARG": expr["PPARG"],
        "CD36": expr["CD36"],
        "time": pd.to_numeric(clinical["OS_time"], errors="coerce"),
        "event": pd.to_numeric(clinical["OS_status"], errors="coerce"),
        "gender": clinical["Gender"],
        "subtype": clinical["Subtype"],
    })
    df = df.dropna(subset=["time", "event"])
    return df


# ------------------------------------------------- 快速 Cox PH（Breslow）
def cox_fit(X, time, event, max_iter=60, tol=1e-10, ties="efron", beta0=None):
    """Cox 部分似然 + Newton-Raphson。ties: 'efron'（R survival::coxph 与
    lifelines 的默认）或 'breslow'。返回 (beta, se, loglik, converged)."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    t = np.asarray(time, dtype=float)
    d = np.asarray(event, dtype=float)
    n, p = X.shape
    beta = np.zeros(p) if beta0 is None else np.asarray(beta0, dtype=float).copy()

    order = np.argsort(t, kind="stable")        # 时间升序
    Xo, to, do = X[order], t[order], d[order]

    # 预先算出每个唯一事件时间的（风险集起点、事件位置）
    uniq = np.unique(to[do > 0])
    noties = (int(do.sum()) == uniq.size)      # 无并列时间 -> 可全向量化
    ev_all = np.where(do > 0)[0]
    ev_groups = []
    for u in uniq:
        s = int(np.searchsorted(to, u, side="left"))     # 风险集 = s..n-1
        ev_pos = np.where((to == u) & (do > 0))[0]
        ev_groups.append((s, ev_pos))

    def obj(b):
        eta = Xo @ b
        # 注意：不可减去 eta.max()——该常数随 beta 变化，会破坏偏似然
        # （协变量已标准化，n 小，exp 不会溢出）
        w = np.exp(eta)
        # 后缀和：S[s] = sum_{j>=s}
        S0 = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])
        S1 = np.concatenate([np.cumsum((w[:, None] * Xo)[::-1], axis=0)[::-1],
                             np.zeros((1, p))], axis=0)
        S2 = np.concatenate(
            [np.cumsum((w[:, None, None] * Xo[:, :, None] * Xo[:, None, :])[::-1],
                       axis=0)[::-1], np.zeros((1, p, p))], axis=0)

        if noties:
            # 全向量化：事件位置即风险集起点
            den = S0[ev_all]
            num1 = S1[ev_all]
            num2 = S2[ev_all]
            ll = eta[ev_all].sum() - np.log(den).sum()
            mu = num1 / den[:, None]
            gr = Xo[ev_all].sum(axis=0) - mu.sum(axis=0)
            hs = (num2 / den[:, None, None]
                  - mu[:, :, None] * mu[:, None, :]).sum(axis=0)
            return ll, gr, hs

        ll = 0.0
        gr = np.zeros(p)
        hs = np.zeros((p, p))
        for s, ev_pos in ev_groups:
            dd = ev_pos.size
            ee = w[ev_pos]
            e0 = ee.sum()
            e1 = (ee[:, None] * Xo[ev_pos]).sum(axis=0)
            e2 = (ee[:, None, None] * Xo[ev_pos][:, :, None]
                  * Xo[ev_pos][:, None, :]).sum(axis=0)
            ll += eta[ev_pos].sum()
            gr += Xo[ev_pos].sum(axis=0)
            if ties == "efron":
                fracs = np.arange(dd) / dd
            else:
                fracs = np.zeros(dd)
            for fr in fracs:
                den = S0[s] - fr * e0
                num1 = S1[s] - fr * e1
                num2 = S2[s] - fr * e2
                ll -= np.log(den)
                mu = num1 / den
                gr -= mu
                hs += num2 / den - np.outer(mu, mu)
        return ll, gr, hs

    ll_old = None
    for _ in range(max_iter):
        ll, gr, hs = obj(beta)
        try:
            step = np.linalg.solve(hs + 1e-10 * np.eye(p), gr)
        except np.linalg.LinAlgError:
            break
        # 步长回溯，保证似然不下降
        s = 1.0
        for _ in range(30):
            nb = beta + s * step
            nll = obj(nb)[0]
            if np.isfinite(nll) and nll >= ll - 1e-12:
                break
            s /= 2
        beta = beta + s * step
        if ll_old is not None and abs(ll - ll_old) < tol * (abs(ll) + 1):
            ll_old = ll
            break
        ll_old = ll

    ll, gr, hs = obj(beta)
    try:
        cov = np.linalg.inv(hs)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    return beta, se, ll, np.all(np.isfinite(beta))


def ols_fit(y, X):
    """y ~ X（无截距列需自行加）。返回 beta, se。"""
    X = np.column_stack([np.ones(len(y)), X])
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    resid = y - X @ b
    dof = max(len(y) - X.shape[1], 1)
    s2 = (resid @ resid) / dof
    se = np.sqrt(np.diag(XtX_inv) * s2)
    return b, se


# ------------------------------------------------------------ 中介效应
def mediation_once(X, M, time, event):
    """返回 (a, b, c, cprime, nie, pm)。X/M 需已标准化。"""
    # a 路径
    ab, _ = ols_fit(M, X[:, None])
    a = ab[1]
    # b 路径（含 X 的 Cox）
    bb, _, _, _ = cox_fit(np.column_stack([X, M]), time, event)
    b = bb[1]
    cprime = bb[0]
    # 总效应
    cb, _, _, _ = cox_fit(X[:, None], time, event)
    c = cb[0]
    nie = a * b
    pm = nie / c if c != 0 else np.nan
    return a, b, c, cprime, nie, pm


def mediation_bootstrap(X, M, time, event, nboot=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(X)
    boot = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        try:
            a, b, c, cp, nie, pm = mediation_once(X[idx], M[idx], time[idx], event[idx])
            if np.all(np.isfinite([a, b, c, nie])):
                boot.append((a, b, c, cp, nie, pm))
        except Exception:
            continue
    boot = np.array(boot)
    return boot


def pct_ci(v, lo=2.5, hi=97.5):
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (np.nan, np.nan)
    return (float(np.percentile(v, lo)), float(np.percentile(v, hi)))


# ------------------------------------------------------------------ 主流程
def main():
    df = load_data()
    n = len(df)
    events = int(df["event"].sum())
    print("=" * 74)
    print("【0】样本与事件")
    print(f"  TCGA-SARC  n = {n}   死亡事件 = {events}  ({events/n:.1%})")
    print(f"  随访时间 min/median/max = {df['time'].min():.0f} / "
          f"{df['time'].median():.0f} / {df['time'].max():.0f} (天)")
    print(f"  性别: {df['gender'].value_counts().to_dict()}")

    # 标准化（z-score）
    Xz = (df["X"] - df["X"].mean()) / df["X"].std(ddof=1)
    t = df["time"].values
    d = df["event"].values

    report = {"n": n, "events": events, "descriptives": {}, "mediators": {}}

    # ---------- 描述统计与相关（a 路径前提）----------
    print()
    print("=" * 74)
    print("【1】评分-表达相关（核对正文 R=0.62 / 0.58）")
    for g in ("PPARG", "CD36"):
        m = df[g].values
        pear = stats.pearsonr(Xz.values, m)
        spear = stats.spearmanr(Xz.values, m)
        r2 = pear.statistic ** 2
        print(f"  {g:6s}  Pearson r = {pear.statistic:.3f}  (P = {pear.pvalue:.3g})"
              f"   Spearman rho = {spear.statistic:.3f}  (P = {spear.pvalue:.3g})"
              f"   R^2 = {r2:.3f}")
        report["descriptives"][g] = {
            "pearson_r": float(pear.statistic), "pearson_p": float(pear.pvalue),
            "spearman_rho": float(spear.statistic), "spearman_p": float(spear.pvalue),
            "r_squared": float(r2),
        }
    # a 路径线性模型
    print()
    for g in ("PPARG", "CD36"):
        m = df[g].values
        b, se = ols_fit(m, Xz.values[:, None])
        # 标准化斜率 = b[1]*sd(X)/sd(M)；X 已 z 化，故为 b[1]/sd(M)
        a_std = b[1] / m.std(ddof=1)
        print(f"  a 路径 {g:6s}: 未标准化斜率 = {b[1]:.4f} (SE {se[1]:.4f}, "
              f"P = {2*(1-stats.t.cdf(abs(b[1]/se[1]), n-2)):.3g})"
              f"   标准化斜率 a = {a_std:.4f}")
        report["descriptives"][g]["a_path_raw"] = float(b[1])
        report["descriptives"][g]["a_path_std"] = float(a_std)

    # ---------- 高/低组表达与生存 ----------
    print()
    print("=" * 74)
    print("【2】高分组 vs 低分组（核对正文 P=0.004 / 0.006 与 PPARG 生存 P=0.019）")
    med = df["X"].median()
    grp = np.where(df["X"] > med, "High", "Low")
    for g in ("PPARG", "CD36"):
        hi = df[g].values[grp == "High"]
        lo = df[g].values[grp == "Low"]
        u = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        tt = stats.ttest_ind(hi, lo, equal_var=False)
        print(f"  {g:6s} High = {hi.mean():.3f}  Low = {lo.mean():.3f}"
              f"   Wilcoxon P = {u.pvalue:.3g}   Welch t P = {tt.pvalue:.3g}")
        report["descriptives"][g].update(
            {"high_mean": float(hi.mean()), "low_mean": float(lo.mean()),
             "wilcoxon_p": float(u.pvalue), "ttest_p": float(tt.pvalue)})

    # PPARG 高低表达的生存（log-rank 手算 + 未调整 Cox）
    for g in ("PPARG", "CD36"):
        m = df[g].values
        gmed = np.median(m)
        gg = (m > gmed).astype(int)
        # log-rank
        lr = logrank(t, d, gg)
        bx, sx, _, _ = cox_fit(gg[:, None].astype(float), t, d)
        hr = np.exp(bx[0])
        lo_ci, hi_ci = np.exp(bx[0] - 1.96 * sx[0]), np.exp(bx[0] + 1.96 * sx[0])
        # 调整评分后的 Cox
        badj, sadj, _, _ = cox_fit(np.column_stack([gg.astype(float), Xz.values]), t, d)
        p_adj = 2 * (1 - stats.norm.cdf(abs(badj[0] / sadj[0])))
        print(f"  {g:6s} 高/低表达  log-rank P = {lr:.4g} | "
              f"未调整 Cox HR = {hr:.3f} ({lo_ci:.3f}-{hi_ci:.3f}) | "
              f"调整评分后 HR = {np.exp(badj[0]):.3f}, P = {p_adj:.4f}")
        report["descriptives"][g].update(
            {"logrank_p_highlow": float(lr), "cox_hr_unadj": float(hr),
             "cox_hr_unadj_lo": float(lo_ci), "cox_hr_unadj_hi": float(hi_ci),
             "cox_hr_adj_score": float(np.exp(badj[0])), "cox_p_adj_score": float(p_adj)})

    # ---------- 中介分析 ----------
    print()
    print("=" * 74)
    print("【3】中介分析重建（bootstrap 1000 次）")
    for g in ("PPARG", "CD36"):
        Mraw = df[g].values
        Mz = (Mraw - Mraw.mean()) / Mraw.std(ddof=1)
        a, b, c, cp, nie, pm = mediation_once(Xz.values, Mz, t, d)
        boot = mediation_bootstrap(Xz.values, Mz, t, d)
        nie_lo, nie_hi = pct_ci(boot[:, 4])
        pm_lo, pm_hi = pct_ci(boot[:, 5])
        a_lo, a_hi = pct_ci(boot[:, 0])
        b_lo, b_hi = pct_ci(boot[:, 1])
        c_lo, c_hi = pct_ci(boot[:, 2])
        cp_lo, cp_hi = pct_ci(boot[:, 3])
        print(f"  --- 中介变量 {g} (n={n}, 事件={events}) ---")
        print(f"    a (X->M, 标准化)      = {a:.4f}   [{a_lo:.4f}, {a_hi:.4f}]")
        print(f"    b (M->Y | X, log-HR)  = {b:.4f}   [{b_lo:.4f}, {b_hi:.4f}]"
              f"   HR_M = {np.exp(b):.3f}")
        print(f"    c   总效应           = {c:.4f}   [{c_lo:.4f}, {c_hi:.4f}]"
              f"   HR_total = {np.exp(c):.3f}")
        print(f"    c'  直接效应         = {cp:.4f}   [{cp_lo:.4f}, {cp_hi:.4f}]")
        print(f"    NIE = a*b             = {nie:.4f}   [{nie_lo:.4f}, {nie_hi:.4f}]"
              f"   {'含 0（不显著）' if nie_lo <= 0 <= nie_hi else '不含 0'}")
        print(f"    PM  中介比例          = {pm:.3f}   [{pm_lo:.3f}, {pm_hi:.3f}]")
        report["mediators"][g] = {
            "a": float(a), "a_ci": [a_lo, a_hi],
            "b": float(b), "b_ci": [b_lo, b_hi], "b_hr": float(np.exp(b)),
            "c": float(c), "c_ci": [c_lo, c_hi], "c_hr": float(np.exp(c)),
            "cprime": float(cp), "cprime_ci": [cp_lo, cp_hi],
            "nie": float(nie), "nie_ci": [nie_lo, nie_hi],
            "nie_crosses_zero": bool(nie_lo <= 0 <= nie_hi),
            "pm": float(pm), "pm_ci": [pm_lo, pm_hi],
            "boot_valid": int(boot.shape[0]),
        }

    # ---------- 反事实检验：CI 上限意味着什么 ----------
    print()
    print("=" * 74)
    print("【4】阴性结论的可排除范围（关键：CI 上限）")
    for g in ("PPARG", "CD36"):
        r = report["mediators"][g]
        hi = r["nie_ci"][1]
        c = r["c"]
        print(f"  {g:6s} NIE 95%CI 上限 = {hi:.4f}  ->  在该上限下中介比例 "
              f"PM_max = {hi/c:.1%}")
        # 反推：若真实 NIE 为此上限，间接路径对应的 HR
        print(f"          对应的间接路径等效 HR = {np.exp(hi):.3f}"
              f"   (即数据无法排除“中介使风险上升 {np.exp(hi)-1:.1%}”以下的中介)")
        r["pm_max_from_ci"] = float(hi / c)

    with open(os.path.join(OUT, "s3_mediation_results.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print()
    print("-> 结果已写出 s3_mediation_results.json")
    return report


def logrank(time, event, group):
    """两组 log-rank 检验。"""
    time = np.asarray(time, float)
    event = np.asarray(event, float)
    group = np.asarray(group)
    times = np.unique(time[event == 1])
    O1 = E1 = V = 0.0
    for tu in times:
        at_risk = time >= tu
        n1 = np.sum(at_risk & (group == 1))
        n0 = np.sum(at_risk & (group == 0))
        n = n1 + n0
        d1 = np.sum((time == tu) & (event == 1) & (group == 1))
        d = np.sum((time == tu) & (event == 1))
        if n <= 1:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan
    chi2 = (O1 - E1) ** 2 / V
    return float(1 - stats.chi2.cdf(chi2, 1))


if __name__ == "__main__":
    main()
