# -*- coding: utf-8 -*-
"""
S3 · 效力分析：n=50（23 例死亡事件）下可检出的最小自然间接效应（MDIE）
--------------------------------------------------------------------------
回答的问题：手稿 3.6 断言 "neither of these genes was a significant mediator"。
本脚本量化该阴性结论的可排除范围——即在 n=50 / 23 事件的条件下，
bootstrap 中介分析实际能够检出多大的间接效应。

设计：
  · 暴露 X（脂质评分，z 化）与随访/删失结构取自真实 TCGA-SARC 数据
  · a 路径固定为实测值（评分->介质的标准化斜率），残余方差按实测 R^2 标定
  · 在真实 NIE 值网格上模拟：M = a*X + sqrt(1-a^2)*N(0,1)
                                h(t) = lam * exp(c'*X + b*M),  b = NIE/a
    c' 取 c_obs - NIE（对数风险尺度可加近似）
  · 删失时间由实测删失个体随访时间重抽样得到，lam 二分校准使事件率≈46%
  · 每个重复运行与手稿相同的 bootstrap（百分位法）检验，记录 CI 是否含 0
  · power = CI 排除 0 的比例；80% 效力对应的 NIE 即 MDIE

用法：python s3_power.py PPARG  (或 CD36)
"""
import sys
import os
import json
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3_mediation import load_data, ols_fit, cox_fit, pct_ci

OUT = os.path.join(os.environ.get("LPS_WORK_ROOT", "."), "_work")
SEED = 20260917

GRID = [0.2, 0.4, 0.6, 0.8, 1.0, 1.3]      # 真实 NIE（log-hazard 尺度）
NREP = 300                                  # 每个网格点模拟重复数
NBOOT = 200                                 # 每次重复内的 bootstrap 次数


# ------------------------------------------------------------------ 工具
def mediation_point(X, M, t, d, beta0=None):
    ab, _ = ols_fit(M, X[:, None])
    a = ab[1]
    bb, _, _, _ = cox_fit(np.column_stack([X, M]), t, d, beta0=beta0)
    b, cprime = bb[1], bb[0]
    cb, _, _, _ = cox_fit(X[:, None], t, d)
    c = cb[0]
    return a, b, c, cprime


def boot_detect(X, M, t, d, nboot, rng, beta0=None):
    """返回 (nie_hat, 是否检出, c_hat)。"""
    a, b, c, cp = mediation_point(X, M, t, d, beta0)
    nie = a * b
    n = len(X)
    vals = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        try:
            aa, bb_, cc, _ = mediation_point(X[idx], M[idx], t[idx], d[idx], beta0)
            if np.isfinite(aa) and np.isfinite(bb_):
                vals.append(aa * bb_)
        except Exception:
            continue
    vals = np.asarray(vals)
    if vals.size < 10:
        return nie, False, c
    lo, hi = pct_ci(vals)
    return nie, bool(lo > 0 or hi < 0), c


def simulate_once(args):
    """单次模拟重复：返回 (是否检出, 实际总效应 c_hat)。"""
    (seed, X, a_obs, c_obs, lam, cens_pool, nie_target, nboot) = args
    rng = np.random.default_rng(seed)
    n = len(X)
    b = nie_target / a_obs
    cp = c_obs - nie_target
    # 介质
    resid_sd = np.sqrt(max(1.0 - a_obs ** 2, 1e-6))
    M = a_obs * X + resid_sd * rng.standard_normal(n)
    # 生存
    eta = cp * X + b * M
    T = rng.exponential(1.0 / (lam * np.exp(eta)))
    C = rng.choice(cens_pool, size=n, replace=True)
    Tobs = np.minimum(T, C)
    ev = (T <= C).astype(float)
    if ev.sum() < 3:
        return (False, np.nan)
    _, det, c_hat = boot_detect(X, M, Tobs, ev, nboot, rng)
    return (det, c_hat)


def calibrate_lambda(X, a_obs, c_obs, nie_target, cens_pool, target=0.46, seed=0):
    """二分校准基线风险，使模拟事件率 ≈ target。"""
    rng = np.random.default_rng(seed)
    n = len(X)
    b = nie_target / a_obs
    cp = c_obs - nie_target
    resid_sd = np.sqrt(max(1.0 - a_obs ** 2, 1e-6))
    reps = []
    for _ in range(300):
        M = a_obs * X + resid_sd * rng.standard_normal(n)
        C = rng.choice(cens_pool, size=n, replace=True)
        reps.append((cp * X + b * M, C))

    def frac(lam):
        tot = 0.0
        for eta, C in reps:
            T = rng.exponential(1.0 / (lam * np.exp(eta)))
            tot += (T <= C).mean()
        return tot / len(reps)

    lo, hi = 1e-6, 10.0
    for _ in range(30):
        mid = np.sqrt(lo * hi)
        if frac(mid) > target:
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


# ------------------------------------------------------------------ 主流程
def main(mediator):
    df = load_data()
    t = df["time"].values.astype(float)
    d = df["event"].values.astype(float)
    X = ((df["X"] - df["X"].mean()) / df["X"].std(ddof=1)).values
    Mraw = df[mediator].values
    M = (Mraw - Mraw.mean()) / Mraw.std(ddof=1)
    n = len(X)

    a_obs, b_obs, c_obs, cp_obs = mediation_point(X, M, t, d)
    cens_pool = t[d == 0]
    print("=" * 74)
    print("效力分析 · 介质 = %s" % mediator)
    print("  n=%d  事件=%d  实测 a=%.4f  b=%.4f  c=%.4f  NIE=%.4f"
          % (n, int(d.sum()), a_obs, b_obs, c_obs, a_obs * b_obs))
    print("  删失个体随访池: %d 例, 中位 %.0f 天" % (len(cens_pool), np.median(cens_pool)))

    results = []
    for nie in GRID:
        t0 = time.time()
        lam = calibrate_lambda(X, a_obs, c_obs, nie, cens_pool)
        jobs = [(SEED + 1000 * int(nie * 10) + i, X, a_obs, c_obs, lam,
                 cens_pool, nie, NBOOT) for i in range(NREP)]
        nw = int(os.environ.get("S3_WORKERS", "4"))
        with ProcessPoolExecutor(max_workers=nw) as ex:
            out = list(ex.map(simulate_once, jobs, chunksize=5))
        det = np.array([o[0] for o in out], dtype=bool)
        cvals = np.array([o[1] for o in out], dtype=float)
        cvals = cvals[np.isfinite(cvals)]
        power = float(det.mean())
        rec = {
            "nie_true": nie,
            "b_true": nie / a_obs,
            "cprime_true": c_obs - nie,
            "lambda": lam,
            "power": power,
            "n_rep": int(det.size),
            "c_hat_mean": float(np.mean(cvals)) if cvals.size else np.nan,
            "pm_true_mean": float(nie / np.mean(cvals)) if cvals.size else np.nan,
            "hr_indirect": float(np.exp(nie)),
            "seconds": round(time.time() - t0, 1),
        }
        results.append(rec)
        print("  NIE_true=%.2f (b=%.3f, HR=%.2f)  power=%.3f  "
              "实测 c_hat=%.3f  PM=%.0f%%  [%.0fs]"
              % (nie, nie / a_obs, np.exp(nie), power,
                 rec["c_hat_mean"], 100 * rec["pm_true_mean"], rec["seconds"]))

    # 线性插值求 80% 效力对应 NIE
    ps = np.array([r["power"] for r in results])
    ns = np.array([r["nie_true"] for r in results])
    if ps.min() >= 0.8:
        mdie = float(ns[np.argmin(ps)])
    elif ps.max() < 0.8:
        mdie = None
    else:
        mdie = float(np.interp(0.8, ps, ns))

    summary = {
        "mediator": mediator,
        "n": n, "events": int(d.sum()),
        "observed": {"a": a_obs, "b": b_obs, "c": c_obs, "cprime": cp_obs,
                     "nie": a_obs * b_obs,
                     "pm": a_obs * b_obs / c_obs},
        "grid": results,
        "mdie_80": mdie,
        "mdie_80_pm": (mdie / c_obs) if mdie else None,
        "mdie_80_hr": float(np.exp(mdie)) if mdie else None,
    }
    with open(os.path.join(OUT, "s3_power_%s.json" % mediator), "w",
              encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print()
    if mdie:
        print("  >>> 80%% 效力对应的最小可检出间接效应 NIE = %.3f"
              "  (间接路径 HR = %.2f, 中介比例 PM = %.0f%%)"
              % (mdie, np.exp(mdie), 100 * mdie / c_obs))
    else:
        print("  >>> 全部网格点效力 < 80%%，需外推")
    return summary


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "PPARG")
