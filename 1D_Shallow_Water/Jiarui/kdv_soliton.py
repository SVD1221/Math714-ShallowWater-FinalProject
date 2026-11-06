"""
Numerical Investigation of Soliton Dynamics in the
Korteweg-de Vries (KdV) Equation

This script investigates numerical solutions of the KdV equation

    u_t + 6*u*u_x + u_xxx = 0

on a periodic spatial domain using a finite-difference
discretization and the classical fourth-order Runge-Kutta
(RK4) method for time integration.

The numerical experiments focus on three aspects:

1. Convergence Analysis
   Evaluate numerical accuracy by comparing the computed
   single-soliton solution with the analytical solution.

2. Numerical Stability
   Investigate the timestep restrictions of the explicit
   RK4 method and compare its stability with forward Euler.

3. Soliton Interaction
   Simulate the interaction of two solitons and examine
   their phase shifts in comparison with theoretical predictions.

"""

import json
import os
import time
 
import numpy as np
import matplotlib
import matplotlib.ticker
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
 
OUT = os.path.dirname(os.path.abspath(__file__))
R_LINEAR = 4.0 * np.sqrt(2.0) / (3.0 * np.sqrt(3.0))  # ~1.089
 
# Colors (reference categorical palette, light mode) and ink
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
 
plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.5, "legend.frameon": False,
    "savefig.dpi": 200, "savefig.bbox": "tight",
})
 
 
# ----------------------------------------------------------------------------
# Core numerics
# ----------------------------------------------------------------------------
def make_grid(L, dx):
    N = int(round(L / dx))
    x = -L / 2 + dx * np.arange(N)
    return x, N
 
 
def wrap(xi, L):
    return (xi + L / 2) % L - L / 2
 
 
def soliton(x, t, c, x0, L):
    xi = wrap(x - x0 - c * t, L)
    return 0.5 * c / np.cosh(0.5 * np.sqrt(c) * xi) ** 2
 
 
def rhs(u, dx):
    p = np.concatenate((u[-2:], u, u[:2]))
    um2, um1, u0, up1, up2 = p[:-4], p[1:-3], p[2:-2], p[3:-1], p[4:]
    nonlin = (up1 + u0 + um1) * (up1 - um1) / dx          # = 6 u u_x (ZK form)
    disp = (up2 - 2.0 * up1 + 2.0 * um1 - um2) / (2.0 * dx ** 3)
    return -nonlin - disp
 
 
def rk4_step(u, dt, dx):
    k1 = rhs(u, dx)
    k2 = rhs(u + 0.5 * dt * k1, dx)
    k3 = rhs(u + 0.5 * dt * k2, dx)
    k4 = rhs(u + dt * k3, dx)
    return u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
 
 
def euler_step(u, dt, dx):
    return u + dt * rhs(u, dx)
 
 
def integrate(u, dx, T, r, stepper=rk4_step, snap_every=None):
    """Integrate to time T with dt <= r dx^3 (adjusted to land exactly on T)."""
    n = int(np.ceil(T / (r * dx ** 3)))
    dt = T / n
    snaps, times = [], []
    if snap_every is not None:
        k_snap = max(1, int(round(snap_every / dt)))
        snaps.append(u.copy()); times.append(0.0)
    for k in range(1, n + 1):
        u = stepper(u, dt, dx)
        if snap_every is not None and k % k_snap == 0:
            snaps.append(u.copy()); times.append(k * dt)
    return u, dt, n, (np.array(times), np.array(snaps))
 
 
def invariants(u, dx):
    return dx * u.sum(), dx * (u ** 2).sum()
 
 
def peaks(u, x, dx, L, thresh):
    """Local maxima above thresh, refined by parabolic interpolation."""
    um1, up1 = np.roll(u, 1), np.roll(u, -1)
    idx = np.where((u > um1) & (u >= up1) & (u > thresh))[0]
    out = []
    for j in idx:
        a, b, c = u[j - 1], u[j], u[(j + 1) % len(u)]
        denom = a - 2 * b + c
        s = 0.5 * (a - c) / denom if denom != 0 else 0.0
        out.append((wrap(x[j] + s * dx, L), b - 0.25 * (a - c) * s))
    return sorted(out, key=lambda p: -p[1])  # tallest first
 
 
# ----------------------------------------------------------------------------
# Study 1: accuracy and convergence of a single soliton
# ----------------------------------------------------------------------------
def study_convergence():
    L, c, x0, T, r = 40.0, 4.0, -10.0, 5.0, 0.4
    rows, profiles = [], {}
    for dx in [0.2, 0.1, 0.05, 0.025]:
        x, N = make_grid(L, dx)
        u0 = soliton(x, 0.0, c, x0, L)
        M0, E0 = invariants(u0, dx)
        t0 = time.time()
        u, dt, n, _ = integrate(u0.copy(), dx, T, r)
        wall = time.time() - t0
        ue = soliton(x, T, c, x0, L)
        err = u - ue
        M1, E1 = invariants(u, dx)
        (xp, amp), = peaks(u, x, dx, L, 0.5 * c / 2)[:1]
        rows.append(dict(
            dx=dx, N=N, dt=dt, steps=n, wall_s=round(wall, 2),
            L2=float(np.sqrt(dx * np.sum(err ** 2))),
            Linf=float(np.max(np.abs(err))),
            peak_pos_err=float(wrap(xp - (x0 + c * T), L)),
            amp_err=float(amp - c / 2),
            mass_drift=float(abs(M1 - M0) / abs(M0)),
            energy_drift=float(abs(E1 - E0) / abs(E0)),
        ))
        profiles[dx] = (x, u, ue)
        print(f"[conv] dx={dx:<6} N={N:<5} steps={n:<7} L2={rows[-1]['L2']:.3e} "
              f"pos_err={rows[-1]['peak_pos_err']:+.3e} ({wall:.1f}s)")
    for i in range(1, len(rows)):
        rows[i]["order_L2"] = float(np.log(rows[i - 1]["L2"] / rows[i]["L2"]) / np.log(2))
 
    # Figure 1: profiles and pointwise error
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
    x, _, ue = profiles[0.2]
    ax[0].plot(profiles[0.2][0], profiles[0.2][1], color=C1, label=r"$\Delta x=0.2$")
    ax[0].plot(profiles[0.05][0], profiles[0.05][1], color=C2, label=r"$\Delta x=0.05$")
    ax[0].plot(x, ue, color=INK, ls="--", lw=1.0, label="exact", zorder=3)
    ax[0].set_xlim(5, 15); ax[0].set_xlabel("$x$"); ax[0].set_ylabel("$u(x,5)$")
    ax[0].set_title("(a) Soliton at $t=5$ (zoom)", loc="left")
    ax[0].legend(loc="upper right")
    for dx, col in [(0.2, C1), (0.1, C2), (0.05, C3)]:
        xx, uu, ee = profiles[dx]
        ax[1].plot(xx, uu - ee, color=col, label=rf"$\Delta x={dx}$")
    ax[1].set_xlim(5, 15); ax[1].set_xlabel("$x$"); ax[1].set_ylabel("$u - u_{\\mathrm{exact}}$")
    ax[1].set_title("(b) Pointwise error at $t=5$", loc="left")
    ax[1].legend(loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_profile_error.pdf")); plt.close(fig)
 
    # Figure 2: convergence
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    dxs = np.array([rw["dx"] for rw in rows]); e2 = np.array([rw["L2"] for rw in rows])
    ax.loglog(dxs, e2, "o-", color=C1, ms=5, label=r"$L^2$ error")
    ref = 0.5 * e2[0] * (dxs / dxs[0]) ** 2
    ax.loglog(dxs, ref, color=INK2, ls=":", lw=1.2, label=r"reference slope 2")
    ax.set_xticks(dxs); ax.set_xticklabels([f"{d:g}" for d in dxs])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$\Delta x$"); ax.set_ylabel(r"$\|u-u_{\mathrm{exact}}\|_{2}$ at $t=5$")
    ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_convergence.pdf")); plt.close(fig)
    return rows
 
 
# ----------------------------------------------------------------------------
# Study 2: stiffness and the time-step limit
# ----------------------------------------------------------------------------
def blows_up(dx, r, stepper, nsteps=20000, L=40.0, c=4.0):
    x, _ = make_grid(L, dx)
    u = soliton(x, 0.0, c, -10.0, L)
    dt = r * dx ** 3
    for k in range(nsteps):
        u = stepper(u, dt, dx)
        if k % 200 == 0 and (not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 10.0):
            return True, k
    return (not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 10.0), nsteps
 
 
def study_stability():
    rows = []
    for dx in [0.2, 0.1, 0.05]:
        lo, hi = 0.5, 2.0
        for _ in range(12):  # bisection on r = dt/dx^3
            mid = 0.5 * (lo + hi)
            if blows_up(dx, mid, rk4_step)[0]:
                hi = mid
            else:
                lo = mid
        rows.append(dict(dx=dx, r_crit=round(0.5 * (lo + hi), 4),
                         r_linear=round(float(R_LINEAR), 4)))
        print(f"[stab] dx={dx}: empirical r* = {rows[-1]['r_crit']}  (linear theory {R_LINEAR:.4f})")
    fe = []
    for r in [0.1, 0.5]:
        bu, k = blows_up(0.1, r, euler_step)
        fe.append(dict(dx=0.1, r=r, blew_up=bool(bu), step=int(k)))
        print(f"[stab] forward Euler dx=0.1 r={r}: blew_up={bu} at step {k}")
    return rows, fe
 
 
# ----------------------------------------------------------------------------
# Study 3: two-soliton collision and phase shift
# ----------------------------------------------------------------------------
def study_collision():
    L, dx, r, T = 80.0, 0.1, 0.4, 12.0
    k1, k2 = 1.0, 0.6
    c1, c2 = 4 * k1 ** 2, 4 * k2 ** 2
    x1, x2 = -30.0, -15.0
    x, N = make_grid(L, dx)
 
    u0 = soliton(x, 0, c1, x1, L) + soliton(x, 0, c2, x2, L)
    M0, E0 = invariants(u0, dx)
    u, dt, n, (tt, snaps) = integrate(u0.copy(), dx, T, r, snap_every=0.1)
    M1, E1 = invariants(u, dx)
    pk = peaks(u, x, dx, L, 0.2)[:2]
    (xa, ampa), (xb, ampb) = pk
 
    # single-soliton reference runs on the same grid (cancels numerical speed error)
    ua, *_ = integrate(soliton(x, 0, c1, x1, L), dx, T, r)
    ub, *_ = integrate(soliton(x, 0, c2, x2, L), dx, T, r)
    xa_ref = peaks(ua, x, dx, L, 0.2)[0][0]
    xb_ref = peaks(ub, x, dx, L, 0.2)[0][0]
 
    th = np.log((k1 + k2) / (k1 - k2))
    res = dict(
        kappa1=k1, kappa2=k2, c1=c1, c2=c2, dx=dx, dt=dt, steps=n,
        shift_large_measured=float(wrap(xa - xa_ref, L)),
        shift_small_measured=float(wrap(xb - xb_ref, L)),
        shift_large_vs_exact_free=float(wrap(xa - (x1 + c1 * T), L)),
        shift_small_vs_exact_free=float(wrap(xb - (x2 + c2 * T), L)),
        shift_large_theory=float(th / k1),
        shift_small_theory=float(-th / k2),
        amp_large_after=float(ampa), amp_small_after=float(ampb),
        mass_drift=float(abs(M1 - M0) / abs(M0)),
        energy_drift=float(abs(E1 - E0) / abs(E0)),
    )
    print("[coll]", json.dumps(res, indent=1))
 
    # Figure 3: space-time plot and snapshots
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
    im = ax[0].imshow(snaps, origin="lower", aspect="auto", cmap="Blues",
                      extent=[x[0], x[-1], tt[0], tt[-1]], vmin=0, vmax=c1 / 2)
    ax[0].set_xlim(-35, 25); ax[0].grid(False)
    ax[0].set_xlabel("$x$"); ax[0].set_ylabel("$t$")
    ax[0].set_title("(a) $u(x,t)$, space-time", loc="left")
    cb = fig.colorbar(im, ax=ax[0], pad=0.02); cb.outline.set_visible(False)
    cb.set_label("$u$")
    cb.ax.tick_params(labelsize=7)
    t_coll = (x2 - x1) / (c1 - c2)
    for tsel, col in [(0.0, C1), (round(t_coll, 1), C2), (T, C3)]:
        j = int(np.argmin(np.abs(tt - tsel)))
        ax[1].plot(x, snaps[j], color=col, label=f"$t={tt[j]:.1f}$")
    ax[1].set_xlim(-35, 25); ax[1].set_xlabel("$x$"); ax[1].set_ylabel("$u$")
    ax[1].set_title("(b) Snapshots", loc="left")
    ax[1].legend(loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_collision.pdf")); plt.close(fig)
    return res
 
 
if __name__ == "__main__":
    results = {}
    results["convergence"] = study_convergence()
    results["stability_rk4"], results["stability_euler"] = study_stability()
    results["collision"] = study_collision()
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("done")
    
