"""
Vizualizace 1D časově závislé Schrödingerovy rovnice
Crank-Nicolsonova metoda
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# ── Styl ───────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
})

COLORS_TIMES = ["#2196F3", "#4CAF50", "#FF9800", "#E91E63", "#9C27B0", "#00BCD4"]
X_MIN, X_MAX = -25.0, 25.0


# ══════════════════════════════════════════════════════════════════════════════
# Pomocná funkce: načte soubor robustně bez ohledu na počet sloupců
# ══════════════════════════════════════════════════════════════════════════════

def _load_prob_file(fname):
    if not os.path.exists(fname):
        raise FileNotFoundError(fname)

    # Přečti počet sloupců z prvního datového řádku
    with open(fname) as f:
        f.readline()                        # přeskoč hlavičku
        first = f.readline().split()
    ncols = len(first)

    if ncols == 7:
        d = np.genfromtxt(fname, skip_header=1, invalid_raise=False,
                          filling_values=np.nan)
        # Pokud genfromtxt vrátil různý počet sloupců, vyfiltruj špatné řádky
        if d.ndim == 1:
            d = d.reshape(1, -1)
        mask = np.isfinite(d).all(axis=1)
        d = d[mask]
        return dict(t=d[:,0], p=d[:,1], j_plus=d[:,2], j_minus=d[:,3],
                    Jp_cumul=d[:,4], Jm_cumul=d[:,5], conserv=d[:,6])

    elif ncols == 4:
        # Starý formát – kumulativní toky musíme dopočítat integrací
        d = np.genfromtxt(fname, skip_header=1, invalid_raise=False,
                          filling_values=np.nan)
        mask = np.isfinite(d).all(axis=1)
        d = d[mask]
        t      = d[:, 0]
        p      = d[:, 1]
        j_plus = d[:, 2]    # okamžitý tok (ven vpravo)
        j_minus = d[:, 3]   # okamžitý tok (ven vlevo)
        dt = np.diff(t, prepend=t[0])
        Jp_cumul  = np.cumsum(j_plus  * dt)
        Jm_cumul  = np.cumsum(j_minus * dt)
        conserv = p + Jp_cumul + Jm_cumul
        return dict(t=t, p=p, j_plus=j_plus, j_minus=j_minus,
                    Jp_cumul=Jp_cumul, Jm_cumul=Jm_cumul, conserv=conserv)

    else:
        raise ValueError(
            f"{fname}: neočekávaný počet sloupců {ncols}. "
            "Očekáváno 4 nebo 7."
        )


# ══════════════════════════════════════════════════════════════════════════════
# 0. Potenciál a počáteční stav
# ══════════════════════════════════════════════════════════════════════════════

def plot_potential(output="fig0_potential.png"):
    x = np.linspace(-10, 10, 1000)
    V = (x**2 - 1.0) * np.exp(-x**2 / 4.0)
    psi0 = np.sqrt(2.0 / np.pi) * np.exp(-2.0 * (x - 1.0)**2)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, V,    color="red",     lw=2.5, label=r"$V(x) = (x^2-1)e^{-x^2/4}$")
    ax.fill_between(x, 0, V, where=(V > 0), alpha=0.10, color="red")
    ax.fill_between(x, 0, V, where=(V < 0), alpha=0.10, color="blue")
    ax.plot(x, psi0, color="#2196F3", lw=2,   ls="--", label=r"$|\psi(x,0)|^2$")
    ax.axhline(0,  color="gray", lw=0.7, ls=":")
    ax.axvline(-7, color="gray", lw=0.8, ls="--", alpha=0.6, label="hranice $a=7$")
    ax.axvline( 7, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel("$x$"); ax.set_ylabel("$V(x)$")
    ax.set_title("Dvojitý potenciál a počáteční vlnová funkce")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {output}")


# ══════════════════════════════════════════════════════════════════════════════
# 1. Vlnová funkce v několika časech + srovnání s volnou částicí
# ══════════════════════════════════════════════════════════════════════════════

def analytical_free(x, t):
    """
    |ψ(x,t)|² pro volnou Gaussovu vlnovou funkci:
    ψ(x,0) = (2/π)^(1/4) exp[-(x-1)²]
    """
    sigma2 = 1.0 + 4.0 * t**2          # (2σ)² s časem roste
    return np.exp(-2.0 * (x - 1.0)**2 / sigma2) / np.sqrt(np.pi * sigma2 / 2.0)


def plot_wave_functions(times=(0, 1, 2, 3, 4, 5), prefix="wave",
                        free_prefix="free", output="fig1_wave_functions.png"):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    axes = axes.flatten()

    x_pot = np.linspace(X_MIN, X_MAX, 1000)
    V_pot = (x_pot**2 - 1.0) * np.exp(-x_pot**2 / 4.0)

    for ax, t, col in zip(axes, times, COLORS_TIMES):
        # Simulace s potenciálem
        fname = f"{prefix}_t_{int(t * 1000)}.dat"
        if os.path.exists(fname):
            d = np.loadtxt(fname)
            ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=r"$|\psi|^2$ sim.")
        else:
            ax.text(0.5, 0.5, f"{fname}\nnenalezen", transform=ax.transAxes,
                    ha="center", va="center", color="red", fontsize=8)

        # Volná částice (numerická nebo analytická)
        if t == 0:
            x_an = np.linspace(X_MIN, X_MAX, 2000)
            ax.plot(x_an, analytical_free(x_an, t), "k--", lw=1.2,
                    label="volná (anal.)")
        else:
            fname_f = f"{free_prefix}_t_{int(t * 1000)}.dat"
            if os.path.exists(fname_f):
                df = np.loadtxt(fname_f)
                ax.plot(df[:, 0], df[:, 1], "k--", lw=1.2, label="volná (num.)")
            x_an = np.linspace(X_MIN, X_MAX, 2000)
            ax.plot(x_an, analytical_free(x_an, t), color="gray", lw=0.8,
                    ls=":", label="volná (anal.)")

        # Potenciál (škálovaný)
        ax.fill_between(x_pot, 0, V_pot * 0.15, alpha=0.12,
                        color="red", label="V(x) [škál.]")
        ax.axvline(-7, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.axvline( 7, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.set_title(f"$t = {t}$", fontsize=12)
        ax.set_xlim(-15, 15); ax.set_ylim(-0.05, 0.75)
        ax.set_xlabel("$x$"); ax.set_ylabel(r"$|\psi(x,t)|^2$")
        ax.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        r"Vlnová funkce $|\psi(x,t)|^2$ – dvojitá bariéra vs. volná částice",
        fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {output}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Časová závislost p(T), J±(T), zákon zachování
# ══════════════════════════════════════════════════════════════════════════════

def plot_time_evolution(prob_file="wave_prob.dat",
                        output="fig2_time_evolution.png"):
    try:
        r = _load_prob_file(prob_file)
    except FileNotFoundError:
        print(f"  ✗  {prob_file} nenalezen"); return

    t      = r["t"];       p      = r["p"]
    Jp_cumul = r["Jp_cumul"];  Jm_cumul = r["Jm_cumul"]
    j_plus = r["j_plus"];  j_minus= r["j_minus"]
    conserv= r["conserv"]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)

    # Panel 1: kumulativní veličiny a zachování
    ax = axes[0]
    ax.plot(t, p,       color="#2196F3", lw=2,   label=r"$p(T)$  v $[-7,7]$")
    ax.plot(t, Jp_cumul,  color="#E91E63", lw=2,   label=r"$J_+(T)$ kum. (pravá hranice)")
    ax.plot(t, Jm_cumul,  color="#4CAF50", lw=2,   label=r"$J_-(T)$ kum. (levá hranice)")
    ax.plot(t, conserv, "k",             lw=1.2, ls="--",
            label=r"$p + J_+ + J_-$  (zachování)")
    ax.axhline(1.0, color="gray", lw=0.7, ls=":")
    ax.set_ylabel("Pravděpodobnost")
    ax.set_title("Pravděpodobnost a kumulativní toky v čase")
    ax.legend(fontsize=9); ax.set_ylim(-0.02, 1.12)

    # Panel 2: okamžitý tok
    ax = axes[1]
    ax.plot(t, j_plus,  color="#E91E63", lw=1.5, label=r"$j(+7,t)$ okamžitý")
    ax.plot(t, j_minus, color="#4CAF50", lw=1.5, label=r"$j(-7,t)$ okamžitý (ven)")
    ax.axhline(0, color="gray", lw=0.7, ls=":")
    ax.set_ylabel("Tok $j$")
    ax.set_title(r"Okamžitý tok $j$ na hranicích $x = \pm 7$")
    ax.legend(fontsize=9)

    # Panel 3: odchylka od zachování
    ax = axes[2]
    ax.plot(t, (conserv - 1.0) * 1e3, color="darkorange", lw=1.5)
    ax.axhline(0, color="gray", lw=0.7, ls=":")
    ax.set_xlabel("Čas $T$")
    ax.set_ylabel(r"$(p+J_++J_-) - 1\ [\times 10^{-3}]$")
    ax.set_title("Odchylka od zákona zachování")

    fig.suptitle("Časový vývoj pravděpodobností a toků",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {output}")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Konvergenční studie: p, J+, J- pro T=1, a=7
# ══════════════════════════════════════════════════════════════════════════════

def plot_convergence(output="fig3_convergence.png"):
    configs = [
        ("conv_dx100_prob.dat", 0.1),
        ("conv_dx50_prob.dat",  0.05),
        ("conv_dx25_prob.dat",  0.025),
    ]

    results = []
    for fname, dx in configs:
        try:
            r = _load_prob_file(fname)
        except FileNotFoundError:
            print(f"  ✗  {fname} nenalezen, přeskakuji"); continue
        except ValueError as e:
            print(f"  ✗  {fname}: {e}"); continue

        idx = np.argmin(np.abs(r["t"] - 1.0))
        results.append({
            "dx":      dx,
            "p":       r["p"][idx],
            "Jp_cumul":  r["Jp_cumul"][idx],
            "Jm_cumul":  r["Jm_cumul"][idx],
            "conserv": r["conserv"][idx],
        })

    if len(results) < 2:
        print("  ✗  Nedostatek konvergenčních dat (potřeba ≥ 2 sítě)"); return

    dxs  = np.array([r["dx"]      for r in results])
    ps   = np.array([r["p"]       for r in results])
    Jps  = np.array([r["Jp_cumul"]  for r in results])
    Jms  = np.array([r["Jm_cumul"]  for r in results])
    cons = np.array([r["conserv"] for r in results])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: hodnoty veličin vs Δx
    ax = axes[0]
    ax.semilogx(dxs, ps,   "o-", color="#2196F3", ms=8, lw=2, label=r"$p(T=1)$")
    ax.semilogx(dxs, Jps,  "s-", color="#E91E63", ms=8, lw=2, label=r"$J_+(T=1)$")
    ax.semilogx(dxs, Jms,  "^-", color="#4CAF50", ms=8, lw=2, label=r"$J_-(T=1)$")
    ax.semilogx(dxs, cons, "D-", color="gray",    ms=8, lw=1.5,
                label=r"$p+J_++J_-$")
    ax.set_xlabel(r"$\Delta x$"); ax.set_ylabel("Hodnota")
    ax.set_title(r"Konvergence $p,\,J_\pm$ při $T=1,\; a=7$")
    ax.legend(fontsize=9); ax.invert_xaxis()

    # Číselné hodnoty nad body
    for dx, p_val in zip(dxs, ps):
        ax.annotate(f"{p_val:.6f}", (dx, p_val),
                    textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=7, color="#2196F3")

    # Panel 2: relativní konvergence (vůči nejjemnější síti)
    ax = axes[1]
    p_ref  = ps[-1];  Jp_ref = Jps[-1];  Jm_ref = Jms[-1]
    eps    = 1e-15
    rel_p  = np.abs(ps  - p_ref)  / (np.abs(p_ref)  + eps)
    rel_Jp = np.abs(Jps - Jp_ref) / (np.abs(Jp_ref) + eps)
    rel_Jm = np.abs(Jms - Jm_ref) / (np.abs(Jm_ref) + eps)

    coarse = dxs[:-1]
    if len(coarse) > 0:
        ax.loglog(coarse, rel_p[:-1],  "o-", color="#2196F3", ms=8, lw=2, label=r"$p$")
        ax.loglog(coarse, rel_Jp[:-1], "s-", color="#E91E63", ms=8, lw=2, label=r"$J_+$")
        ax.loglog(coarse, rel_Jm[:-1], "^-", color="#4CAF50", ms=8, lw=2, label=r"$J_-$")
        # Referenční O(Δx²) přímka
        ref_line = 5e-2 * (coarse / coarse[0])**2
        ax.loglog(coarse, ref_line, "k--", lw=1, label=r"$O(\Delta x^2)$")

    ax.set_xlabel(r"$\Delta x$"); ax.set_ylabel("Relativní chyba")
    ax.set_title("Řád konvergence (vůči nejjemnější síti)")
    ax.legend(fontsize=9); ax.invert_xaxis()

    fig.suptitle(
        r"Konvergenční studie: $p, J_+, J_-$ při $T=1$ pro různá $\Delta x$",
        fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  {output}")


# ══════════════════════════════════════════════════════════════════════════════
# main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Generuji obrázky…\n")
    plot_potential()
    plot_wave_functions(times=[0, 1, 2, 3, 4, 5])
    plot_time_evolution()
    plot_convergence()
    print("\nHotovo.")
