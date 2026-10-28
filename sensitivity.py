"""Sensitivity studies for the NGL recovery model. Writes CSVs and PNG charts to results/."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ngl_plant import simulate, Inputs

os.makedirs("results", exist_ok=True)
C2, C3, INK, MUTED, GRID = "#3B5FD0", "#F14A00", "#12131C", "#5B5F72", "#E1E3EC"

def sweep(param, values):
    rows = []
    for v in values:
        r = simulate(Inputs(**{param: float(v)}))
        rows.append(dict(x=float(v), C2=r["recovery_pct"]["C2"], C3=r["recovery_pct"]["C3"],
                         net_kW=r["net_compression_kW"], ngl_tpd=r["ngl_tpd"],
                         T_exp=r["T_expander_out_C"]))
    with open(f"results/sweep_{param}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    return rows

def style(ax, xlabel, ylabel):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=MUTED, fontsize=10); ax.set_ylabel(ylabel, color=MUTED, fontsize=10)

def panel(rows, xlabel, title, fname, base_x, invert=False):
    x = [r["x"] for r in rows]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=160, gridspec_kw=dict(wspace=0.35))
    for key, col in (("C2", C2), ("C3", C3)):
        y = [r[key] for r in rows]
        a1.plot(x, y, color=col, lw=2, marker="o", ms=4)
        j = int(np.argmin(x)) if invert else int(np.argmax(x))   # right-hand end as drawn
        a1.annotate("Ethane (C2)" if key == "C2" else "Propane (C3)", (x[j], y[j]), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=INK)
    a1.set_ylim(0, 105); a1.set_xlim(*(sorted([min(x), max(x)])))
    a1.margins(x=0); lo, hi = min(x), max(x); span = hi - lo
    a1.set_xlim((hi + 0.04*span, lo - 0.24*span) if invert else (lo - 0.04*span, hi + 0.24*span)); style(a1, xlabel, "Recovery to NGL (%)")
    a1.set_title("Ethane and propane recovery", loc="left", fontsize=11, color=INK, fontweight="bold")
    a1.axvline(base_x, color=MUTED, lw=0.8, ls=(0, (3, 3))); a1.text(base_x, 3, " base case", color=MUTED, fontsize=8)
    y = [r["net_kW"] / 1000 for r in rows]
    a2.plot(x, y, color=INK, lw=2, marker="o", ms=4)
    style(a2, xlabel, "Net compression power (MW)")
    a2.set_title("Recompression power after expander credit", loc="left", fontsize=11, color=INK, fontweight="bold")
    a2.axvline(base_x, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    a2.set_ylim(0, max(y) * 1.15)
    if invert: a2.invert_xaxis()
    a2.margins(x=0.04)
    fig.suptitle(title, x=0.01, ha="left", fontsize=13, color=INK, fontweight="bold")
    for a in (a1, a2): a.set_xticks(x[::2])
    fig.tight_layout(rect=(0, 0, 1, 0.9)); fig.subplots_adjust(top=0.8); fig.savefig(f"results/{fname}.png", facecolor="white"); plt.close(fig)

rp = sweep("P_column_bar", np.arange(14, 36.1, 2))
panel(rp, "Demethanizer pressure (bar)", "Lower column pressure recovers more ethane but costs more compression",
      "sensitivity_column_pressure", 22, invert=True)
rt = sweep("T_cold_C", np.arange(-15, -52, -5))
panel(rt, "Cold separator temperature (°C)", "Colder inlet gas sharply raises ethane recovery for a modest rise in compression",
      "sensitivity_cold_temperature", -35, invert=True)
for r in rp: print("P", r["x"], round(r["C2"],1), round(r["C3"],1), round(r["net_kW"]))
for r in rt: print("T", r["x"], round(r["C2"],1), round(r["C3"],1), round(r["net_kW"]))
