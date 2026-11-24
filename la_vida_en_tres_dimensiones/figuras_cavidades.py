#!/usr/bin/env python3
"""Figuras de los Desafíos I, III, IV y V a partir de los CSV de resultados/.

No recalcula nada: lee lo que dejó cavidades.py, así que corre offline.

    python figuras_cavidades.py
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

import cavidades as cv  # noqa: E402

RES = AQUI / "resultados"
GRIS, DESTACADA, LIGANDO = "#B8C2CC", "#C0392B", "#2E86C1"


def leer_interes():
    with open(RES / "residuos_interes.csv") as f:
        filas = list(csv.DictReader(f))
    return {(r["estructura"], int(r["residuo_pdb"])): r for r in filas}


def barras_volumen(ax, clave, residuo_destacado=None, top=10, titulo=None):
    cavs = cv.leer_cavidades_csv(RES / f"cavidades_{clave}.csv")[:top]
    colores = []
    for c in cavs:
        nums = {int(r[1:]) for r in c["residuos"]}
        if residuo_destacado is not None and residuo_destacado in nums:
            colores.append(DESTACADA)
        elif c["contacto_ligando"]:
            colores.append(LIGANDO)
        else:
            colores.append(GRIS)
    ax.bar([str(c["cavidad"]) for c in cavs], [c["volumen"] for c in cavs], color=colores)
    ax.set_xlabel("Cavidad (ordenadas por volumen)")
    ax.set_ylabel("Volumen (Å³)")
    ax.set_title(titulo or clave, fontsize=10)
    ax.grid(axis="y", alpha=0.3)


def leyenda(ax, destacada):
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=DESTACADA, label=destacada),
                       Patch(color=LIGANDO, label="toca el ligando del cristal"),
                       Patch(color=GRIS, label="otras")], fontsize=8, loc="upper right")


def figura_egfr(interes):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw={"width_ratios": [1.2, 1.2, 0.8]})
    barras_volumen(axes[0], "1M14_A", 813, titulo="1M14 (activo) - Asp837 UniProt = ASP 813 del PDB")
    barras_volumen(axes[1], "3W32_A", 837, titulo="3W32 (inactivo, sin el inhibidor) - ASP 837")
    leyenda(axes[1], "contiene la Asp837 (UniProt)")
    a, b = interes[("1M14_A", 813)], interes[("3W32_A", 837)]
    etiquetas = ["1M14\n(activo)", "3W32\n(inactivo)"]
    pkas = [float(a["pka"]), float(b["pka"])]
    ax = axes[2]
    ax.bar(etiquetas, pkas, color=["#27AE60", "#8E44AD"])
    ax.axhline(3.8, ls="--", color="k", lw=1)
    ax.text(1.45, 3.9, "pKa modelo Asp (3.8)", ha="right", fontsize=8)
    for i, (p, r) in enumerate(zip(pkas, (a, b))):
        ax.text(i, p + 0.1, f"{p:.2f}\n({r['enterrado']}% enterrada)", ha="center", fontsize=8)
    ax.set_ylim(0, 9)
    ax.set_ylabel("pKa (PROPKA)")
    ax.set_title("Asp837 (HRD) en cada confórmero", fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "fig_cavidades_egfr.png", dpi=150)
    plt.close(fig)


def figura_albuminas(interes):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw={"width_ratios": [1.1, 1.1, 1.2]})
    barras_volumen(axes[0], "1AO6_A", 199, titulo="HSA 1AO6_A")
    leyenda(axes[0], "contiene la Lys catalítica")
    barras_volumen(axes[1], "4JK4_A", 221, titulo="BSA 4JK4_A (sin el ligando DIU)")
    pares = [("Lys195", 195, "Arg194", 194), ("Lys199", 199, "Arg198", 198),
             ("Arg222", 222, "Lys221", 221), ("Arg410", 410, "Arg409", 409),
             ("Tyr411", 411, "Tyr410", 410), ("Cys34", 34, "Cys34", 34)]
    ax = axes[2]
    x = range(len(pares))
    sh = [float(interes[("1AO6_A", h)]["shift"]) for _, h, _, _ in pares]
    sb = [float(interes[("4JK4_A", b)]["shift"]) for _, _, _, b in pares]
    ax.bar([i - 0.2 for i in x], sh, width=0.4, label="HSA", color="#E67E22")
    ax.bar([i + 0.2 for i in x], sb, width=0.4, label="BSA", color="#16A085")
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{h}\n{b}" for h, _, b, _ in pares], fontsize=8)
    ax.axhline(0, color="k", lw=0.8)
    for lim in (-1, 1):
        ax.axhline(lim, color="k", lw=0.8, ls=":")
    ax.set_ylabel("pKa predicho - pKa modelo")
    ax.set_title("Corrimiento de pKa (posiciones equivalentes HSA / BSA)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RES / "fig_cavidades_albuminas.png", dpi=150)
    plt.close(fig)


def figura_akr1c4():
    fig, ax = plt.subplots(figsize=(7, 4))
    barras_volumen(ax, "2FVL_B", 55, titulo="AKR1C4 2FVL_B (sin el NADP+)")
    leyenda(ax, "contiene la Tyr55 (sitio activo)")
    fig.tight_layout()
    fig.savefig(RES / "fig_cavidades_akr1c4.png", dpi=150)
    plt.close(fig)


def main():
    interes = leer_interes()
    figura_egfr(interes)
    figura_albuminas(interes)
    figura_akr1c4()
    print("Figuras en", RES)
    return 0


if __name__ == "__main__":
    sys.exit(main())
