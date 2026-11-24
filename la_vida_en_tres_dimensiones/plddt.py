#!/usr/bin/env python3
"""DESAFÍO VI - Graficar el pLDDT por posición de un modelo de AlphaFold2.

Acepta tres formatos de entrada (se detecta solo por la extensión / contenido):

  1. JSON de ColabFold:  <jobname>_scores_rank_001_alphafold2_ptm_model_X_seed_000.json
     -> clave "plddt" (una lista con un valor por residuo; también trae "pae", "ptm").
  2. JSON de confianza de AlphaFold DB: AF-...-confidence_v1.json
     -> clave "confidenceScore".
  3. Un PDB de AlphaFold/ColabFold: el pLDDT viene guardado en la columna del
     factor B (se toma el del carbono alfa de cada residuo).

Además, si se pasa --pae con el JSON de PAE (formato AlphaFold DB o el mismo
scores.json de ColabFold), grafica el mapa de error alineado por pares.

Ejemplos:
    python plddt.py                                    # modelo de E7 guardado en datos/alphafold/
    python plddt.py MiCorrida_scores_rank_001_alphafold2_ptm_model_1_seed_000.json
    python plddt.py MiCorrida_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb --salida e7.png
"""

import argparse
import json
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

import explorar_pdb as ep  # noqa: E402

DIR_AF = AQUI / "datos" / "alphafold"
MODELO_E7 = DIR_AF / "AF-0000000365760270-model_v1.pdb"
PAE_E7 = DIR_AF / "AF-0000000365760270-predicted_aligned_error_v1.json"

# Bandas y colores oficiales de AlphaFold
CATEGORIAS = [  # (mínimo, etiqueta, color)
    (90, "muy alta (>90)", "#0053D6"),
    (70, "confiable (70-90)", "#65CBF3"),
    (50, "baja (50-70)", "#FFDB13"),
    (0, "muy baja (<50)", "#FF7D45"),
]


def leer_plddt(ruta):
    """Devuelve la lista de pLDDT (un valor por residuo) desde JSON de ColabFold/AFDB o PDB."""
    ruta = Path(ruta)
    if ruta.suffix.lower() == ".json":
        with open(ruta) as f:
            d = json.load(f)
        if isinstance(d, list):
            d = d[0]
        for clave in ("plddt", "confidenceScore"):
            if clave in d:
                return [float(x) for x in d[clave]]
        raise ValueError(f"{ruta.name}: no tiene clave 'plddt' ni 'confidenceScore'")
    # PDB: factor B del CA
    plddt = []
    for a in ep.leer_atomos(ruta, solo_proteina=True):
        if a["nombre"] == "CA":
            plddt.append(a["bfactor"])
    if not plddt:
        raise ValueError(f"{ruta.name}: no encontré carbonos alfa")
    return plddt


def leer_pae(ruta):
    """Matriz de PAE (lista de listas) desde el JSON de AFDB o el scores.json de ColabFold."""
    with open(ruta) as f:
        d = json.load(f)
    if isinstance(d, list):
        d = d[0]
    for clave in ("predicted_aligned_error", "pae"):
        if clave in d:
            return d[clave]
    raise ValueError(f"{Path(ruta).name}: no tiene PAE")


def categoria(valor):
    for minimo, etiqueta, _ in CATEGORIAS:
        if valor > minimo or minimo == 0:
            return etiqueta
    return CATEGORIAS[-1][1]


def resumen(plddt):
    """Media y fracción de residuos en cada categoría de confianza."""
    n = len(plddt)
    frac = {etq: 0 for _, etq, _ in CATEGORIAS}
    for v in plddt:
        frac[categoria(v)] += 1
    return {"n": n, "media": round(sum(plddt) / n, 2),
            "fracciones": {k: round(v / n, 2) for k, v in frac.items()}}


def segmentos(plddt, umbral=70.0, largo_min=3):
    """Tramos contiguos (1-indexados) con pLDDT >= umbral ('alta') o < umbral ('baja').

    Los tramos más cortos que largo_min se absorben en el tramo anterior para no
    fragmentar por un residuo suelto.
    """
    crudos = []
    for i, v in enumerate(plddt, start=1):
        tipo = "alta" if v >= umbral else "baja"
        if crudos and crudos[-1][0] == tipo:
            crudos[-1][2] = i
        else:
            crudos.append([tipo, i, i])
    limpios = []
    for t, ini, fin in crudos:
        if limpios and (fin - ini + 1) < largo_min:
            limpios[-1][2] = fin
        elif limpios and limpios[-1][0] == t:
            limpios[-1][2] = fin
        else:
            limpios.append([t, ini, fin])
    return [tuple(s) for s in limpios]


def graficar(plddt, salida, titulo="pLDDT por posición", pae=None, secuencia=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos = list(range(1, len(plddt) + 1))
    ncols = 2 if pae is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(13 if pae is not None else 10, 4.6),
                             gridspec_kw={"width_ratios": [1.7, 1]} if pae is not None else None)
    ax = axes[0] if pae is not None else axes
    limites = [100, 90, 70, 50, 0]
    for (minimo, etq, color), techo in zip(CATEGORIAS, limites):
        ax.axhspan(minimo, techo, color=color, alpha=0.15, lw=0)
    colores = [next(c for m, _, c in CATEGORIAS if v > m or m == 0) for v in plddt]
    ax.plot(pos, plddt, color="#222222", lw=1.3, zorder=2)
    ax.scatter(pos, plddt, c=colores, s=14, zorder=3, edgecolors="none")
    for minimo, etq, color in CATEGORIAS:
        ax.scatter([], [], c=color, s=30, label=etq)
    ax.set_xlim(1, len(plddt))
    ax.set_ylim(0, 100)
    ax.set_xlabel("Posición del residuo")
    ax.set_ylabel("pLDDT")
    ax.set_title(titulo)
    ax.grid(alpha=0.3, ls=":")
    ax.legend(loc="lower right", fontsize=8, title="Confianza", title_fontsize=8)
    if secuencia and len(secuencia) == len(plddt):
        sec = ax.secondary_xaxis("top")
        sec.set_xticks(pos[::5])
        sec.set_xticklabels(secuencia[::5], fontsize=7, family="monospace")

    if pae is not None:
        axp = axes[1]
        im = axp.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, origin="upper",
                        extent=(0.5, len(pae) + 0.5, len(pae) + 0.5, 0.5))
        axp.set_title("PAE (Å)")
        axp.set_xlabel("Residuo alineado")
        axp.set_ylabel("Residuo evaluado")
        fig.colorbar(im, ax=axp, fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(salida, dpi=150)
    plt.close(fig)
    return salida


def main(argv=None):
    p = argparse.ArgumentParser(description="Grafica el pLDDT por posición (AlphaFold2/ColabFold).")
    p.add_argument("entrada", nargs="?", default=str(MODELO_E7),
                   help="JSON de ColabFold (scores_rank_001), JSON de confianza de AFDB o PDB")
    p.add_argument("--pae", default=None, help="JSON con la matriz de PAE (opcional)")
    p.add_argument("--salida", default=str(AQUI / "resultados" / "fig_plddt_E7.png"))
    p.add_argument("--titulo", default="Proteína E7 de HPV16 (P03129) - pLDDT por posición")
    a = p.parse_args(argv)

    plddt = leer_plddt(a.entrada)
    pae = None
    if a.pae:
        pae = leer_pae(a.pae)
    elif a.entrada == str(MODELO_E7):
        pae = leer_pae(PAE_E7)
    secuencia = None
    if Path(a.entrada).suffix.lower() == ".pdb":
        secuencia = ep.secuencia_ca(a.entrada)

    r = resumen(plddt)
    print(f"{Path(a.entrada).name}: {r['n']} residuos, pLDDT medio {r['media']}")
    for k, v in r["fracciones"].items():
        print(f"  {k:<18} {v:.0%}")
    for t, ini, fin in segmentos(plddt):
        vals = plddt[ini - 1:fin]
        print(f"  tramo {ini:>3}-{fin:<3} confianza {t:<4}  (media {sum(vals) / len(vals):.1f})")
    Path(a.salida).parent.mkdir(parents=True, exist_ok=True)
    graficar(plddt, a.salida, a.titulo, pae, secuencia)
    print(f"Figura: {a.salida}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
