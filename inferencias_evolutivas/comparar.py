#!/usr/bin/env python3
"""DESAFÍO I - Comparación de bartmosca con humano y mosca (y con el resto).

Trabaja sobre un alineamiento múltiple ya hecho (Clustal Omega o MAFFT, en
FASTA alineado) y calcula:

  1. Identidad y p-distancia de a pares.
       identidad = idénticos / columnas donde ninguna de las dos tiene gap
       p-distancia = 1 - identidad   (proporción de sitios distintos)
  2. Para bartmosca, una clasificación columna por columna: ¿el residuo
     coincide con el humano, con la mosca, con ambos o con ninguno?
  3. Una ventana deslizante de identidad de bartmosca contra humano y contra
     mosca, para ver si el parecido es parejo o si hay regiones (mosaico).

Salidas en resultados/:
  - identidad_<set>.csv, pdist_<set>.csv
  - mosaico_bartmosca.csv  (columna a columna)
  - fig_identidad_<set>.png, fig_mosaico_bartmosca.png

Uso:
    python comparar.py                                  # usa los alineamientos guardados
    python comparar.py --aln resultados/aln_mafft_set3.fasta --set set3
"""

import argparse
import csv
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

from secuencias import leer_fasta  # noqa: E402

RESULTADOS = AQUI / "resultados"
GAP = "-"
HUMANO, MOSCA, BART = "Homo_sapiens", "Drosophila_melanogaster", "bartmosca"
DISTAL = "Drosophila_distal"


# ---------- Métricas de a pares ----------
def identidad(a, b):
    """Proporción de columnas idénticas entre dos secuencias alineadas
    (sin contar columnas donde alguna tiene gap)."""
    if len(a) != len(b):
        raise ValueError("Las secuencias alineadas deben tener el mismo largo")
    columnas = [(x, y) for x, y in zip(a, b) if x != GAP and y != GAP]
    if not columnas:
        return 0.0
    return sum(x == y for x, y in columnas) / len(columnas)


def p_distancia(a, b):
    return 1 - identidad(a, b)


def matriz(alineamiento, funcion):
    nombres = [n for n, _ in alineamiento]
    secs = dict(alineamiento)
    return nombres, [[funcion(secs[i], secs[j]) for j in nombres] for i in nombres]


def diferencias(a, b):
    """Lista de (columna, residuo_a, residuo_b) donde difieren (columnas sin gap)."""
    return [(k + 1, x, y) for k, (x, y) in enumerate(zip(a, b))
            if x != GAP and y != GAP and x != y]


# ---------- Mosaico ----------
def clasificar_columnas(bart, humano, mosca):
    """Para cada columna donde bartmosca tiene residuo, dice a quién se parece."""
    filas = []
    posicion = 0
    for col, (b, h, m) in enumerate(zip(bart, humano, mosca), start=1):
        if b == GAP:
            continue
        posicion += 1
        if b == h and b == m:
            clase = "ambos"
        elif b == h:
            clase = "solo_humano"
        elif b == m:
            clase = "solo_mosca"
        else:
            clase = "ninguno"
        filas.append({"columna": col, "posicion_bart": posicion,
                      "bart": b, "humano": h, "mosca": m, "clase": clase})
    return filas


def ventana_identidad(a, b, ancho=15):
    """Identidad en ventanas deslizantes (centradas) sobre el alineamiento."""
    valores = []
    for inicio in range(0, len(a) - ancho + 1):
        valores.append((inicio + ancho // 2 + 1, identidad(a[inicio:inicio + ancho], b[inicio:inicio + ancho])))
    return valores


def punto_de_corte(quimera, izquierda, derecha):
    """Modelo de quimera de un solo corte: los residuos 1..k vienen de
    `izquierda` y los k+1..fin de `derecha`. Para cada k cuenta cuántos
    residuos de la quimera quedan explicados y devuelve
    (mejor_puntaje, total_de_residuos, [posiciones k con el mejor puntaje]).
    Las posiciones son de la quimera (sin contar gaps)."""
    columnas = [(q, i, d) for q, i, d in zip(quimera, izquierda, derecha) if q != GAP]
    total = len(columnas)
    puntajes = []
    for k in range(total + 1):
        izq = sum(q == i for q, i, _ in columnas[:k])
        der = sum(q == d for q, _, d in columnas[k:])
        puntajes.append(izq + der)
    mejor = max(puntajes)
    return mejor, total, [k for k, p in enumerate(puntajes) if p == mejor]


def resumen_regiones(filas):
    """Busca el último residuo 'solo_mosca': a partir de ahí bartmosca es humano."""
    solo_mosca = [f["posicion_bart"] for f in filas if f["clase"] == "solo_mosca"]
    corte = max(solo_mosca) if solo_mosca else 0
    n_term = [f for f in filas if f["posicion_bart"] <= corte]
    resto = [f for f in filas if f["posicion_bart"] > corte]

    def ident(grupo, clave):
        return sum(f["bart"] == f[clave] for f in grupo) / len(grupo) if grupo else 0.0

    return {
        "corte": corte,
        "n_term_largo": len(n_term),
        "n_term_id_humano": ident(n_term, "humano"),
        "n_term_id_mosca": ident(n_term, "mosca"),
        "resto_largo": len(resto),
        "resto_id_humano": ident(resto, "humano"),
        "resto_id_mosca": ident(resto, "mosca"),
    }


# ---------- Salidas ----------
def guardar_matriz(nombres, valores, ruta, decimales=3):
    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + nombres)
        for n, fila in zip(nombres, valores):
            w.writerow([n] + [round(v, decimales) for v in fila])


def figura_identidad(nombres, valores, ruta, titulo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.0 * len(nombres) + 2.5, 0.9 * len(nombres) + 1.5))
    im = ax.imshow([[100 * v for v in fila] for fila in valores], cmap="viridis", vmin=40, vmax=100)
    cortos = [n.replace("_", " ") for n in nombres]
    ax.set_xticks(range(len(nombres)), cortos, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(nombres)), cortos, fontsize=8)
    for i, fila in enumerate(valores):
        for j, v in enumerate(fila):
            ax.text(j, i, f"{100 * v:.0f}", ha="center", va="center", fontsize=7,
                    color="white" if v < 0.8 else "black")
    fig.colorbar(im, ax=ax, label="% identidad")
    ax.set_title(titulo, fontsize=10)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def figura_mosaico(filas, ventanas_h, ventanas_m, ruta, nombre_mosca="mosca"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colores = {"ambos": "#9e9e9e", "solo_humano": "#1f77b4", "solo_mosca": "#d62728", "ninguno": "#ffbf00"}
    etiquetas = {"ambos": "igual a ambos", "solo_humano": "solo humano",
                 "solo_mosca": "solo mosca", "ninguno": "ninguno"}

    from matplotlib.patches import Patch

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 2]})
    handles = []
    for clase, color in colores.items():
        xs = [f["columna"] for f in filas if f["clase"] == clase]
        ax1.bar(xs, [1] * len(xs), width=1.0, color=color)
        handles.append(Patch(color=color, label=f"{etiquetas[clase]} ({len(xs)})"))
    for f in filas:
        if f["clase"] in ("solo_mosca", "ninguno"):
            ax1.text(f["columna"], 1.08, f["bart"], ha="center", fontsize=7)
    ax1.set_yticks([])
    ax1.set_ylim(0, 1.3)
    ax1.legend(handles=handles, ncol=4, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle(f"bartmosca, columna por columna: ¿humano o mosca ({nombre_mosca})?", fontsize=11)

    ax2.plot(*zip(*ventanas_h), color=colores["solo_humano"], lw=2, label="bartmosca vs humano")
    ax2.plot(*zip(*ventanas_m), color=colores["solo_mosca"], lw=2, label=f"bartmosca vs {nombre_mosca}")
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("identidad (ventana 15)")
    ax2.set_xlabel("columna del alineamiento")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def analizar(ruta_aln, nombre_set, figuras=True):
    alineamiento = leer_fasta(ruta_aln)
    RESULTADOS.mkdir(exist_ok=True)

    nombres, ident = matriz(alineamiento, identidad)
    _, pdist = matriz(alineamiento, p_distancia)
    guardar_matriz(nombres, ident, RESULTADOS / f"identidad_{nombre_set}.csv")
    guardar_matriz(nombres, pdist, RESULTADOS / f"pdist_{nombre_set}.csv")
    if figuras:
        figura_identidad(nombres, ident, RESULTADOS / f"fig_identidad_{nombre_set}.png",
                         f"Identidad de a pares ({nombre_set}, {Path(ruta_aln).name})")

    print(f"\n== {nombre_set} ({Path(ruta_aln).name}) ==")
    i_b = nombres.index(BART)
    for j, n in sorted(enumerate(nombres), key=lambda t: -ident[i_b][t[0]]):
        if n != BART:
            print(f"  bartmosca vs {n:<26} identidad {100 * ident[i_b][j]:5.1f}%   p-dist {pdist[i_b][j]:.3f}")

    secs = dict(alineamiento)
    # set3: mosaico contra la mosca del archivo del curso (isoforma proximal).
    # set4: mosaico contra la isoforma distal (el parálogo).
    mosca = {"set3": MOSCA, "set4": DISTAL}.get(nombre_set)
    if mosca and {HUMANO, mosca, BART} <= set(secs):
        sufijo = "" if nombre_set == "set3" else "_distal"
        filas = clasificar_columnas(secs[BART], secs[HUMANO], secs[mosca])
        with open(RESULTADOS / f"mosaico_bartmosca{sufijo}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0]))
            w.writeheader()
            w.writerows(filas)
        reg = resumen_regiones(filas)
        conteo = {c: sum(f["clase"] == c for f in filas) for c in ("ambos", "solo_humano", "solo_mosca", "ninguno")}
        print(f"  mosaico vs {mosca}: {conteo}")
        print(f"  N-terminal (pos 1-{reg['corte']}, {reg['n_term_largo']} aa): "
              f"id humano {100 * reg['n_term_id_humano']:.1f}%  id mosca {100 * reg['n_term_id_mosca']:.1f}%")
        print(f"  resto (pos {reg['corte'] + 1}-, {reg['resto_largo']} aa): "
              f"id humano {100 * reg['resto_id_humano']:.1f}%  id mosca {100 * reg['resto_id_mosca']:.1f}%")
        mejor, total, cortes = punto_de_corte(secs[BART], secs[mosca], secs[HUMANO])
        print(f"  quimera de un corte (mosca -> humano): explica {mejor}/{total} residuos, "
              f"corte después de la posición {cortes[0]}..{cortes[-1]}")
        pos_bart = {f["columna"]: f["posicion_bart"] for f in filas}
        print("  difiere del humano en (posición en bartmosca, bart, humano):",
              [(pos_bart[c], b, h) for c, b, h in diferencias(secs[BART], secs[HUMANO])])
        if figuras:
            figura_mosaico(filas, ventana_identidad(secs[BART], secs[HUMANO]),
                           ventana_identidad(secs[BART], secs[mosca]),
                           RESULTADOS / f"fig_mosaico_bartmosca{sufijo}.png",
                           mosca.replace("_", " "))
    return nombres, ident


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--aln", help="Alineamiento en FASTA (si no, usa los guardados)")
    parser.add_argument("--set", default="set3", help="Nombre del set (set3 / set9)")
    args = parser.parse_args(argv)

    if args.aln:
        analizar(args.aln, args.set)
        return 0
    for nombre in ("set3", "set4", "set9"):
        candidatos = [RESULTADOS / f"aln_{nombre}.fasta", RESULTADOS / f"aln_mafft_{nombre}.fasta"]
        ruta = next((c for c in candidatos if c.exists()), None)
        if ruta is None:
            print(f"No hay alineamiento para {nombre}: corré clustal_ebi.py o mafft_local.py")
            continue
        analizar(ruta, nombre)
    return 0


if __name__ == "__main__":
    sys.exit(main())
