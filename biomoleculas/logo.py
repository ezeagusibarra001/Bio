#!/usr/bin/env python3
"""DESAFÍO VIII - Sequence logo de un conjunto de secuencias alineadas (logomaker).

Un logo muestra, para cada posición del alineamiento, qué aminoácidos aparecen
y con qué frecuencia. En modo 'informacion' (bits) la altura total de cada
columna indica cuán conservada está la posición: una posición idéntica en
todas las secuencias es alta, una muy variable es baja.

Las secuencias tienen que estar alineadas (mismo largo). Como son largas, el
logo se parte en varias filas.

Ejemplos:
    python logo.py datos/secuencias.fasta
    python logo.py datos/secuencias.fasta --tipo informacion --por-fila 50 -o logo_bits.png
    python logo.py datos/secuencias.fasta --solo-variables
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")  # permite generar la imagen sin pantalla
import matplotlib.pyplot as plt  # noqa: E402
import logomaker  # noqa: E402

from estructura_secundaria import leer_fasta  # noqa: E402

TIPOS = {"frecuencia": "probability", "informacion": "information", "conteo": "counts"}
ETIQUETAS_Y = {"frecuencia": "Frecuencia", "informacion": "Bits", "conteo": "Conteo"}


def validar_alineamiento(secuencias):
    largos = {len(s) for s in secuencias}
    if len(largos) != 1:
        raise ValueError(f"Las secuencias deben tener el mismo largo (hay largos {sorted(largos)}). "
                         "Alinealas primero (por ejemplo con Clustal Omega o MAFFT).")


def posiciones_variables(secuencias):
    """Posiciones (base 1) donde no todas las secuencias tienen el mismo residuo."""
    return [i + 1 for i, columna in enumerate(zip(*secuencias)) if len(set(columna)) > 1]


def matriz_logo(secuencias, tipo="frecuencia", pseudoconteo=0):
    """Matriz posición x aminoácido. Ojo: logomaker usa pseudocount=1 por defecto,
    lo que con pocas secuencias inventa aminoácidos que no están; por eso va en 0."""
    return logomaker.alignment_to_matrix(secuencias, to_type=TIPOS[tipo],
                                         pseudocount=pseudoconteo,
                                         characters_to_ignore="-.X")


def dibujar(secuencias, tipo="frecuencia", por_fila=60, salida="logo.png", titulo=None):
    matriz = matriz_logo(secuencias, tipo)
    matriz.index = matriz.index + 1  # posiciones en base 1
    variables = set(posiciones_variables(secuencias))
    n_filas = -(-len(matriz) // por_fila)

    figura, ejes = plt.subplots(n_filas, 1, figsize=(max(10, por_fila * 0.28), 2.2 * n_filas),
                                squeeze=False)
    for fila, eje in enumerate(ejes[:, 0]):
        tramo = matriz.iloc[fila * por_fila:(fila + 1) * por_fila]
        logo = logomaker.Logo(tramo, ax=eje, color_scheme="chemistry",
                              font_name="DejaVu Sans Mono", vpad=0.05)
        for posicion in variables & set(tramo.index):
            logo.highlight_position(p=posicion, color="gold", alpha=0.35)
        logo.style_spines(visible=False)
        logo.style_spines(spines=["left", "bottom"], visible=True)
        eje.set_ylabel(ETIQUETAS_Y[tipo])
        eje.set_xticks([p for p in tramo.index if p % 10 == 0 or p == tramo.index[0]])
    ejes[-1, 0].set_xlabel("Posición en el alineamiento")
    figura.suptitle(titulo or f"Logo de {len(secuencias)} secuencias "
                    f"(en amarillo: posiciones variables)", fontsize=12)
    figura.tight_layout()
    figura.savefig(salida, dpi=150)
    plt.close(figura)
    return salida


def main(argv=None):
    parser = argparse.ArgumentParser(description="Genera un sequence logo con logomaker.")
    parser.add_argument("fasta", help="FASTA con las secuencias alineadas")
    parser.add_argument("--tipo", choices=TIPOS, default="frecuencia")
    parser.add_argument("--por-fila", type=int, default=60, help="Posiciones por fila del gráfico")
    parser.add_argument("-o", "--salida", default="logo.png")
    parser.add_argument("--solo-variables", action="store_true",
                        help="Solo imprimir las posiciones variables, sin graficar")
    args = parser.parse_args(argv)

    with open(args.fasta, encoding="utf-8") as archivo:
        registros = leer_fasta(archivo.read())
    secuencias = [s.upper() for _, s in registros]
    try:
        validar_alineamiento(secuencias)
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    variables = posiciones_variables(secuencias)
    print(f"{len(secuencias)} secuencias de {len(secuencias[0])} residuos.")
    print(f"Posiciones variables ({len(variables)}):")
    for p in variables:
        print(f"  {p:>4}: " + " / ".join(f"{n}={s[p - 1]}" for n, s in
                                         zip((n for n, _ in registros), secuencias)))
    if not args.solo_variables:
        print(f"Logo guardado en {dibujar(secuencias, args.tipo, args.por_fila, args.salida)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
