#!/usr/bin/env python3
"""DESAFÍOS V y VI - Alineamiento global por programación dinámica (Needleman-Wunsch).

Lleva "a lo concreto" los pasos que se piden esquematizar en el Desafío V y
permite reproducir lo que muestra la herramienta de Freiburg del Desafío VI:
cambiando el gap penalty cambia la matriz, el score y el alineamiento óptimo.

Idea del algoritmo:
  1. Se arma una matriz F de (len(sec1)+1) x (len(sec2)+1).
  2. La primera fila y la primera columna se inicializan con penalidades de gap
     acumuladas (alinear contra "nada" cuesta un gap por residuo).
  3. Cada celda F[i][j] se calcula como el MÁXIMO de tres caminos posibles:
       - diagonal : F[i-1][j-1] + (match si coinciden, mismatch si no)
       - arriba   : F[i-1][j]   + gap   (gap en la secuencia horizontal)
       - izquierda: F[i][j-1]   + gap   (gap en la secuencia vertical)
     Por eso es programación dinámica: cada celda reutiliza resultados ya
     calculados de subproblemas más chicos, en lugar de recalcularlos.
  4. El score óptimo global queda en la esquina inferior derecha.
  5. El alineamiento se reconstruye con un *traceback* desde esa esquina hasta
     el origen, siguiendo el camino que generó cada valor.

Por defecto usa match=+1, mismatch=-1, gap=-2 (esquema simple, típico de
nucleótidos). Para proteínas lo correcto es usar una matriz de sustitución
(PAM/BLOSUM); acá se deja el esquema simple para ver la mecánica.

Ejemplos:
    python needleman_wunsch.py BANANA MANZANA
    python needleman_wunsch.py ANA ANANA --gap -1 --matriz
    python needleman_wunsch.py GATTACA GCATGCU --match 1 --mismatch -1 --gap -2
"""

import argparse
import sys

GAP = "-"


def alinear(sec1, sec2, match=1, mismatch=-1, gap=-2):
    """Alinea sec1 (filas) y sec2 (columnas) y devuelve un dict con resultados.

    Claves: 'matriz' (la tabla F), 'score' (score óptimo global),
    'alineamiento' (tupla (a1, a2)) e 'identidad' (% de columnas idénticas).
    """
    sec1, sec2 = sec1.upper(), sec2.upper()
    n, m = len(sec1), len(sec2)

    # Paso 1-2: matriz con la primera fila y columna inicializadas con gaps.
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap

    # Paso 3: se completa cada celda con el máximo de los tres caminos.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            puntaje = match if sec1[i - 1] == sec2[j - 1] else mismatch
            diagonal = F[i - 1][j - 1] + puntaje
            arriba = F[i - 1][j] + gap
            izquierda = F[i][j - 1] + gap
            F[i][j] = max(diagonal, arriba, izquierda)

    # Paso 5: traceback desde la esquina inferior derecha hasta el origen.
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        puntaje = match if (i > 0 and j > 0 and sec1[i - 1] == sec2[j - 1]) else mismatch
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + puntaje:
            a1.append(sec1[i - 1])
            a2.append(sec2[j - 1])
            i, j = i - 1, j - 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1.append(sec1[i - 1])
            a2.append(GAP)
            i -= 1
        else:
            a1.append(GAP)
            a2.append(sec2[j - 1])
            j -= 1

    a1 = "".join(reversed(a1))
    a2 = "".join(reversed(a2))
    return {
        "matriz": F,
        "score": F[n][m],
        "alineamiento": (a1, a2),
        "identidad": identidad(a1, a2),
    }


def identidad(a1, a2):
    """Porcentaje de columnas idénticas (ignora las columnas con gap)."""
    columnas = [(x, y) for x, y in zip(a1, a2) if x != GAP and y != GAP]
    if not columnas:
        return 0.0
    iguales = sum(1 for x, y in columnas if x == y)
    return round(100 * iguales / len(columnas), 1)


def formato_alineamiento(a1, a2):
    """Devuelve las dos secuencias alineadas con una línea de marcas en el medio."""
    marcas = "".join("|" if x == y and x != GAP else " " for x, y in zip(a1, a2))
    return f"{a1}\n{marcas}\n{a2}"


def formato_matriz(F, sec1, sec2):
    """Dibuja la matriz F con las secuencias como encabezados (para el Desafío VI)."""
    sec1, sec2 = sec1.upper(), sec2.upper()
    ancho = max(4, max(len(str(v)) for fila in F for v in fila) + 1)
    encabezado = ["", ""] + list(sec2)
    lineas = ["".join(c.rjust(ancho) for c in encabezado)]
    etiquetas_filas = [""] + list(sec1)
    for etiqueta, fila in zip(etiquetas_filas, F):
        celdas = [etiqueta] + [str(v) for v in fila]
        lineas.append("".join(c.rjust(ancho) for c in celdas))
    return "\n".join(lineas)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Alineamiento global de dos secuencias (Needleman-Wunsch).")
    parser.add_argument("sec1", nargs="?", default="BANANA", help="Primera secuencia")
    parser.add_argument("sec2", nargs="?", default="MANZANA", help="Segunda secuencia")
    parser.add_argument("--match", type=int, default=1, help="Puntaje por coincidencia")
    parser.add_argument("--mismatch", type=int, default=-1, help="Puntaje por discordancia")
    parser.add_argument("--gap", type=int, default=-2, help="Penalidad por gap")
    parser.add_argument("--matriz", action="store_true", help="Imprime la matriz de programación dinámica")
    args = parser.parse_args(argv)

    resultado = alinear(args.sec1, args.sec2, args.match, args.mismatch, args.gap)
    a1, a2 = resultado["alineamiento"]

    if args.matriz:
        print(formato_matriz(resultado["matriz"], args.sec1, args.sec2))
        print()
    print(f"match={args.match}  mismatch={args.mismatch}  gap={args.gap}")
    print(formato_alineamiento(a1, a2))
    print(f"score = {resultado['score']}   identidad = {resultado['identidad']}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
