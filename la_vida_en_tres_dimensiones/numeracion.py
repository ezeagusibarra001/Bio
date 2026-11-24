#!/usr/bin/env python3
"""Numeración de residuos: PDB vs UniProt y equivalencias entre proteínas homólogas.

Es la trampa más común del TP: el mismo residuo puede tener números distintos
según la fuente.

- UniProt numera la proteína *precursora* completa, desde la Met inicial.
- Muchos PDB numeran la proteína *madura* (sin péptido señal / propéptido).
  Ejemplos de este TP (sale del registro DBREF de cada archivo):
    1M14 (EGFR) : UniProt = PDB + 24   -> la Asp837 de UniProt es la ASP 813 del PDB
    3W32 (EGFR) : UniProt = PDB        -> la Asp837 de UniProt es la ASP 837 del PDB
    1AO6 (HSA)  : UniProt = PDB + 24   -> la "Lys199" de la bibliografía es la LYS 199 del PDB
    4JK4 (BSA)  : UniProt = PDB + 24
- Entre especies (albúmina humana vs bovina) la numeración tampoco coincide:
  hace falta un alineamiento para saber qué residuo ocupa la posición equivalente.

Uso:
    python numeracion.py            # equivalencias HSA (1AO6) -> BSA (4JK4)
"""

import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

import explorar_pdb as ep  # noqa: E402

GAP = "-"


def alinear_global(s1, s2, match=2, mismatch=-1, gap=-3):
    """Needleman-Wunsch con gap lineal. Devuelve (a1, a2).

    Alcanza para pares muy parecidos (ortólogos, misma proteína con variantes);
    para relaciones lejanas habría que usar BLOSUM y gaps afines.
    """
    n, m = len(s1), len(s2)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    for i in range(1, n + 1):
        fila, ant = F[i], F[i - 1]
        c1 = s1[i - 1]
        for j in range(1, m + 1):
            d = ant[j - 1] + (match if c1 == s2[j - 1] else mismatch)
            u = ant[j] + gap
            lft = fila[j - 1] + gap
            fila[j] = d if d >= u and d >= lft else (u if u >= lft else lft)
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1.append(s1[i - 1]); a2.append(GAP); i -= 1
        else:
            a1.append(GAP); a2.append(s2[j - 1]); j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def identidad(a1, a2):
    cols = [(x, y) for x, y in zip(a1, a2) if x != GAP and y != GAP]
    return round(100 * sum(x == y for x, y in cols) / len(cols), 1) if cols else 0.0


def mapa_equivalencias(nums1, seq1, nums2, seq2):
    """Alinea dos cadenas y devuelve {num_residuo_1: (num_residuo_2, aa1, aa2)}.

    nums1/nums2 son los números de residuo (del PDB) en el mismo orden que seq1/seq2.
    """
    a1, a2 = alinear_global(seq1, seq2)
    mapa, i, j = {}, 0, 0
    for x, y in zip(a1, a2):
        if x != GAP and y != GAP:
            mapa[nums1[i]] = (nums2[j], x, y)
        if x != GAP:
            i += 1
        if y != GAP:
            j += 1
    return mapa, (a1, a2)


def cadena_numerada(ruta_pdb, cadena):
    """(lista de números, secuencia 1 letra) de los residuos con CA de una cadena."""
    res = ep.residuos_ca(ruta_pdb, cadena)
    return list(res.keys()), ep.tres_a_una(res.values())


def equivalencias_albuminas(residuos_hsa=(34, 195, 199, 222, 410, 411)):
    hsa = AQUI / "datos" / "estructuras" / "1AO6.pdb"
    bsa = AQUI / "datos" / "estructuras" / "4JK4.pdb"
    n1, s1 = cadena_numerada(hsa, "A")
    n2, s2 = cadena_numerada(bsa, "A")
    mapa, aln = mapa_equivalencias(n1, s1, n2, s2)
    filas = [(r, s1[n1.index(r)], *mapa.get(r, (None, None, None))[::2]) for r in residuos_hsa]
    return filas, identidad(*aln), mapa


def main():
    filas, ident, mapa = equivalencias_albuminas()
    print(f"HSA 1AO6_A vs BSA 4JK4_A: identidad {ident}% (columnas sin gap)")
    print("HSA (PDB)  ->  BSA (PDB)")
    for r_h, aa_h, r_b, aa_b in filas:
        print(f"  {aa_h}{r_h:<5} ->  {aa_b}{r_b}")
    inverso = {v[0]: k for k, v in mapa.items()}
    print(f"BSA K221 se alinea con HSA {mapa[inverso[221]][1]}{inverso[221]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
