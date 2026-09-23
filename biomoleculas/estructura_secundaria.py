#!/usr/bin/env python3
"""DESAFÍO V - Predicción de estructura secundaria (método de Chou-Fasman simplificado).

Para cada residuo predice:  H = hélice alfa,  B = hoja beta,  L = loop/bucle.

Idea: cada aminoácido tiene una *propensión* a formar parte de una hélice (Pα)
o de una hoja beta (Pβ), calculada por Chou y Fasman (1974) a partir de
estructuras de proteínas conocidas. Un valor > 1 indica que el aminoácido
aparece en esa estructura más seguido de lo esperado por azar.

Algoritmo:
  1. Para cada residuo se promedian Pα y Pβ en una ventana centrada
     (6 residuos para hélice, 5 para hoja beta, como en el método original).
  2. Si <Pα> ≥ umbral_H y <Pα> ≥ <Pβ>  → H
     si <Pβ> ≥ umbral_B y <Pβ> > <Pα>  → B
     si no                              → L
  3. Se suavizan los segmentos: una hélice de menos de 4 residuos o una hoja
     de menos de 3 se convierten en L (no son estables tan cortas).

Inputs:  secuencia por argumento, archivo FASTA (--fasta) o stdin.
Outputs: texto alineado (default), fasta, csv o json.

Ejemplos:
    python estructura_secundaria.py MLPGLALLLLAAWTMRALEVPTDGNAPLLVEPQIAMFCGR
    python estructura_secundaria.py --fasta datos/secuencias.fasta --formato csv
"""

import argparse
import csv
import json
import sys

# Propensiones de Chou-Fasman (valores × 100): (P_alfa, P_beta)
PROPENSIONES = {
    "A": (142, 83), "R": (98, 93), "N": (67, 89), "D": (101, 54), "C": (70, 119),
    "Q": (111, 110), "E": (151, 37), "G": (57, 75), "H": (100, 87), "I": (108, 160),
    "L": (121, 130), "K": (114, 74), "M": (145, 105), "F": (113, 138), "P": (57, 55),
    "S": (77, 75), "T": (83, 119), "W": (108, 137), "Y": (69, 147), "V": (106, 170),
}


def validar(secuencia):
    secuencia = "".join(secuencia.split()).upper()
    invalidos = sorted(set(secuencia) - set(PROPENSIONES))
    if not secuencia:
        raise ValueError("La secuencia está vacía.")
    if invalidos:
        raise ValueError(f"Aminoácidos no reconocidos: {', '.join(invalidos)}")
    return secuencia


def promedio_en_ventana(valores, i, ancho):
    """Promedio de `valores` en una ventana de `ancho` centrada en i (recortada en bordes)."""
    desde = max(0, i - ancho // 2)
    hasta = min(len(valores), desde + ancho)
    desde = max(0, hasta - ancho)
    ventana = valores[desde:hasta]
    return sum(ventana) / len(ventana)


def suavizar(prediccion, minimo={"H": 4, "B": 3}):
    """Convierte en 'L' los segmentos H/B más cortos que el mínimo."""
    resultado, i = list(prediccion), 0
    while i < len(resultado):
        j = i
        while j < len(resultado) and resultado[j] == resultado[i]:
            j += 1
        if resultado[i] in minimo and j - i < minimo[resultado[i]]:
            resultado[i:j] = "L" * (j - i)
        i = j
    return "".join(resultado)


def predecir(secuencia, umbral_h=1.03, umbral_b=1.05, ventana_h=6, ventana_b=5):
    """Devuelve (prediccion, detalle) donde detalle tiene <Pα> y <Pβ> por residuo."""
    secuencia = validar(secuencia)
    p_alfa = [PROPENSIONES[aa][0] / 100 for aa in secuencia]
    p_beta = [PROPENSIONES[aa][1] / 100 for aa in secuencia]

    estados, detalle = [], []
    for i, aa in enumerate(secuencia):
        alfa = promedio_en_ventana(p_alfa, i, ventana_h)
        beta = promedio_en_ventana(p_beta, i, ventana_b)
        if alfa >= umbral_h and alfa >= beta:
            estado = "H"
        elif beta >= umbral_b and beta > alfa:
            estado = "B"
        else:
            estado = "L"
        estados.append(estado)
        detalle.append({"posicion": i + 1, "residuo": aa,
                        "p_alfa": round(alfa, 3), "p_beta": round(beta, 3)})

    prediccion = suavizar("".join(estados))
    for fila, estado in zip(detalle, prediccion):
        fila["estructura"] = estado
    return prediccion, detalle


def leer_fasta(texto):
    registros, nombre, partes = [], None, []
    for linea in texto.splitlines():
        linea = linea.strip()
        if linea.startswith(">"):
            if partes:
                registros.append((nombre, "".join(partes)))
            nombre, partes = linea[1:].split()[0], []
        elif linea:
            partes.append(linea)
    if partes:
        registros.append((nombre or "secuencia", "".join(partes)))
    return registros


def resumen(prediccion):
    total = len(prediccion)
    return {e: round(100 * prediccion.count(e) / total, 1) for e in "HBL"}


def imprimir_texto(nombre, secuencia, prediccion, ancho=60):
    print(f"> {nombre} ({len(secuencia)} aa)")
    for i in range(0, len(secuencia), ancho):
        print(f"{i + 1:>5}  {secuencia[i:i + ancho]}")
        print(f"{'':>5}  {prediccion[i:i + ancho]}")
    porcentajes = resumen(prediccion)
    print(f"       Hélice: {porcentajes['H']}%  Hoja β: {porcentajes['B']}%  Loop: {porcentajes['L']}%\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Predice estructura secundaria (H/B/L) por Chou-Fasman.")
    parser.add_argument("secuencia", nargs="?", help="Secuencia proteica (código de 1 letra)")
    parser.add_argument("--fasta", help="Archivo FASTA con una o más proteínas")
    parser.add_argument("--formato", choices=["texto", "fasta", "csv", "json"], default="texto")
    parser.add_argument("--umbral-h", type=float, default=1.03)
    parser.add_argument("--umbral-b", type=float, default=1.05)
    args = parser.parse_args(argv)

    if args.fasta:
        with open(args.fasta, encoding="utf-8") as archivo:
            registros = leer_fasta(archivo.read())
    elif args.secuencia:
        registros = [("secuencia", args.secuencia)]
    elif not sys.stdin.isatty():
        texto = sys.stdin.read()
        registros = leer_fasta(texto) if texto.lstrip().startswith(">") else [("stdin", texto)]
    else:
        parser.error("Pasá una secuencia, un --fasta o algo por stdin.")

    resultados = []
    for nombre, secuencia in registros:
        try:
            prediccion, detalle = predecir(secuencia, args.umbral_h, args.umbral_b)
        except ValueError as error:
            print(f"[{nombre}] Error: {error}", file=sys.stderr)
            return 1
        resultados.append((nombre, validar(secuencia), prediccion, detalle))

    if args.formato == "texto":
        for nombre, secuencia, prediccion, _ in resultados:
            imprimir_texto(nombre, secuencia, prediccion)
    elif args.formato == "fasta":
        for nombre, _, prediccion, _ in resultados:
            print(f">{nombre}_estructura_secundaria\n{prediccion}")
    elif args.formato == "csv":
        escritor = csv.writer(sys.stdout)
        escritor.writerow(["proteina", "posicion", "residuo", "p_alfa", "p_beta", "estructura"])
        for nombre, _, _, detalle in resultados:
            for fila in detalle:
                escritor.writerow([nombre, *fila.values()])
    else:
        print(json.dumps([{"proteina": n, "secuencia": s, "prediccion": p, "porcentajes": resumen(p)}
                          for n, s, p, _ in resultados], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
