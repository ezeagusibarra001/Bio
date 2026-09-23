#!/usr/bin/env python3
"""DESAFÍO III - Buscador de regiones promotoras delimitadas por cajas TATA.

Consigna: tomando como input un archivo con una secuencia de ADN, identificar
las regiones promotoras de un gen, considerando que tal región COMIENZA y
TERMINA con la caja TATA ('TATAAA').

Criterio: se buscan todas las apariciones de la caja y se reporta como región
cada tramo que va desde el inicio de una caja hasta el final de la caja
siguiente (ambas incluidas).
  --modo consecutivas (default): regiones entre cajas vecinas (1-2, 2-3, 3-4...)
  --modo pares: cajas agrupadas de a dos sin compartir (1-2, 3-4...)

El archivo puede ser FASTA (una o varias secuencias) o texto plano.
Las posiciones se informan en base 1, como en las bases de datos biológicas.

Ejemplos:
    python cajas_tata.py datos/ejemplo_adn.fasta
    python cajas_tata.py datos/ejemplo_adn.fasta --modo pares --ambas-hebras
"""

import argparse
import sys

CAJA_TATA = "TATAAA"
COMPLEMENTO = str.maketrans("ACGTN", "TGCAN")


def leer_secuencias(ruta):
    """Lee un archivo FASTA o de texto plano. Devuelve lista de (nombre, secuencia)."""
    secuencias, nombre, partes = [], None, []
    with open(ruta, encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                if partes:
                    secuencias.append((nombre or "secuencia", "".join(partes)))
                nombre, partes = linea[1:].split()[0] if len(linea) > 1 else "secuencia", []
            else:
                partes.append(linea.upper())
    if partes:
        secuencias.append((nombre or "secuencia", "".join(partes)))
    return secuencias


def validar_adn(secuencia):
    invalidos = sorted(set(secuencia) - set("ACGTN"))
    if invalidos:
        raise ValueError(f"La secuencia tiene caracteres que no son ADN: {', '.join(invalidos)}")


def reverso_complementario(secuencia):
    return secuencia.translate(COMPLEMENTO)[::-1]


def posiciones_motivo(secuencia, motivo=CAJA_TATA):
    """Índices (base 0) donde empieza el motivo, incluyendo solapamientos."""
    posiciones, inicio = [], secuencia.find(motivo)
    while inicio != -1:
        posiciones.append(inicio)
        inicio = secuencia.find(motivo, inicio + 1)
    return posiciones


def regiones_promotoras(secuencia, motivo=CAJA_TATA, modo="consecutivas"):
    """Lista de dicts con inicio/fin (base 1, inclusivos) y la secuencia de cada región."""
    cajas = posiciones_motivo(secuencia, motivo)
    paso = 1 if modo == "consecutivas" else 2
    regiones = []
    for i in range(0, len(cajas) - 1, paso):
        inicio, fin = cajas[i], cajas[i + 1] + len(motivo)  # fin exclusivo, base 0
        regiones.append({"inicio": inicio + 1, "fin": fin, "largo": fin - inicio,
                         "secuencia": secuencia[inicio:fin]})
    return regiones


def main(argv=None):
    parser = argparse.ArgumentParser(description="Identifica regiones promotoras delimitadas por cajas TATA.")
    parser.add_argument("archivo", help="Archivo FASTA o de texto con la secuencia de ADN")
    parser.add_argument("--modo", choices=["consecutivas", "pares"], default="consecutivas")
    parser.add_argument("--motivo", default=CAJA_TATA, help="Motivo a buscar (default: TATAAA)")
    parser.add_argument("--ambas-hebras", action="store_true",
                        help="Buscar también en la hebra reversa complementaria")
    parser.add_argument("--min-largo", type=int, default=0, help="Descartar regiones más cortas")
    args = parser.parse_args(argv)

    try:
        secuencias = leer_secuencias(args.archivo)
    except FileNotFoundError:
        print(f"Error: no existe el archivo '{args.archivo}'", file=sys.stderr)
        return 1
    if not secuencias:
        print("Error: el archivo no contiene secuencias.", file=sys.stderr)
        return 1

    motivo = args.motivo.upper()
    total = 0
    for nombre, secuencia in secuencias:
        try:
            validar_adn(secuencia)
        except ValueError as error:
            print(f"[{nombre}] {error}", file=sys.stderr)
            continue

        hebras = [("+", secuencia)]
        if args.ambas_hebras:
            hebras.append(("-", reverso_complementario(secuencia)))

        print(f"== {nombre} ({len(secuencia)} pb) ==")
        for signo, hebra in hebras:
            cajas = posiciones_motivo(hebra, motivo)
            regiones = [r for r in regiones_promotoras(hebra, motivo, args.modo)
                        if r["largo"] >= args.min_largo]
            print(f"  Hebra {signo}: {len(cajas)} caja(s) {motivo} en posiciones "
                  f"{[p + 1 for p in cajas] or '-'}")
            for n, region in enumerate(regiones, 1):
                print(f"    Región {n}: {region['inicio']}-{region['fin']} "
                      f"({region['largo']} pb)\n      {region['secuencia']}")
            if len(cajas) < 2:
                print("    No hay suficientes cajas para delimitar una región.")
            total += len(regiones)
    print(f"\nTotal de regiones encontradas: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
