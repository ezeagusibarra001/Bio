#!/usr/bin/env python3
"""DESAFÍO II - De proteína a ARN.

Toma una secuencia proteica (código de una letra) e imprime una cadena de
ARNm que la codifica. Como el código genético es *degenerado* (varios codones
para un mismo aminoácido) no hay una única respuesta: se puede elegir el
primer codón de la tabla o uno al azar.

Ejemplos:
    python proteina_a_arn.py ATVEKGGKHKTGPNEKGKKIFVQKCSQCHTVLHGLFGRKTGQA
    python proteina_a_arn.py ATVEK --modo aleatorio --semilla 42 --inicio --stop
    python proteina_a_arn.py MW --todas
"""

import argparse
import random
import sys
from itertools import product

from codigo_genetico import (AMINOACIDO_A_CODONES, CODON_INICIO,
                             CODONES_STOP, traducir)

AMINOACIDOS_VALIDOS = set(AMINOACIDO_A_CODONES) - {"*"}


def validar_proteina(secuencia):
    """Normaliza la secuencia y levanta ValueError si tiene letras inválidas."""
    secuencia = secuencia.strip().upper().replace(" ", "")
    invalidos = sorted(set(secuencia) - AMINOACIDOS_VALIDOS)
    if not secuencia:
        raise ValueError("La secuencia está vacía.")
    if invalidos:
        raise ValueError(f"Caracteres no válidos en la secuencia: {', '.join(invalidos)}")
    return secuencia


def proteina_a_arn(secuencia, modo="primero", inicio=False, stop=False, rng=None):
    """Devuelve una cadena de ARNm que codifica `secuencia`.

    modo: 'primero' usa siempre el primer codón de la tabla,
          'aleatorio' elige uno al azar entre los sinónimos.
    inicio: agrega AUG al principio si la proteína no empieza con M.
    stop: agrega un codón de terminación al final.
    """
    secuencia = validar_proteina(secuencia)
    rng = rng or random.Random()
    elegir = (lambda codones: codones[0]) if modo == "primero" else rng.choice

    codones = [elegir(AMINOACIDO_A_CODONES[aa]) for aa in secuencia]
    if inicio and secuencia[0] != "M":
        codones.insert(0, CODON_INICIO)
    if stop:
        codones.append(elegir(CODONES_STOP))
    return "".join(codones)


def cantidad_de_arns_posibles(secuencia):
    """Cuántos ARNm distintos codifican la proteína (sin contar inicio/stop)."""
    total = 1
    for aa in validar_proteina(secuencia):
        total *= len(AMINOACIDO_A_CODONES[aa])
    return total


def todas_las_opciones(secuencia):
    """Generador con todos los ARNm posibles (¡crece exponencialmente!)."""
    opciones = [AMINOACIDO_A_CODONES[aa] for aa in validar_proteina(secuencia)]
    for combinacion in product(*opciones):
        yield "".join(combinacion)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Genera un ARNm que codifica una proteína dada.")
    parser.add_argument("proteina", help="Secuencia proteica en código de una letra")
    parser.add_argument("--modo", choices=["primero", "aleatorio"], default="primero",
                        help="Cómo elegir entre codones sinónimos (default: primero)")
    parser.add_argument("--semilla", type=int, help="Semilla para el modo aleatorio")
    parser.add_argument("--inicio", action="store_true", help="Agregar codón de inicio AUG")
    parser.add_argument("--stop", action="store_true", help="Agregar codón de stop")
    parser.add_argument("--todas", action="store_true",
                        help="Listar TODOS los ARNm posibles (solo para péptidos cortos)")
    args = parser.parse_args(argv)

    try:
        if args.todas:
            n = cantidad_de_arns_posibles(args.proteina)
            if n > 10_000:
                parser.error(f"Hay {n:,} combinaciones, demasiadas para listar.")
            for arn in todas_las_opciones(args.proteina):
                print(arn)
            return 0

        arn = proteina_a_arn(args.proteina, args.modo, args.inicio, args.stop,
                             random.Random(args.semilla))
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(arn)
    n = cantidad_de_arns_posibles(args.proteina)
    print(f"\n# Largo: {len(arn)} nt | Verificación (traducción): {traducir(arn)}", file=sys.stderr)
    print(f"# Existen {n:.3e} ARNm distintos que codifican esta proteína.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
