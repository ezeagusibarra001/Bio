#!/usr/bin/env python3
"""DESAFÍO VII - Búsqueda de motivos proteicos en secuencias de UniProt
(basado en Rosalind "MPRT": https://rosalind.info/problems/mprt/).

Dado: hasta 15 identificadores de UniProt (uno por línea en un archivo o
      como argumentos).
Retorna: para cada proteína que tenga el motivo de N-glicosilación
         N{P}[ST]{P}, su identificador y las posiciones (base 1) donde aparece.

Notación de motivos:
    [XY] -> "X o Y"       {X} -> "cualquier aminoácido excepto X"
Se traduce a una expresión regular. Para encontrar apariciones SOLAPADAS
(importante en este problema) se usa un lookahead: (?=(...)).

Nota: la URL vieja www.uniprot.org/uniprot/<id>.fasta que menciona el TP hoy
redirige a la API REST https://rest.uniprot.org/uniprotkb/<id>.fasta, que es
la que se usa acá. Los IDs estilo Rosalind "P07204_TRBM_HUMAN" se recortan
al accession ("P07204") para la consulta, pero se imprimen completos.

Ejemplos:
    python motivos_uniprot.py --archivo datos/ids_ejemplo.txt
    python motivos_uniprot.py B5ZC00 P07204_TRBM_HUMAN
    python motivos_uniprot.py B5ZC00 --motivo "[ST]x[RK]"      # otro motivo
    python motivos_uniprot.py --fasta-local datos/mis_prot.fasta # sin internet
"""

import argparse
import re
import sys
import time
import urllib.error
import urllib.request

MOTIVO_N_GLICOSILACION = "N{P}[ST]{P}"
URL_UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"
MAX_IDS = 15


def motivo_a_regex(motivo):
    """Convierte la notación de motivos a una regex. 'x' (minúscula) = cualquiera.

    >>> motivo_a_regex("N{P}[ST]{P}")
    'N[^P][ST][^P]'
    """
    patron, i = [], 0
    while i < len(motivo):
        caracter = motivo[i]
        if caracter in "[{":
            cierre = "]" if caracter == "[" else "}"
            fin = motivo.index(cierre, i)
            letras = motivo[i + 1:fin].upper()
            if not letras.isalpha():
                raise ValueError(f"Grupo inválido en el motivo: {motivo[i:fin + 1]}")
            patron.append(f"[{letras}]" if caracter == "[" else f"[^{letras}]")
            i = fin + 1
        elif caracter == "x":
            patron.append(".")
            i += 1
        elif caracter.isalpha():
            patron.append(caracter.upper())
            i += 1
        else:
            raise ValueError(f"Carácter inesperado en el motivo: {caracter!r}")
    return "".join(patron)


def buscar_motivo(secuencia, motivo=MOTIVO_N_GLICOSILACION):
    """Posiciones (base 1) de todas las apariciones, incluyendo solapadas."""
    regex = re.compile(f"(?=({motivo_a_regex(motivo)}))")
    return [m.start() + 1 for m in regex.finditer(secuencia.upper())]


def parsear_fasta(texto):
    """Devuelve dict {id: secuencia}. El id es el accession si el header es de UniProt."""
    secuencias, actual = {}, None
    for linea in texto.splitlines():
        linea = linea.strip()
        if linea.startswith(">"):
            header = linea[1:].split()[0]
            partes = header.split("|")  # formato UniProt: sp|P07204|TRBM_HUMAN
            actual = partes[1] if len(partes) >= 3 else header
            secuencias[actual] = ""
        elif actual:
            secuencias[actual] += linea
    return secuencias


def accession(identificador):
    return identificador.split("_")[0]


def descargar_secuencia(identificador, reintentos=3):
    url = URL_UNIPROT.format(accession(identificador))
    for intento in range(reintentos):
        try:
            with urllib.request.urlopen(url, timeout=20) as respuesta:
                texto = respuesta.read().decode("utf-8")
            secuencias = parsear_fasta(texto)
            if not secuencias:
                raise ValueError(f"UniProt no devolvió secuencia para {identificador}")
            return next(iter(secuencias.values()))
        except urllib.error.HTTPError as error:
            if error.code in (400, 404):
                raise ValueError(f"{identificador} no existe en UniProt") from error
            ultimo = error
        except urllib.error.URLError as error:
            ultimo = error
        time.sleep(1 + intento)
    raise ConnectionError(f"No se pudo descargar {identificador}: {ultimo}")


def leer_ids(ruta):
    with open(ruta, encoding="utf-8") as archivo:
        return [linea.strip() for linea in archivo if linea.strip()]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Busca motivos proteicos en proteínas de UniProt.")
    parser.add_argument("ids", nargs="*", help="Identificadores de UniProt")
    parser.add_argument("--archivo", help="Archivo con un identificador por línea")
    parser.add_argument("--fasta-local", help="Usar un FASTA local en vez de descargar de UniProt")
    parser.add_argument("--motivo", default=MOTIVO_N_GLICOSILACION,
                        help="Motivo en notación Rosalind (default: N{P}[ST]{P})")
    args = parser.parse_args(argv)

    try:
        regex = motivo_a_regex(args.motivo)
    except ValueError as error:
        parser.error(str(error))

    if args.fasta_local:
        with open(args.fasta_local, encoding="utf-8") as archivo:
            secuencias = parsear_fasta(archivo.read())
        ids = args.ids or list(secuencias)
    else:
        ids = args.ids + (leer_ids(args.archivo) if args.archivo else [])
        if not ids:
            parser.error("Indicá identificadores o un --archivo.")
        if len(ids) > MAX_IDS:
            print(f"Aviso: el problema admite hasta {MAX_IDS} IDs; se procesarán todos igual.",
                  file=sys.stderr)
        secuencias = {}

    print(f"# Motivo {args.motivo} -> regex {regex}", file=sys.stderr)
    codigo_salida = 0
    for identificador in ids:
        try:
            secuencia = (secuencias.get(identificador) or secuencias.get(accession(identificador))
                         if args.fasta_local else descargar_secuencia(identificador))
            if not secuencia:
                raise ValueError(f"{identificador} no está en el FASTA local")
        except (ValueError, ConnectionError) as error:
            print(f"Error: {error}", file=sys.stderr)
            codigo_salida = 1
            continue
        posiciones = buscar_motivo(secuencia, args.motivo)
        if posiciones:
            print(identificador)
            print(" ".join(map(str, posiciones)))
    return codigo_salida


if __name__ == "__main__":
    sys.exit(main())
