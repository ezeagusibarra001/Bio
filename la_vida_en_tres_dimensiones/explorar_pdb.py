#!/usr/bin/env python3
"""PARA PENSAR (ubiquitina, 1UBQ) - Leer un archivo PDB "a mano".

El formato PDB es texto plano en columnas fijas (heredado de las tarjetas
perforadas de 80 columnas). Cada línea empieza con un tipo de registro de 6
caracteres que dice qué hay en el resto de la línea:

  HEADER / TITLE / COMPND / SOURCE  -> metadatos de la entrada
  EXPDTA                            -> método experimental
  REMARK   2 RESOLUTION.            -> resolución (cristalografía)
  REMARK   3 ... R VALUE            -> calidad del refinamiento
  DBREF                             -> a qué entrada de UniProt/GenBank corresponde
  SEQRES                            -> secuencia completa de la cadena (3 letras)
  ATOM / HETATM                     -> coordenadas de cada átomo (proteína / no proteína)

En las líneas ATOM/HETATM cada dato ocupa columnas fijas:

  cols  7-11 número de átomo     cols 13-16 nombre del átomo (CA, N, O, ...)
  col  17    altloc              cols 18-20 residuo (MET, GLN, ...)
  col  22    cadena              cols 23-26 número de residuo
  cols 31-54 coordenadas X Y Z   cols 55-60 ocupación
  cols 61-66 factor B            cols 77-78 elemento

Este módulo se reutiliza en el resto del TP (cavidades, pLDDT).

Ejemplos:
    python explorar_pdb.py                       # resume datos/1UBQ.pdb
    python explorar_pdb.py datos/estructuras/1M14.pdb
"""

import argparse
import sys
from collections import Counter, OrderedDict
from pathlib import Path

AQUI = Path(__file__).resolve().parent
PDB_UBIQUITINA = AQUI / "datos" / "1UBQ.pdb"

TRES_A_UNA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # variantes que aparecen en algunos PDB
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def tres_a_una(residuos):
    """Convierte una lista de códigos de 3 letras a un string de 1 letra ('X' si no se conoce)."""
    return "".join(TRES_A_UNA.get(r.upper(), "X") for r in residuos)


def leer_atomos(ruta, cadena=None, solo_proteina=False, primer_modelo=True):
    """Devuelve una lista de dicts, uno por átomo (ATOM y HETATM).

    - cadena: si se indica, filtra por esa cadena.
    - solo_proteina: descarta los HETATM (ligandos, iones, aguas).
    - primer_modelo: en archivos con varios modelos (RMN) se queda con el primero.
    Las posiciones alternativas (altloc) se reducen a la primera ('' o 'A').
    """
    atomos = []
    with open(ruta) as f:
        for linea in f:
            registro = linea[:6].strip()
            if registro == "ENDMDL" and primer_modelo:
                break
            if registro not in ("ATOM", "HETATM"):
                continue
            if solo_proteina and registro == "HETATM":
                continue
            altloc = linea[16]
            if altloc not in (" ", "A"):
                continue
            if cadena is not None and linea[21] != cadena:
                continue
            atomos.append({
                "registro": registro,
                "serial": int(linea[6:11]),
                "nombre": linea[12:16].strip(),
                "residuo": linea[17:20].strip(),
                "cadena": linea[21],
                "num": int(linea[22:26]),
                "icode": linea[26].strip(),
                "x": float(linea[30:38]),
                "y": float(linea[38:46]),
                "z": float(linea[46:54]),
                "ocupacion": float(linea[54:60] or 0),
                "bfactor": float(linea[60:66] or 0),
                "elemento": linea[76:78].strip(),
            })
    return atomos


def residuos_ca(ruta, cadena=None):
    """OrderedDict {numero_residuo: nombre_3_letras} a partir de los carbonos alfa."""
    res = OrderedDict()
    for a in leer_atomos(ruta, cadena=cadena, solo_proteina=True):
        if a["nombre"] == "CA" and a["num"] not in res:
            res[a["num"]] = a["residuo"]
    return res


def secuencia_ca(ruta, cadena=None):
    """Secuencia (1 letra) de los residuos que efectivamente tienen coordenadas."""
    return tres_a_una(residuos_ca(ruta, cadena).values())


def secuencia_seqres(ruta, cadena=None):
    """Secuencia (1 letra) declarada en los registros SEQRES (incluye residuos no resueltos)."""
    por_cadena = OrderedDict()
    with open(ruta) as f:
        for linea in f:
            if linea.startswith("SEQRES"):
                c = linea[11]
                por_cadena.setdefault(c, []).extend(linea[19:].split())
    if cadena is None:
        cadena = next(iter(por_cadena))
    return tres_a_una(por_cadena.get(cadena, []))


def leer_encabezado(ruta):
    """Extrae metadatos básicos: id, título, método, resolución, R-work/R-free, DBREF."""
    info = {"titulo": "", "dbref": []}
    with open(ruta) as f:
        for linea in f:
            reg = linea[:6]
            if reg == "HEADER":
                info["clasificacion"] = linea[10:50].strip()
                info["fecha"] = linea[50:59].strip()
                info["id"] = linea[62:66].strip()
            elif reg == "TITLE ":
                info["titulo"] = (info["titulo"] + " " + linea[10:80].strip()).strip()
            elif reg == "EXPDTA":
                info["metodo"] = linea[10:79].strip()
            elif linea.startswith("REMARK   2 RESOLUTION."):
                valor = linea[22:].split()[0]
                try:
                    info["resolucion"] = float(valor)
                except ValueError:
                    info["resolucion"] = None
            elif linea.startswith("REMARK   3   R VALUE") and "WORKING" in linea and "r_work" not in info:
                try:
                    info["r_work"] = float(linea.split(":")[-1])
                except ValueError:
                    pass
            elif linea.startswith("REMARK   3   FREE R VALUE  ") and "r_free" not in info:
                try:
                    info["r_free"] = float(linea.split(":")[-1])
                except ValueError:
                    pass
            elif reg == "DBREF ":
                p = linea.split()
                # DBREF idcode chain seqBegin seqEnd db dbAccession dbIdCode dbseqBegin dbseqEnd
                info["dbref"].append({
                    "cadena": linea[12],
                    "pdb_inicio": int(linea[14:18]),
                    "pdb_fin": int(linea[20:24]),
                    "base": linea[26:32].strip(),
                    "acceso": linea[33:41].strip(),
                    "db_inicio": int(linea[55:60]),
                    "db_fin": int(linea[62:67]),
                })
            elif reg in ("ATOM  ", "HETATM"):
                break
    return info


def offset_uniprot(ruta, cadena):
    """Diferencia (numeración UniProt - numeración del PDB) según el DBREF de la cadena."""
    for d in leer_encabezado(ruta)["dbref"]:
        if d["cadena"] == cadena:
            return d["db_inicio"] - d["pdb_inicio"]
    raise ValueError(f"No hay DBREF para la cadena {cadena}")


TIPOS_HELICE = {1: "hélice α", 5: "hélice 3-10"}


def estructura_secundaria(ruta):
    """Lee los registros HELIX y SHEET: lista de (tipo, cadena, inicio, fin)."""
    elementos = []
    with open(ruta) as f:
        for linea in f:
            if linea.startswith("HELIX "):
                clase = int(linea[38:40] or 1)
                elementos.append((TIPOS_HELICE.get(clase, f"hélice clase {clase}"),
                                  linea[19], int(linea[21:25]), int(linea[33:37])))
            elif linea.startswith("SHEET "):
                elementos.append(("hebra β", linea[21], int(linea[22:26]), int(linea[33:37])))
    return sorted(elementos, key=lambda e: e[2])


def bfactor_por_residuo(ruta, cadena=None):
    """{num_residuo: factor B del CA}. En AlphaFold esa columna guarda el pLDDT."""
    return {a["num"]: a["bfactor"] for a in leer_atomos(ruta, cadena=cadena, solo_proteina=True)
            if a["nombre"] == "CA"}


def contar_registros(ruta):
    """Cuenta ATOM, HETATM y los HETATM que son agua (HOH)."""
    atomos = leer_atomos(ruta)
    c = Counter(a["registro"] for a in atomos)
    aguas = sum(1 for a in atomos if a["registro"] == "HETATM" and a["residuo"] == "HOH")
    elementos = Counter(a["elemento"] for a in atomos if a["registro"] == "ATOM")
    return {"ATOM": c["ATOM"], "HETATM": c["HETATM"], "aguas": aguas, "elementos": dict(elementos)}


def resumen(ruta):
    info = leer_encabezado(ruta)
    regs = contar_registros(ruta)
    cadenas = sorted({a["cadena"] for a in leer_atomos(ruta, solo_proteina=True)})
    seq_seqres = secuencia_seqres(ruta, cadenas[0])
    seq_ca = secuencia_ca(ruta, cadenas[0])
    lineas = [
        f"Entrada      : {info.get('id')}  ({info.get('clasificacion')}, {info.get('fecha')})",
        f"Título       : {info['titulo']}",
        f"Método       : {info.get('metodo')}",
        f"Resolución   : {info.get('resolucion')} Å",
        f"R-work/R-free: {info.get('r_work')} / {info.get('r_free')}",
        f"DBREF        : " + "; ".join(
            f"cadena {d['cadena']} {d['pdb_inicio']}-{d['pdb_fin']} -> {d['base']} {d['acceso']} {d['db_inicio']}-{d['db_fin']}"
            for d in info["dbref"]),
        f"Cadenas      : {', '.join(cadenas)}",
        f"Átomos       : {regs['ATOM']} ATOM, {regs['HETATM']} HETATM ({regs['aguas']} de agua)",
        f"Elementos    : {regs['elementos']}",
        f"SEQRES ({cadenas[0]})  : {len(seq_seqres)} aa  {seq_seqres}",
        f"Desde CA ({cadenas[0]}): {len(seq_ca)} aa  {seq_ca}",
    ]
    ss = estructura_secundaria(ruta)
    if ss:
        lineas.append("Estructura 2ria: " + "; ".join(f"{t} {i}-{f}" for t, _, i, f in ss))
    b = bfactor_por_residuo(ruta, cadenas[0])
    if b:
        mas_movil = max(b, key=b.get)
        lineas.append(f"Factor B (CA)  : medio {sum(b.values()) / len(b):.1f} Å², "
                      f"máximo {b[mas_movil]:.1f} Å² en el residuo {mas_movil}")
    return "\n".join(lineas)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Resumen de un archivo PDB.")
    parser.add_argument("pdb", nargs="?", default=str(PDB_UBIQUITINA))
    args = parser.parse_args(argv)
    print(resumen(args.pdb))
    return 0


if __name__ == "__main__":
    sys.exit(main())
