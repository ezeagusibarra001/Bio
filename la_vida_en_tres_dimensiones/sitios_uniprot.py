#!/usr/bin/env python3
"""DESAFÍO II - Sitios biológicamente relevantes de una proteína según UniProt.

Lee la entrada de UniProt en JSON (bajada con la API REST,
https://rest.uniprot.org/uniprotkb/<ACC>.json, y guardada en datos/uniprot/)
y lista las anotaciones de posición: sitio activo, sitios de unión, "Site",
variantes naturales. Si se le pasa un PDB, chequea qué residuo hay en esa
posición en la estructura (usando el offset del DBREF).

Ejemplos:
    python sitios_uniprot.py                          # AKR1C4 (P17516)
    python sitios_uniprot.py P17516 --pdb datos/estructuras/2FVL.pdb --cadena B
    python sitios_uniprot.py P00533 --tipos "Active site" "Binding site"
"""

import argparse
import json
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

import explorar_pdb as ep  # noqa: E402

DIR_UNIPROT = AQUI / "datos" / "uniprot"
TIPOS_POR_DEFECTO = ("Active site", "Binding site", "Site", "Natural variant")


def leer_entrada(acceso):
    with open(DIR_UNIPROT / f"{acceso}.json") as f:
        return json.load(f)


def sitios(entrada, tipos=TIPOS_POR_DEFECTO):
    """Lista de dicts con las features de los tipos pedidos (posiciones en numeración UniProt)."""
    seq = entrada["sequence"]["value"]
    salida = []
    for f in entrada["features"]:
        if f["type"] not in tipos:
            continue
        ini = f["location"]["start"]["value"]
        fin = f["location"]["end"]["value"]
        alt = f.get("alternativeSequence", {})
        salida.append({
            "tipo": f["type"],
            "inicio": ini,
            "fin": fin,
            "residuos": seq[ini - 1:fin],
            "descripcion": f.get("description", ""),
            "ligando": f.get("ligand", {}).get("name", ""),
            "variante": (alt.get("originalSequence", "") + ">" + "/".join(alt.get("alternativeSequences", []))) if alt else "",
        })
    return salida


def posiciones(entrada, tipos=("Active site", "Binding site", "Site")):
    """Conjunto de posiciones UniProt cubiertas por las features pedidas."""
    pos = set()
    for s in sitios(entrada, tipos):
        pos.update(range(s["inicio"], s["fin"] + 1))
    return pos


def residuo_en_pdb(ruta_pdb, cadena, pos_uniprot):
    """Qué residuo (3 letras) hay en el PDB en la posición UniProt dada, y su número en el PDB."""
    off = ep.offset_uniprot(ruta_pdb, cadena)
    num_pdb = pos_uniprot - off
    return num_pdb, ep.residuos_ca(ruta_pdb, cadena).get(num_pdb)


def main(argv=None):
    p = argparse.ArgumentParser(description="Sitios anotados en UniProt.")
    p.add_argument("acceso", nargs="?", default="P17516")
    p.add_argument("--tipos", nargs="+", default=list(TIPOS_POR_DEFECTO))
    p.add_argument("--pdb", help="PDB para chequear los residuos en la estructura")
    p.add_argument("--cadena", default="A")
    a = p.parse_args(argv)

    e = leer_entrada(a.acceso)
    nombre = e["proteinDescription"]["recommendedName"]["fullName"]["value"]
    gen = e["genes"][0]["geneName"]["value"]
    print(f"{a.acceso} | {gen} | {nombre} | {e['sequence']['length']} aa | {e['organism']['scientificName']}")
    for s in sitios(e, a.tipos):
        rango = f"{s['inicio']}" if s["inicio"] == s["fin"] else f"{s['inicio']}-{s['fin']}"
        extra = " ".join(x for x in (s["descripcion"], s["ligando"], s["variante"]) if x)
        linea = f"  {s['tipo']:<15} {rango:<8} {s['residuos']:<12} {extra}"
        if a.pdb and s["tipo"] != "Natural variant":
            nums = [residuo_en_pdb(a.pdb, a.cadena, i) for i in range(s["inicio"], s["fin"] + 1)]
            linea += "   | PDB: " + " ".join(f"{ep.TRES_A_UNA.get(r, '?') if r else '-'}{n}" for n, r in nums)
        print(linea)
    return 0


if __name__ == "__main__":
    sys.exit(main())
