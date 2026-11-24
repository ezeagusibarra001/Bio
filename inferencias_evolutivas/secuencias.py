#!/usr/bin/env python3
"""Lectura/escritura de FASTA y armado de los sets de secuencias del TP.

Los encabezados del archivo del curso (`SecuenciasCytocromoC.fasta`) traen
comas, corchetes y paréntesis (p. ej. `..._proximal,_isoform_B_[Drosophila...]`)
que rompen el formato Newick de los árboles. Por eso cada secuencia se
renombra con una etiqueta corta (el organismo) y se guarda la correspondencia.

Sets que se arman:
  - set3: humano, mosca y bartmosca (Desafío I/II, "dadas las secuencias de
    Mosca, humano y Moscahumano").
  - set9: las 9 secuencias del archivo (Desafío I, "el resto de las secuencias").

  - set4: set3 + el citocromo c *distal* de Drosophila (NP_477164.1), la otra
    copia (parálogo) del gen en la mosca. No viene en el archivo del curso
    (que trae la isoforma *proximal*); se baja de NCBI con --descargar-distal.

Uso:
    python secuencias.py                      # regenera set3, set9 (y set4 si ya está el distal)
    python secuencias.py --descargar-distal   # baja NP_477164.1 de NCBI y arma set4
"""

import re
import ssl
import sys
import urllib.request
from pathlib import Path

AQUI = Path(__file__).resolve().parent
DATOS = AQUI / "datos"
ORIGINAL = DATOS / "SecuenciasCytocromoC.fasta"
DISTAL = DATOS / "NP_477164.1.fasta"
ETIQUETA_DISTAL = "Drosophila_distal"

SET3 = ["Homo_sapiens", "Drosophila_melanogaster", "bartmosca"]
SET4 = ["Homo_sapiens", "Drosophila_melanogaster", ETIQUETA_DISTAL, "bartmosca"]


def descargar_ncbi(accesion, destino):
    """Baja una proteína de NCBI (efetch, formato FASTA) y la guarda tal cual."""
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=protein&id={accesion}&rettype=fasta&retmode=text&tool=tp_bioinfo")
    sistema = Path("/etc/ssl/cert.pem")
    contexto = ssl.create_default_context(cafile=str(sistema) if sistema.exists() else None)
    req = urllib.request.Request(url, headers={"User-Agent": "tp-bioinfo-client"})
    with urllib.request.urlopen(req, timeout=60, context=contexto) as r:
        texto = r.read().decode()
    if not texto.startswith(">"):
        raise RuntimeError(f"Respuesta inesperada de NCBI para {accesion}: {texto[:200]}")
    Path(destino).write_text(texto)
    return texto


def leer_fasta(ruta):
    """Devuelve una lista de tuplas (encabezado, secuencia), en el orden del archivo."""
    registros = []
    for linea in Path(ruta).read_text().splitlines():
        linea = linea.strip()
        if not linea:
            continue
        if linea.startswith(">"):
            registros.append([linea[1:].strip(), ""])
        elif registros:
            registros[-1][1] += linea.replace(" ", "").upper()
    return [tuple(r) for r in registros]


def escribir_fasta(registros, ruta, ancho=60):
    with open(ruta, "w") as f:
        for nombre, sec in registros:
            f.write(f">{nombre}\n")
            for i in range(0, len(sec), ancho):
                f.write(sec[i:i + ancho] + "\n")


def etiqueta_corta(encabezado):
    """`NP_061820.1_cytochrome_c_[Homo_sapiens]` -> `Homo_sapiens`.

    Toma lo que está entre corchetes (el organismo) y deja solo letras, números
    y guiones bajos. Si no hay corchetes (caso `bartmosca`) usa el encabezado.
    """
    m = re.search(r"\[([^\]]+)\]", encabezado)
    base = m.group(1) if m else encabezado
    base = re.sub(r"[^A-Za-z0-9_]+", "_", base).strip("_")
    return base


def armar_sets():
    originales = leer_fasta(ORIGINAL)
    renombrados = [(etiqueta_corta(h), s) for h, s in originales]
    etiquetas = [e for e, _ in renombrados]
    if len(set(etiquetas)) != len(etiquetas):
        raise ValueError(f"Etiquetas repetidas: {etiquetas}")

    escribir_fasta(renombrados, DATOS / "set9.fasta")
    por_nombre = dict(renombrados)
    escribir_fasta([(e, por_nombre[e]) for e in SET3], DATOS / "set3.fasta")
    if DISTAL.exists():
        (_, sec_distal), = leer_fasta(DISTAL)
        por_nombre[ETIQUETA_DISTAL] = sec_distal
        escribir_fasta([(e, por_nombre[e]) for e in SET4], DATOS / "set4.fasta")

    with open(DATOS / "etiquetas.tsv", "w") as f:
        f.write("etiqueta\tencabezado_original\tlargo\n")
        for (h, s), (e, _) in zip(originales, renombrados):
            f.write(f"{e}\t{h}\t{len(s)}\n")
    return renombrados


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--descargar-distal" in argv:
        descargar_ncbi("NP_477164.1", DISTAL)
        print(f"Bajado {DISTAL.name} de NCBI")
    for etiqueta, sec in armar_sets():
        print(f"{etiqueta:<26} {len(sec):>4} aa")
    print(f"\nEscritos: {DATOS / 'set3.fasta'}, {DATOS / 'set9.fasta'}, {DATOS / 'etiquetas.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
