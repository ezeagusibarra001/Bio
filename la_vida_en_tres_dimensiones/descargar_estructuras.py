#!/usr/bin/env python3
"""Baja (o vuelve a bajar) todos los datos del TP a datos/.

- Estructuras del PDB (RCSB):     1UBQ, 1M14, 3W32, 2FVL, 1AO6, 4JK4
- Entradas de UniProt (JSON):     P17516 (AKR1C4), P00533 (EGFR), P02768 (HSA),
                                  P02769 (BSA), P03129 (E7 de HPV16)
- Modelos de AlphaFold DB de E7:  PDB + pLDDT + PAE (ambos modelos que lista la API)

Los archivos ya están en el repo; esto es solo para reproducir.

    python descargar_estructuras.py
"""

import gzip
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

AQUI = Path(__file__).resolve().parent
DATOS = AQUI / "datos"
UA = {"User-Agent": "Mozilla/5.0 (tp-bioinformatica)"}

PDBS = ["1M14", "3W32", "2FVL", "1AO6", "4JK4"]
UNIPROT = ["P17516", "P00533", "P02768", "P02769", "P03129"]


def _leer_url(url):
    """urllib primero; si el Python no tiene certificados raíz (pasa con el instalador
    de python.org en macOS), usa curl, que trae los del sistema."""
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read()
    except urllib.error.URLError as e:
        if "CERTIFICATE_VERIFY_FAILED" not in str(e):
            raise
        return subprocess.run(["curl", "-sSfL", "-m", "120", "-A", UA["User-Agent"], url],
                              check=True, capture_output=True).stdout


def bajar(url, destino):
    datos = _leer_url(url)
    if datos[:2] == b"\x1f\x8b":          # algunos archivos de AlphaFold DB vienen gzipeados
        datos = gzip.decompress(datos)
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_bytes(datos)
    print(f"  {destino.relative_to(AQUI)}  ({len(datos):,} bytes)")


def main():
    bajar("https://files.rcsb.org/download/1UBQ.pdb", DATOS / "1UBQ.pdb")
    for p in PDBS:
        bajar(f"https://files.rcsb.org/download/{p}.pdb", DATOS / "estructuras" / f"{p}.pdb")
    for acc in UNIPROT:
        bajar(f"https://rest.uniprot.org/uniprotkb/{acc}.json", DATOS / "uniprot" / f"{acc}.json")

    api = DATOS / "alphafold" / "api_P03129.json"
    bajar("https://alphafold.ebi.ac.uk/api/prediction/P03129", api)
    for entrada in json.loads(api.read_text()):
        eid = entrada["entryId"]
        for url in (entrada["pdbUrl"], entrada["plddtDocUrl"], entrada["paeDocUrl"]):
            bajar(url, DATOS / "alphafold" / f"{eid}-{url.rsplit(eid + '-', 1)[-1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
