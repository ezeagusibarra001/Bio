#!/usr/bin/env python3
"""Descarga de datos del TP (necesita internet; el análisis NO lo necesita).

1. Baja las 9 secuencias de citocromo c de la tabla del enunciado desde NCBI
   Protein (E-utilities, efetch) y las guarda en datos/citocromo_c.fasta.
2. Opcional (--clustalo): manda esas secuencias a Clustal Omega en el EBI
   (REST de Job Dispatcher), espera el resultado y guarda el alineamiento
   (formato clustal_num) y los árboles que devuelve el servidor.

El EBI exige un email en cada job; se usa uno genérico (se puede cambiar con
--email). Respetá su política de uso: un job a la vez y nada de reintentos en
loop.

Ejemplos:
    python descargar_datos.py
    python descargar_datos.py --clustalo
"""

import argparse
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

CARPETA = Path(__file__).resolve().parent
DATOS = CARPETA / "datos"

ACCESIONES = [
    "NP_061820.1",     # Homo sapiens
    "NP_001072946.1",  # Gallus gallus
    "NP_001065289.1",  # Pan troglodytes
    "NP_001157486.1",  # Equus caballus
    "NP_001183974.1",  # Canis lupus familiaris
    "AEP27192.1",      # Gorilla gorilla
    "XP_024245566.1",  # Oncorhynchus tshawytscha
    "NP_001086101.1",  # Xenopus laevis
    "NP_477164.1",     # Drosophila melanogaster
]

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CLUSTALO = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


def url_efetch(accesiones):
    params = {"db": "protein", "id": ",".join(accesiones), "rettype": "fasta",
              "retmode": "text", "tool": "generic"}
    return f"{EFETCH}?{urllib.parse.urlencode(params)}"


def descargar_fasta(destino=DATOS / "citocromo_c.fasta"):
    with urllib.request.urlopen(url_efetch(ACCESIONES), timeout=60) as resp:
        texto = resp.read().decode()
    texto = "\n".join(l for l in texto.splitlines() if l.strip()) + "\n"
    destino.write_text(texto)
    print(f"{texto.count('>')} secuencias -> {destino}")
    return destino


def _get(url):
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode()


def correr_clustalo(fasta, email, prefijo=DATOS / "citocromo_c_clustalo"):
    datos = urllib.parse.urlencode({
        "email": email, "stype": "protein", "outfmt": "clustal_num",
        "title": "citocromo_c", "sequence": Path(fasta).read_text(),
    }).encode()
    with urllib.request.urlopen(f"{CLUSTALO}/run", data=datos, timeout=60) as resp:
        job = resp.read().decode().strip()
    print("job:", job)
    while (estado := _get(f"{CLUSTALO}/status/{job}").strip()) in ("RUNNING", "QUEUED"):
        time.sleep(10)
    if estado != "FINISHED":
        raise RuntimeError(f"Clustal Omega terminó con estado {estado}")
    tipos = {"aln-clustal_num": ".aln", "pim": ".pim", "tree": ".guide.dnd", "phylotree": ".phylo.ph"}
    disponibles = _get(f"{CLUSTALO}/resulttypes/{job}")
    for tipo, ext in tipos.items():
        if f"<identifier>{tipo}</identifier>" in disponibles:
            Path(f"{prefijo}{ext}").write_text(_get(f"{CLUSTALO}/result/{job}/{tipo}"))
            print(f"{tipo} -> {prefijo}{ext}")
    return job


def main(argv=None):
    parser = argparse.ArgumentParser(description="Descarga las secuencias y (opcional) corre Clustal Omega.")
    parser.add_argument("--clustalo", action="store_true", help="Además, alinear con Clustal Omega en el EBI")
    parser.add_argument("--email", default="anonymous@example.org", help="Email que exige el EBI")
    args = parser.parse_args(argv)
    DATOS.mkdir(exist_ok=True)
    fasta = descargar_fasta()
    if args.clustalo:
        correr_clustalo(fasta, args.email)
    return 0


if __name__ == "__main__":
    sys.exit(main())
