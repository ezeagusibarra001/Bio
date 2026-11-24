#!/usr/bin/env python3
"""DESAFÍO II - Alineamiento múltiple con Clustal Omega (servicio REST del EBI).

Hace lo mismo que la página https://www.ebi.ac.uk/Tools/msa/clustalo/ pero
desde la terminal: sube el FASTA, espera a que el trabajo termine y baja los
resultados.

  - `resultados/aln_<set>.clustal_num`  alineamiento en formato Clustal (con numeración)
  - `resultados/aln_<set>.fasta`        el mismo alineamiento en FASTA (lo usa IQ-TREE)
  - `resultados/pim_<set>.txt`          matriz de identidad que calcula Clustal
  - `resultados/ebi_jobs.tsv`           id del trabajo en el EBI, para poder rastrearlo

El EBI exige un email en cada pedido; se usa uno genérico.

Uso:
    python clustal_ebi.py set3 set9
"""

import ssl
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

AQUI = Path(__file__).resolve().parent
DATOS = AQUI / "datos"
RESULTADOS = AQUI / "resultados"
BASE = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
EMAIL = "anonymous@example.org"


def _contexto_ssl():
    """El Python de python.org en macOS no trae certificados raíz: si existe el
    bundle del sistema se usa ese (la verificación sigue activa)."""
    sistema = Path("/etc/ssl/cert.pem")
    return ssl.create_default_context(cafile=str(sistema) if sistema.exists() else None)


def _pedir(url, datos=None):
    req = urllib.request.Request(url, data=datos, headers={"User-Agent": "tp-bioinfo-client"})
    with urllib.request.urlopen(req, timeout=60, context=_contexto_ssl()) as r:
        return r.read().decode()


def enviar(fasta_texto, outfmt="clustal_num"):
    datos = urllib.parse.urlencode({
        "email": EMAIL,
        "sequence": fasta_texto,
        "stype": "protein",
        "outfmt": outfmt,
    }).encode()
    return _pedir(f"{BASE}/run", datos).strip()


def esperar(job_id, intervalo=5, maximo=600):
    inicio = time.time()
    while time.time() - inicio < maximo:
        estado = _pedir(f"{BASE}/status/{job_id}").strip()
        if estado == "FINISHED":
            return
        if estado not in ("RUNNING", "QUEUED", "PENDING"):
            raise RuntimeError(f"El trabajo {job_id} terminó en estado {estado}")
        time.sleep(intervalo)
    raise TimeoutError(f"El trabajo {job_id} no terminó en {maximo} s")


def resultado(job_id, tipo):
    return _pedir(f"{BASE}/result/{job_id}/{tipo}")


def alinear_set(nombre):
    fasta = (DATOS / f"{nombre}.fasta").read_text()
    RESULTADOS.mkdir(exist_ok=True)

    job = enviar(fasta, "clustal_num")
    print(f"[{nombre}] trabajo {job} enviado...")
    esperar(job)
    (RESULTADOS / f"aln_{nombre}.clustal_num").write_text(resultado(job, "aln-clustal_num"))
    (RESULTADOS / f"pim_{nombre}.txt").write_text(resultado(job, "pim"))

    # Segundo pedido con salida FASTA (la entrada de IQ-TREE).
    job_fa = enviar(fasta, "fa")
    esperar(job_fa)
    (RESULTADOS / f"aln_{nombre}.fasta").write_text(resultado(job_fa, "aln-fasta"))

    with open(RESULTADOS / "ebi_jobs.tsv", "a") as f:
        f.write(f"{nombre}\t{job}\t{job_fa}\t{time.strftime('%Y-%m-%d')}\n")
    print(f"[{nombre}] listo: resultados/aln_{nombre}.clustal_num y .fasta")


def main(argv=None):
    sets = (argv if argv is not None else sys.argv[1:]) or ["set3", "set9"]
    for nombre in sets:
        alinear_set(nombre)
    return 0


if __name__ == "__main__":
    sys.exit(main())
