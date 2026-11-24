#!/usr/bin/env python3
"""DESAFÍO II (control) - Alineamiento múltiple con MAFFT local.

Se usa como segundo alineador para comprobar que el resultado no depende de
la herramienta (y como respaldo si el servicio del EBI no responde).
Modo L-INS-i (`--localpair --maxiterate 1000`), el más preciso de MAFFT para
pocas secuencias cortas.

  - `resultados/aln_mafft_<set>.fasta`
  - `resultados/aln_mafft_<set>.clustal`   (mismo alineamiento, formato Clustal)

Requiere MAFFT instalado (`brew install mafft` / `conda install -c bioconda mafft`).

Uso:
    python mafft_local.py set3 set9
"""

import shutil
import subprocess
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
DATOS = AQUI / "datos"
RESULTADOS = AQUI / "resultados"
OPCIONES = ["--localpair", "--maxiterate", "1000", "--quiet"]


def correr_mafft(entrada, extra=()):
    if shutil.which("mafft") is None:
        raise SystemExit("No se encontró `mafft` en el PATH (brew install mafft).")
    salida = subprocess.run(["mafft", *OPCIONES, *extra, str(entrada)],
                            check=True, capture_output=True, text=True)
    return salida.stdout


def main(argv=None):
    sets = (argv if argv is not None else sys.argv[1:]) or ["set3", "set9"]
    RESULTADOS.mkdir(exist_ok=True)
    for nombre in sets:
        entrada = DATOS / f"{nombre}.fasta"
        (RESULTADOS / f"aln_mafft_{nombre}.fasta").write_text(correr_mafft(entrada))
        (RESULTADOS / f"aln_mafft_{nombre}.clustal").write_text(correr_mafft(entrada, ["--clustalout"]))
        print(f"[{nombre}] resultados/aln_mafft_{nombre}.fasta y .clustal")
    return 0


if __name__ == "__main__":
    sys.exit(main())
