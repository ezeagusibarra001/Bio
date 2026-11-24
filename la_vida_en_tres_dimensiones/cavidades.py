#!/usr/bin/env python3
"""DESAFÍOS I, III, IV y V - Cavidades proteicas y residuos activados (pKa anómalos).

El TP pide usar CaviDB (https://www.cavidb.org). Según su paper
(Velez Rueda et al. 2023, Database, doi:10.1093/database/baad010), CaviDB:
  - predice cavidades con Fpocket (parámetros por defecto),
  - marca una cavidad como "drogable" si el drug score de Fpocket es > 0.5,
  - estima el pKa de los residuos ionizables con PROPKA y, por cavidad, el
    "net pKa shift" = suma de |pKa predicho - pKa modelo| de sus residuos.

CaviDB no respondía desde esta máquina (timeout) y Fpocket no se pudo instalar,
así que este script reproduce la misma lógica con herramientas instalables por pip:
  - cavidades: pyKVFinder (Guerra et al. 2021, BMC Bioinformatics), que da
    volumen, área, profundidad e hidropatía de cada cavidad y los residuos que la tapizan;
  - pKa: PROPKA 3 (la misma herramienta que usa CaviDB).
Los números NO son idénticos a los de CaviDB (otro detector, otros parámetros),
pero las preguntas del TP se pueden responder igual y los links a CaviDB quedan
en el README para comparar.

Flujo por estructura:
  1. extraer la cadena pedida, solo proteína (sin ligandos ni aguas), altloc A;
  2. pyKVFinder -> tabla de cavidades;
  3. PROPKA -> pKa por residuo (el .pka se guarda en resultados/propka/);
  4. cruzar: pKa medio, rango y net shift por cavidad, residuos activados,
     solapamiento con los residuos de interés y con el ligando cristalográfico.

Uso:
    python cavidades.py                 # corre todo (necesita pyKVFinder y propka)
    python cavidades.py --solo 1AO6_A   # una sola estructura
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("KMP_WARNINGS", "0")  # silencia avisos de OpenMP de pyKVFinder

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

import explorar_pdb as ep  # noqa: E402

DIR_ESTRUCTURAS = AQUI / "datos" / "estructuras"
DIR_CADENAS = AQUI / "datos" / "cadenas"
DIR_RESULTADOS = AQUI / "resultados"
DIR_PROPKA = DIR_RESULTADOS / "propka"

# pKa "modelo" (residuo libre en agua) que usa PROPKA 3
PKA_MODELO = {"ASP": 3.80, "GLU": 4.50, "HIS": 6.50, "CYS": 9.00,
              "TYR": 10.00, "LYS": 10.50, "ARG": 12.50, "N+": 8.00, "C-": 3.20}
UMBRAL_ACTIVADO = 1.0  # |pKa - pKa modelo| >= 1 unidad de pH -> lo considero "activado"

# Estructuras del TP. Los residuos de interés van en numeración del PDB
# (ver numeracion.py para la conversión desde UniProt).
ESTRUCTURAS = {
    "1M14_A": {"pdb": "1M14", "cadena": "A", "proteina": "EGFR (confórmero activo)",
               "interes": {813: "Asp837 UniProt (HRD, aceptor de protones)",
                           721: "Lys745 UniProt (une ATP)",
                           831: "Asp855 UniProt (DFG, une ATP)",
                           837: "residuo 837 *del PDB* (= UniProt 861)"},
               "ligandos": []},
    "3W32_A": {"pdb": "3W32", "cadena": "A", "proteina": "EGFR (confórmero inactivo)",
               "interes": {837: "Asp837 UniProt (HRD, aceptor de protones)",
                           745: "Lys745 UniProt (une ATP)",
                           855: "Asp855 UniProt (DFG, une ATP)"},
               "ligandos": ["W32"]},
    "2FVL_B": {"pdb": "2FVL", "cadena": "B", "proteina": "AKR1C4",
               "interes": {50: "Asp50 (une NADP)", 55: "Tyr55 (sitio activo)",
                           84: "Lys84 (baja el pKa de Tyr55)", 117: "His117 (tétrada AKR)"},
               "ligandos": ["NAP"]},
    "1AO6_A": {"pdb": "1AO6", "cadena": "A", "proteina": "albúmina humana (HSA)",
               "interes": {34: "Cys34", 195: "Lys195", 199: "Lys199", 222: "Arg222",
                           410: "Arg410", 411: "Tyr411"},
               "ligandos": []},
    "4JK4_A": {"pdb": "4JK4", "cadena": "A", "proteina": "albúmina bovina (BSA)",
               "interes": {221: "Lys221 (equivale a Arg222 de HSA)",
                           198: "Arg198 (equivale a Lys199 de HSA)",
                           194: "Arg194 (equivale a Lys195 de HSA)",
                           409: "Arg409 (equivale a Arg410 de HSA)",
                           410: "Tyr410 (equivale a Tyr411 de HSA)",
                           34: "Cys34"},
               "ligandos": ["DIU"]},
}


# ---------------------------------------------------------------- preparación
def extraer_cadena(origen, cadena, destino):
    """Escribe solo los ATOM de una cadena (primer modelo, altloc ''/A)."""
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(origen) as f, open(destino, "w") as g:
        for linea in f:
            if linea.startswith("ENDMDL"):
                break
            if linea.startswith("ATOM") and linea[21] == cadena and linea[16] in " A":
                g.write(linea[:16] + " " + linea[17:])
                n += 1
        g.write("END\n")
    return n


def residuos_cerca_de_ligando(ruta_pdb, cadena, ligandos, corte=4.0):
    """Residuos de la cadena con algún átomo a <= corte Å de un átomo pesado del ligando
    (solo ligandos asignados a esa misma cadena en el PDB)."""
    if not ligandos:
        return set()
    atomos = ep.leer_atomos(ruta_pdb)
    lig = [(a["x"], a["y"], a["z"]) for a in atomos
           if a["registro"] == "HETATM" and a["residuo"] in ligandos
           and a["cadena"] == cadena and a["elemento"] != "H"]
    prot = [a for a in atomos if a["registro"] == "ATOM" and a["cadena"] == cadena]
    c2 = corte * corte
    cerca = set()
    for a in prot:
        for (x, y, z) in lig:
            if (a["x"] - x) ** 2 + (a["y"] - y) ** 2 + (a["z"] - z) ** 2 <= c2:
                cerca.add(a["num"])
                break
    return cerca


# ---------------------------------------------------------------- cavidades
def detectar_cavidades(ruta_cadena):
    """Corre pyKVFinder (parámetros por defecto) y devuelve una lista de cavidades."""
    import pyKVFinder  # import diferido: los tests no lo necesitan
    r = pyKVFinder.run_workflow(str(ruta_cadena), include_depth=True, include_hydropathy=True)
    if r is None:
        return []
    cavs = []
    for cid in sorted(r.volume, key=lambda c: -r.volume[c]):
        residuos = sorted({(int(num), nombre) for num, _, nombre in r.residues[cid]})
        cavs.append({
            "id_kv": cid,
            "volumen": r.volume[cid],
            "area": r.area[cid],
            "prof_max": r.max_depth[cid],
            "prof_media": r.avg_depth[cid],
            "hidropatia": r.avg_hydropathy[cid],
            "residuos": residuos,
        })
    for i, c in enumerate(cavs, start=1):   # numeradas de mayor a menor volumen, como en CaviDB
        c["cavidad"] = i
    return cavs


# ---------------------------------------------------------------- PROPKA
def correr_propka(ruta_cadena, dir_salida=DIR_PROPKA):
    """Corre propka3 sobre la cadena y devuelve la ruta del .pka generado."""
    dir_salida = Path(dir_salida)
    dir_salida.mkdir(parents=True, exist_ok=True)
    propka = Path(sys.executable).with_name("propka3")
    subprocess.run([str(propka), "--quiet", str(Path(ruta_cadena).resolve())],
                   cwd=dir_salida, check=True, capture_output=True)
    return dir_salida / (Path(ruta_cadena).stem + ".pka")


def leer_propka(ruta_pka):
    """Parsea un .pka de PROPKA 3.

    Devuelve {num_residuo: {"residuo", "cadena", "pka", "modelo", "shift", "enterrado"}}
    para los grupos ionizables de la proteína (se ignoran N+/C- terminales y las
    Cys de puentes disulfuro, que PROPKA marca con pKa 99.99).
    'enterrado' es el % de burial que PROPKA informa en la tabla de determinantes.
    """
    lineas = Path(ruta_pka).read_text().splitlines()
    enterrado = {}
    for l in lineas:
        p = l.split()
        # tabla de determinantes: "ASP 813 A   4.68    71 %  ..."
        if len(p) >= 6 and p[0] in PKA_MODELO and p[5] == "%" and p[1].lstrip("-").isdigit():
            enterrado.setdefault(int(p[1]), int(p[4]))
    res = {}
    try:
        ini = next(i for i, l in enumerate(lineas) if "SUMMARY OF THIS PREDICTION" in l)
    except StopIteration:
        raise ValueError(f"{ruta_pka} no tiene sección SUMMARY")
    for l in lineas[ini + 2:]:
        if l.startswith("---") or not l.strip():
            break
        p = l.split()
        if len(p) < 5 or p[0] not in PKA_MODELO or p[0] in ("N+", "C-"):
            continue
        num = int(p[1])
        pka, modelo = float(p[3]), float(p[4])
        if pka >= 99:   # PROPKA pone 99.99 a las Cys en puentes disulfuro: no son ionizables
            continue
        res[num] = {"residuo": p[0], "cadena": p[2], "pka": pka, "modelo": modelo,
                    "shift": round(pka - modelo, 2), "enterrado": enterrado.get(num)}
    return res


def anotar_cavidades(cavidades, pkas, ligando=frozenset()):
    """Agrega a cada cavidad: residuos ionizables, pKa medio/rango, net shift, activados."""
    for c in cavidades:
        nums = [n for n, _ in c["residuos"]]
        ion = [(n, pkas[n]) for n in nums if n in pkas]
        c["n_residuos"] = len(nums)
        c["n_ionizables"] = len(ion)
        c["pka_medio"] = round(sum(d["pka"] for _, d in ion) / len(ion), 2) if ion else None
        c["pka_min"] = min((d["pka"] for _, d in ion), default=None)
        c["pka_max"] = max((d["pka"] for _, d in ion), default=None)
        c["net_shift"] = round(sum(abs(d["shift"]) for _, d in ion), 2)
        c["activados"] = [f"{d['residuo']}{n}({d['pka']:.2f})" for n, d in ion
                          if abs(d["shift"]) >= UMBRAL_ACTIVADO]
        c["contacto_ligando"] = len(set(nums) & set(ligando))
    return cavidades


def cavidades_de_residuo(cavidades, num):
    return [c["cavidad"] for c in cavidades if any(n == num for n, _ in c["residuos"])]


# ---------------------------------------------------------------- salida
COLUMNAS = ["cavidad", "volumen", "area", "prof_max", "prof_media", "hidropatia",
            "n_residuos", "n_ionizables", "pka_medio", "pka_min", "pka_max", "net_shift",
            "contacto_ligando", "activados", "residuos"]


def guardar_cavidades(cavidades, ruta):
    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNAS)
        for c in cavidades:
            fila = []
            for k in COLUMNAS:
                v = c.get(k)
                if k == "residuos":
                    v = " ".join(f"{ep.TRES_A_UNA.get(r, 'X')}{n}" for n, r in v)
                elif k == "activados":
                    v = " ".join(v)
                fila.append(v)
            w.writerow(fila)


def leer_cavidades_csv(ruta):
    """Lee un CSV generado por guardar_cavidades (para tests y figuras)."""
    with open(ruta) as f:
        filas = list(csv.DictReader(f))
    for r in filas:
        for k in ("cavidad", "n_residuos", "n_ionizables", "contacto_ligando"):
            r[k] = int(r[k])
        for k in ("volumen", "area", "prof_max", "prof_media", "hidropatia", "net_shift"):
            r[k] = float(r[k])
        r["residuos"] = r["residuos"].split()
        r["activados"] = r["activados"].split()
    return filas


def procesar(clave, usar_pka_guardado=False):
    cfg = ESTRUCTURAS[clave]
    origen = DIR_ESTRUCTURAS / f"{cfg['pdb']}.pdb"
    cadena_pdb = DIR_CADENAS / f"{clave}.pdb"
    extraer_cadena(origen, cfg["cadena"], cadena_pdb)

    cavs = detectar_cavidades(cadena_pdb)
    pka_path = DIR_PROPKA / f"{clave}.pka"
    if not (usar_pka_guardado and pka_path.exists()):
        pka_path = correr_propka(cadena_pdb)
    pkas = leer_propka(pka_path)
    lig = residuos_cerca_de_ligando(origen, cfg["cadena"], cfg["ligandos"])
    anotar_cavidades(cavs, pkas, lig)
    guardar_cavidades(cavs, DIR_RESULTADOS / f"cavidades_{clave}.csv")

    interes = []
    for num, desc in cfg["interes"].items():
        d = pkas.get(num, {})
        interes.append({"estructura": clave, "residuo_pdb": num,
                        "aa": ep.residuos_ca(cadena_pdb, cfg["cadena"]).get(num),
                        "descripcion": desc, "pka": d.get("pka"), "modelo": d.get("modelo"),
                        "shift": d.get("shift"), "enterrado": d.get("enterrado"),
                        "cavidades": " ".join(map(str, cavidades_de_residuo(cavs, num)))})
    resumen = {"estructura": clave, "proteina": cfg["proteina"], "n_cavidades": len(cavs),
               "vol_mayor": round(cavs[0]["volumen"], 1) if cavs else None,
               "vol_total": round(sum(c["volumen"] for c in cavs), 1),
               "residuos_ligando": len(lig),
               "cav_con_ligando": " ".join(str(c["cavidad"]) for c in cavs if c["contacto_ligando"])}
    return resumen, interes, cavs


def main(argv=None):
    p = argparse.ArgumentParser(description="Cavidades (pyKVFinder) + pKa (PROPKA).")
    p.add_argument("--solo", choices=sorted(ESTRUCTURAS), help="Procesar una sola estructura")
    p.add_argument("--usar-pka-guardado", action="store_true",
                   help="No vuelve a correr PROPKA si ya está el .pka en resultados/propka/")
    a = p.parse_args(argv)
    DIR_RESULTADOS.mkdir(exist_ok=True)

    claves = [a.solo] if a.solo else list(ESTRUCTURAS)
    resumenes, interes = [], []
    for clave in claves:
        r, i, cavs = procesar(clave, a.usar_pka_guardado)
        resumenes.append(r)
        interes.extend(i)
        print(f"{clave:7s} {r['proteina']:28s} {r['n_cavidades']:3d} cavidades | mayor {r['vol_mayor']} Å3")
        for c in cavs[:3]:
            print(f"        cav {c['cavidad']}: {c['volumen']:.0f} Å3, prof. máx {c['prof_max']:.1f} Å, "
                  f"hidropatía {c['hidropatia']:+.2f}, pKa {c['pka_min']}-{c['pka_max']}, "
                  f"activados: {', '.join(c['activados']) or '-'}")

    if not a.solo:
        with open(DIR_RESULTADOS / "cavidades_resumen.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(resumenes[0]))
            w.writeheader(); w.writerows(resumenes)
        with open(DIR_RESULTADOS / "residuos_interes.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(interes[0]))
            w.writeheader(); w.writerows(interes)
        with open(DIR_RESULTADOS / "parametros_cavidades.json", "w") as f:
            import pyKVFinder
            import propka
            json.dump({"detector": "pyKVFinder", "version_pykvfinder": pyKVFinder.__version__,
                       "propka": getattr(propka, "__version__", "3.x"),
                       "parametros": {"step": 0.6, "probe_in": 1.4, "probe_out": 4.0,
                                      "removal_distance": 2.4, "volume_cutoff": 5.0},
                       "umbral_activado": UMBRAL_ACTIVADO}, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
