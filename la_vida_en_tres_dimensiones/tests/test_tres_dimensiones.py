import json
import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import cavidades as cv  # noqa: E402
import explorar_pdb as ep  # noqa: E402
import numeracion as nu  # noqa: E402
import plddt as pl  # noqa: E402
import sitios_uniprot as su  # noqa: E402

DATOS = RAIZ / "datos"
RES = RAIZ / "resultados"
UBQ = DATOS / "1UBQ.pdb"
EST = DATOS / "estructuras"


# ---------- Parser PDB (ubiquitina) ----------
def test_tres_a_una():
    assert ep.tres_a_una(["MET", "GLN", "ILE", "PHE"]) == "MQIF"
    assert ep.tres_a_una(["XYZ"]) == "X"


def test_encabezado_1ubq():
    info = ep.leer_encabezado(UBQ)
    assert info["id"] == "1UBQ"
    assert info["metodo"] == "X-RAY DIFFRACTION"
    assert info["resolucion"] == pytest.approx(1.8)
    assert info["dbref"][0]["acceso"] == "P62988"


def test_secuencia_seqres_igual_a_ca():
    seqres = ep.secuencia_seqres(UBQ, "A")
    ca = ep.secuencia_ca(UBQ, "A")
    assert len(seqres) == 76
    assert seqres == ca
    assert seqres.startswith("MQIFVKTLTG")
    assert seqres.endswith("LRLRGG")


def test_conteo_de_registros():
    c = ep.contar_registros(UBQ)
    assert c["ATOM"] == 602
    assert c["HETATM"] == c["aguas"] == 58      # todos los HETATM son aguas


def test_estructura_secundaria_y_factor_b():
    ss = ep.estructura_secundaria(UBQ)
    assert ("hélice α", "A", 23, 34) in ss
    assert sum(1 for t, *_ in ss if t == "hebra β") == 5
    b = ep.bfactor_por_residuo(UBQ, "A")
    assert max(b, key=b.get) == 76             # la cola C-terminal es lo más móvil


# ---------- Numeración PDB vs UniProt ----------
def test_offsets_dbref():
    assert ep.offset_uniprot(EST / "1M14.pdb", "A") == 24
    assert ep.offset_uniprot(EST / "3W32.pdb", "A") == 0
    assert ep.offset_uniprot(EST / "1AO6.pdb", "A") == 24


def test_asp837_de_uniprot_en_cada_pdb():
    egfr = su.leer_entrada("P00533")
    assert egfr["sequence"]["value"][836] == "D"                       # D837 en UniProt
    assert su.residuo_en_pdb(EST / "1M14.pdb", "A", 837) == (813, "ASP")
    assert su.residuo_en_pdb(EST / "3W32.pdb", "A", 837) == (837, "ASP")
    # la trampa: el "837" del archivo 1M14 es otro residuo
    assert ep.residuos_ca(EST / "1M14.pdb", "A")[837] == "LEU"


def test_alineamiento_identico():
    assert nu.alinear_global("ACDEF", "ACDEF") == ("ACDEF", "ACDEF")
    a1, a2 = nu.alinear_global("ACDEF", "ACEF")
    assert a2.count("-") == 1 and a2.replace("-", "") == "ACEF"


def test_equivalencias_hsa_bsa():
    filas, ident, mapa = nu.equivalencias_albuminas()
    tabla = {r_h: (aa_h, r_b, aa_b) for r_h, aa_h, r_b, aa_b in filas}
    assert tabla[199] == ("K", 198, "R")     # la Lys199 humana es Arg en la bovina
    assert tabla[222] == ("R", 221, "K")     # ...y la Lys221 bovina ocupa el lugar de Arg222
    assert 70 < ident < 82


# ---------- UniProt ----------
def test_sitios_akr1c4():
    e = su.leer_entrada("P17516")
    s = su.sitios(e, ("Active site", "Site"))
    assert {(x["inicio"], x["residuos"]) for x in s} == {(55, "Y"), (84, "K")}
    nadp = [x for x in su.sitios(e, ("Binding site",)) if x["inicio"] == 20]
    assert nadp[0]["residuos"] == "GFGTY" and nadp[0]["ligando"] == "NADP(+)"
    assert su.residuo_en_pdb(EST / "2FVL.pdb", "B", 55) == (55, "TYR")


# ---------- PROPKA y cavidades (sobre archivos guardados) ----------
def test_leer_propka_hsa():
    pk = cv.leer_propka(RES / "propka" / "1AO6_A.pka")
    assert pk[199]["residuo"] == "LYS"
    assert pk[199]["pka"] == pytest.approx(6.78)
    assert pk[199]["shift"] == pytest.approx(6.78 - 10.5, abs=0.01)
    assert pk[199]["enterrado"] == 100
    assert all(d["pka"] < 99 for d in pk.values())      # sin Cys de puentes disulfuro


def test_anotar_cavidades_net_shift():
    cavs = [{"residuos": [(1, "LYS"), (2, "ALA"), (3, "ASP")]}]
    pkas = {1: {"residuo": "LYS", "pka": 8.0, "shift": -2.5},
            3: {"residuo": "ASP", "pka": 4.3, "shift": 0.5}}
    c = cv.anotar_cavidades(cavs, pkas, ligando={3})[0]
    assert c["n_ionizables"] == 2
    assert c["net_shift"] == pytest.approx(3.0)          # |−2.5| + |0.5|
    assert c["activados"] == ["LYS1(8.00)"]               # solo |shift| >= 1
    assert c["contacto_ligando"] == 1


def test_extraer_cadena_sin_ligandos(tmp_path):
    destino = tmp_path / "2FVL_B.pdb"
    n = cv.extraer_cadena(EST / "2FVL.pdb", "B", destino)
    lineas = [l for l in destino.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert n == len(lineas) > 2000
    assert all(l.startswith("ATOM") and l[21] == "B" for l in lineas)


def test_ligando_nadp_toca_sitios_de_uniprot():
    cerca = cv.residuos_cerca_de_ligando(EST / "2FVL.pdb", "B", ["NAP"])
    assert {50, 55, 84, 117}.issubset(cerca)


def test_resultados_guardados_coherentes():
    hsa = cv.leer_cavidades_csv(RES / "cavidades_1AO6_A.csv")
    assert hsa[0]["cavidad"] == 1 and "K199" in hsa[0]["residuos"]
    assert [c["volumen"] for c in hsa] == sorted((c["volumen"] for c in hsa), reverse=True)
    bsa = cv.leer_cavidades_csv(RES / "cavidades_4JK4_A.csv")
    assert any("K221" in c["residuos"] and "LYS221(5.39)" in c["activados"] for c in bsa)


# ---------- pLDDT ----------
def test_plddt_pdb_igual_a_json_de_confianza():
    pdb = pl.leer_plddt(DATOS / "alphafold" / "AF-0000000365760270-model_v1.pdb")
    js = pl.leer_plddt(DATOS / "alphafold" / "AF-0000000365760270-confidence_v1.json")
    assert len(pdb) == len(js) == 98
    assert pdb == pytest.approx(js)


def test_plddt_formato_colabfold(tmp_path):
    f = tmp_path / "x_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"
    f.write_text(json.dumps({"plddt": [95.0, 80.0, 40.0], "pae": [[0, 1, 2]] * 3, "ptm": 0.5}))
    assert pl.leer_plddt(f) == [95.0, 80.0, 40.0]
    assert pl.leer_pae(f)[0] == [0, 1, 2]


def test_resumen_y_segmentos_e7():
    v = pl.leer_plddt(DATOS / "alphafold" / "AF-0000000365760270-confidence_v1.json")
    r = pl.resumen(v)
    assert r["media"] == pytest.approx(81.34, abs=0.01)
    assert sum(r["fracciones"].values()) == pytest.approx(1.0, abs=0.02)
    assert ("alta", 45, 98) in pl.segmentos(v)            # dominio C-terminal de dedo de zinc
    assert min(v[27:44]) < 50                             # tramo ácido desordenado


def test_categorias():
    assert pl.categoria(95) == "muy alta (>90)"
    assert pl.categoria(75) == "confiable (70-90)"
    assert pl.categoria(30) == "muy baja (<50)"


def test_grafico_se_genera(tmp_path):
    salida = tmp_path / "plddt.png"
    pl.graficar([90, 60, 30, 80], salida, pae=[[0, 5, 9, 2]] * 4)
    assert salida.exists() and salida.stat().st_size > 1000
