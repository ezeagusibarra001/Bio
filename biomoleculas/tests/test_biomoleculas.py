import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import estructura_secundaria as es  # noqa: E402
import logo  # noqa: E402
import motivos_uniprot as mu  # noqa: E402


# ---------- Desafío V ----------
def test_prediccion_tiene_un_estado_por_residuo():
    sec = "MLPGLALLLLAAWTMRALEVPTDGNAPLLVEPQIAMFCGR"
    prediccion, detalle = es.predecir(sec)
    assert len(prediccion) == len(sec) == len(detalle)
    assert set(prediccion) <= set("HBL")


def test_formadores_de_helice_hoja_y_loop():
    assert es.predecir("EEEEAAAAMMMMEEEE")[0] == "H" * 16
    assert es.predecir("VVIIVVYYIIVV")[0] == "B" * 12
    assert es.predecir("GPGPGPNGPG")[0] == "L" * 10


def test_suavizado():
    assert es.suavizar("LLHHLLBBBL") == "LLLLLLBBBL"
    assert es.suavizar("HHHHLBL") == "HHHHLLL"


def test_secuencia_invalida():
    with pytest.raises(ValueError):
        es.predecir("MKB1")


# ---------- Desafío VII ----------
def test_motivo_a_regex():
    assert mu.motivo_a_regex("N{P}[ST]{P}") == "N[^P][ST][^P]"
    assert mu.motivo_a_regex("[ST]x[RK]") == "[ST].[RK]"
    with pytest.raises(ValueError):
        mu.motivo_a_regex("N{P}[S1]")


def test_busqueda_con_solapamientos():
    # N2 (NASN), N5 (NNTS), N6 (NTSA); N11 falla por la P siguiente.
    assert mu.buscar_motivo("MNASNNTSAKNPSTNGT") == [2, 5, 6]
    assert mu.buscar_motivo("NPST") == []
    assert mu.buscar_motivo("NSSP") == []   # la 4ta posición no puede ser P


def test_parseo_fasta_uniprot():
    texto = ">sp|P07204|TRBM_HUMAN Thrombomodulin\nMLG\nVL\n"
    assert mu.parsear_fasta(texto) == {"P07204": "MLGVL"}
    assert mu.accession("P07204_TRBM_HUMAN") == "P07204"


# ---------- Desafío VIII ----------
def test_posiciones_variables_del_tp():
    registros = es.leer_fasta((RAIZ / "datos" / "secuencias.fasta").read_text())
    secuencias = [s for _, s in registros]
    logo.validar_alineamiento(secuencias)
    assert logo.posiciones_variables(secuencias) == [15, 27, 30, 103, 129]


def test_matriz_sin_pseudoconteos():
    matriz = logo.matriz_logo(["MAK", "MGK", "MGR"])
    assert matriz.loc[0, "M"] == 1.0
    assert matriz.loc[1, "G"] == pytest.approx(2 / 3)


def test_alineamiento_invalido():
    with pytest.raises(ValueError):
        logo.validar_alineamiento(["MAK", "MA"])
