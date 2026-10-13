import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import needleman_wunsch as nw  # noqa: E402


# ---------- Matriz de programación dinámica ----------
def test_dimensiones_de_la_matriz():
    F = nw.alinear("ANA", "ANANA")["matriz"]
    assert len(F) == 4            # len(sec1) + 1
    assert all(len(fila) == 6 for fila in F)   # len(sec2) + 1


def test_inicializacion_con_gaps():
    F = nw.alinear("ANA", "ANANA", gap=-2)["matriz"]
    assert F[0] == [0, -2, -4, -6, -8, -10]
    assert [fila[0] for fila in F] == [0, -2, -4, -6]


# ---------- Score y alineamiento ----------
def test_secuencias_identicas():
    r = nw.alinear("GATTACA", "GATTACA")
    assert r["alineamiento"] == ("GATTACA", "GATTACA")
    assert r["identidad"] == 100.0
    assert r["score"] == 7        # 7 matches de +1


def test_banana_manzana():
    r = nw.alinear("BANANA", "MANZANA")
    a1, a2 = r["alineamiento"]
    assert len(a1) == len(a2)                 # mismo largo tras alinear
    assert a1.replace("-", "") == "BANANA"    # no se pierden residuos
    assert a2.replace("-", "") == "MANZANA"
    assert r["score"] == 2


def test_un_gap_necesario():
    # ANA dentro de ANANA: hacen falta 2 gaps y las 3 columnas sin gap coinciden.
    a1, a2 = nw.alinear("ANA", "ANANA")["alineamiento"]
    assert a1.count("-") == 2
    assert a1.replace("-", "") == "ANA"


# ---------- Efecto del gap penalty (Desafío VI) ----------
def test_gap_mas_barato_sube_el_score():
    suave = nw.alinear("ANA", "ANANA", gap=-1)["score"]
    duro = nw.alinear("ANA", "ANANA", gap=-4)["score"]
    assert suave > duro


# ---------- Cálculo de identidad ----------
def test_identidad_ignora_columnas_con_gap():
    # 3 columnas sin gap, todas iguales -> 100%, aunque haya gaps.
    assert nw.identidad("A-NA", "AANA") == 100.0


def test_identidad_cuenta_coincidencias():
    assert nw.identidad("GATTACA", "GCATGCU") == pytest.approx(42.9, abs=0.1)
