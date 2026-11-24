import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import analisis_citocromo_c as cc  # noqa: E402
import descargar_datos as dd  # noqa: E402

FASTA = RAIZ / "datos" / "citocromo_c.fasta"
ALN = RAIZ / "datos" / "citocromo_c_mafft.aln"


# ---------- Lectura de datos ----------
def test_fasta_tiene_las_nueve_secuencias_de_la_tabla():
    seqs = cc.leer_fasta(FASTA)
    assert list(seqs) == dd.ACCESIONES
    assert set(seqs) == set(cc.ESPECIES)


def test_fasta_humano_coincide_con_el_enunciado():
    humano = cc.leer_fasta(FASTA)["NP_061820.1"]
    assert humano.startswith("MGDVEKGKKIFIMKCSQCHTVEKGGKHK")
    assert len(humano) == 105  # 104 de la proteína madura + la Met inicial


def test_alineamiento_conserva_las_secuencias_originales():
    seqs = cc.leer_fasta(FASTA)
    aln, _ = cc.leer_clustal(ALN)
    assert set(aln) == set(seqs)
    largos = {len(s) for s in aln.values()}
    assert len(largos) == 1                     # todas del mismo largo
    for acc, s in aln.items():
        assert s.replace("-", "") == seqs[acc]  # alinear no pierde residuos


# ---------- Línea de conservación ----------
@pytest.mark.parametrize("columna, simbolo", [
    ("KKKK", "*"),     # todos iguales
    ("IIVV", ":"),     # I y V comparten el grupo fuerte MILV
    ("SSAA", ":"),     # grupo fuerte STA
    ("SSGG", "."),     # S y G solo comparten grupos débiles (SAG, SGND)
    ("KKDW", " "),     # nada en común
    ("KK-K", " "),     # una columna con gap no se marca
])
def test_simbolo_conservacion(columna, simbolo):
    assert cc.simbolo_conservacion(list(columna)) == simbolo


def test_linea_recalculada_igual_a_la_del_programa():
    aln, cons_programa = cc.leer_clustal(ALN)
    assert cc.linea_conservacion(aln) == cons_programa


# ---------- Identidad y diferencias ----------
def test_identidad_ignora_gaps():
    assert cc.identidad("A-CD", "AACD") == 100.0
    assert cc.identidad("ACDE", "ACDF") == 75.0


def test_humano_chimpance_y_gorila_son_identicos():
    r = cc.analizar(ALN)
    a = r["alineamiento"]
    assert a["Humano"] == a["Chimpancé"] == a["Gorila"]


def test_humano_vs_gallina_13_diferencias():
    r = cc.analizar(ALN)
    difs = r["humano_vs_gallina"]
    assert len(difs) == 13
    assert (13, "M", "Q") in difs
    i = r["nombres"].index("Humano")
    j = r["nombres"].index("Gallina")
    assert r["identidad"][i][j] == pytest.approx(100 * 92 / 105)


def test_la_mosca_es_la_mas_distinta_de_todas():
    r = cc.analizar(ALN)
    nombres, mid = r["nombres"], r["identidad"]
    k = nombres.index("Mosca de la fruta")
    for i, n in enumerate(nombres):
        if n == "Mosca de la fruta":
            continue
        otros = [mid[i][j] for j in range(len(nombres)) if j not in (i, k)]
        assert mid[i][k] < min(otros)


# ---------- Árboles ----------
D_JUGUETE = [  # árbol aditivo conocido: ((A:1,B:2):1,(C:3,D:1)) sin raíz
    [0, 3, 5, 3],
    [3, 0, 6, 4],
    [5, 6, 0, 4],
    [3, 4, 4, 0],
]


def test_neighbor_joining_recupera_un_arbol_aditivo():
    arbol = cc.neighbor_joining(list("ABCD"), D_JUGUETE)
    # NJ debe juntar A con B y C con D; la suma de ramas del árbol real es 8
    def suma(n):
        return sum(l + suma(h) for h, l in n.hijos)
    assert suma(arbol) == pytest.approx(8)
    newick = arbol.newick(1)
    assert "(A:1.0,B:2.0)" in newick


def test_upgma_es_ultrametrico():
    arbol = cc.upgma(list("ABCD"), D_JUGUETE)

    def alturas(n, acum=0.0):
        if n.es_hoja():
            return [acum]
        return [a for h, l in n.hijos for a in alturas(h, acum + l)]
    hs = alturas(arbol)
    assert max(hs) - min(hs) == pytest.approx(0)


def test_enraizar_deja_al_grupo_externo_solo_de_un_lado():
    r = cc.analizar(ALN)
    raiz = r["nj"]
    assert len(raiz.hijos) == 2
    lado_a, lado_b = (h.hojas() for h, _ in raiz.hijos)
    assert lado_a == ["Mosca de la fruta"]
    assert len(lado_b) == 8


def test_primates_forman_un_grupo_en_ambos_arboles():
    r = cc.analizar(ALN)
    for arbol in (r["upgma"], r["nj"]):
        assert _hay_clado(arbol, {"Humano", "Chimpancé", "Gorila"})
        assert _hay_clado(arbol, {"Caballo", "Perro"})


def _hay_clado(nodo, conjunto):
    if set(nodo.hojas()) == conjunto:
        return True
    return any(_hay_clado(h, conjunto) for h, _ in nodo.hijos)


# ---------- Descarga (sin red) ----------
def test_url_efetch_pide_fasta_de_proteinas():
    url = dd.url_efetch(["NP_061820.1", "NP_477164.1"])
    assert "db=protein" in url and "rettype=fasta" in url
    assert "NP_061820.1%2CNP_477164.1" in url
