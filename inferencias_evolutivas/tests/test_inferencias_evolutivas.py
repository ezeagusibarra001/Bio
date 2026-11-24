import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import arboles  # noqa: E402
import comparar  # noqa: E402
import secuencias  # noqa: E402

BART = ("MGSGDAENGKKIFVQKCAQCHTYEVGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLE"
        "NPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE")


# ---------- Secuencias ----------
def test_etiqueta_corta_saca_el_organismo():
    assert secuencias.etiqueta_corta("NP_061820.1_cytochrome_c_[Homo_sapiens]") == "Homo_sapiens"
    assert secuencias.etiqueta_corta(
        "NP_001285984.1_cytochrome_c_proximal,_isoform_B_[Drosophila_melanogaster]") == "Drosophila_melanogaster"
    assert secuencias.etiqueta_corta("bartmosca") == "bartmosca"


def test_archivo_del_curso_tiene_9_secuencias_y_bartmosca():
    registros = dict(secuencias.leer_fasta(secuencias.ORIGINAL))
    assert len(registros) == 9
    assert registros["bartmosca"] == BART


def test_bartmosca_comparte_los_primeros_30_residuos_con_la_mosca_distal():
    (_, distal), = secuencias.leer_fasta(secuencias.DISTAL)
    assert BART[:30] == distal[:30]
    assert BART[30] != distal[30]


def test_bartmosca_es_humano_desde_la_posicion_26():
    humano = dict(secuencias.leer_fasta(secuencias.ORIGINAL))["NP_061820.1_cytochrome_c_[Homo_sapiens]"]
    # El humano tiene 2 residuos menos al principio: bart[25:] == humano[23:]
    assert BART[25:] == humano[23:]
    assert len(BART[25:]) == 82


# ---------- Identidad / distancia ----------
def test_identidad_ignora_gaps():
    assert comparar.identidad("AC-GT", "ACTGA") == pytest.approx(3 / 4)
    assert comparar.p_distancia("ACGT", "ACGT") == 0


def test_identidad_exige_mismo_largo():
    with pytest.raises(ValueError):
        comparar.identidad("ACG", "AC")


def test_clasificar_columnas():
    filas = comparar.clasificar_columnas("AXCD-", "AYCZE", "AXWZE")
    assert [f["clase"] for f in filas] == ["ambos", "solo_mosca", "solo_humano", "ninguno"]
    assert [f["posicion_bart"] for f in filas] == [1, 2, 3, 4]


def test_punto_de_corte_quimera_perfecta():
    izquierda = "AAAAAAAAAA"
    derecha = "BBBBBBBBBB"
    quimera = "AAAABBBBBB"
    mejor, total, cortes = comparar.punto_de_corte(quimera, izquierda, derecha)
    assert (mejor, total, cortes) == (10, 10, [4])


def test_punto_de_corte_ambiguo_si_las_fuentes_coinciden():
    # En las posiciones 3-4 las dos fuentes son iguales: el corte puede estar en 2, 3 o 4.
    mejor, total, cortes = comparar.punto_de_corte("AACCBB", "AACCAA", "BBCCBB")
    assert mejor == total == 6
    assert cortes == [2, 3, 4]


def test_alineamiento_guardado_set4_explica_bartmosca_entera():
    aln = dict(secuencias.leer_fasta(RAIZ / "resultados" / "aln_mafft_set4.fasta"))
    mejor, total, _ = comparar.punto_de_corte(aln["bartmosca"], aln["Drosophila_distal"], aln["Homo_sapiens"])
    assert mejor == total == len(BART)


# ---------- Árboles ----------
NEWICK = "(A:0.1,(B:0.2,(C:0.3,D:0.4)80/95:0.05)70/60:0.06,E:0.5);"


def test_leer_newick():
    raiz = arboles.leer_newick(NEWICK)
    assert sorted(raiz.hojas()) == ["A", "B", "C", "D", "E"]
    assert arboles.soporte_de(raiz, ["C", "D"]) == "80/95"
    assert arboles.soporte_de(raiz, ["B", "C", "D"]) == "70/60"


def _largo_total(nodo):
    return sum(n.largo for n in arboles._todos(nodo) if n.padre is not None)


def test_enraizar_conserva_hojas_largos_y_soportes():
    raiz = arboles.leer_newick(NEWICK)
    enraizado = arboles.enraizar(raiz, "E")
    assert sorted(enraizado.hojas()) == ["A", "B", "C", "D", "E"]
    assert _largo_total(enraizado) == pytest.approx(_largo_total(raiz))
    # La raíz separa al outgroup del resto.
    assert sorted(len(h.hojas()) for h in enraizado.hijos) == [1, 4]
    # El soporte de (C,D) viaja con su rama.
    assert arboles.soporte_de(enraizado, ["C", "D"]) == "80/95"


def test_enraizar_mueve_soporte_al_clado_complementario():
    # Sin raíz, '70/60' es el split {B,C,D} | {A,E}. Enraizando en B, ese split
    # aparece como el clado {A,E} y debe conservar el soporte.
    enraizado = arboles.enraizar(arboles.leer_newick(NEWICK), "B")
    assert arboles.soporte_de(enraizado, ["A", "E"]) == "70/60"


def test_arbol_guardado_ubica_a_bartmosca_con_el_humano():
    ruta = RAIZ / "resultados" / "iqtree" / "set9_mfp.treefile"
    arbol = arboles.enraizar(arboles.leer_newick(ruta.read_text()), arboles.OUTGROUP)
    assert arboles.hermanas(arbol, "bartmosca") == ["Homo_sapiens"]
