import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cajas_tata  # noqa: E402
import juego_rpg  # noqa: E402
import proteina_a_arn as p2a  # noqa: E402
from codigo_genetico import CODON_A_AMINOACIDO, traducir  # noqa: E402

SEC1 = "ATVEKGGKHKTGPNEKGKKIFVQKCSQCHTVLHGLFGRKTGQA"


# ---------- Código genético ----------
def test_codigo_genetico_completo():
    assert len(CODON_A_AMINOACIDO) == 64
    assert sum(aa == "*" for aa in CODON_A_AMINOACIDO.values()) == 3
    assert len(set(CODON_A_AMINOACIDO.values()) - {"*"}) == 20


# ---------- Desafío II ----------
@pytest.mark.parametrize("modo", ["primero", "aleatorio"])
def test_arn_codifica_la_proteina(modo):
    arn = p2a.proteina_a_arn(SEC1, modo=modo, rng=random.Random(0))
    assert len(arn) == 3 * len(SEC1)
    assert set(arn) <= set("ACGU")
    assert traducir(arn) == SEC1


def test_inicio_y_stop():
    arn = p2a.proteina_a_arn("ATV", inicio=True, stop=True)
    assert arn.startswith("AUG")
    assert CODON_A_AMINOACIDO[arn[-3:]] == "*"
    assert traducir(arn) == "MATV"


def test_secuencia_invalida():
    with pytest.raises(ValueError):
        p2a.proteina_a_arn("AB1")


def test_cantidad_de_opciones():
    assert p2a.cantidad_de_arns_posibles("MW") == 1
    assert p2a.cantidad_de_arns_posibles("L") == 6
    assert len(list(p2a.todas_las_opciones("LS"))) == 36


# ---------- Desafío III ----------
def test_regiones_consecutivas_y_pares():
    sec = "GGTATAAACCCTATAAAGGGTATAAAC"
    assert cajas_tata.posiciones_motivo(sec) == [2, 11, 20]
    consecutivas = cajas_tata.regiones_promotoras(sec)
    assert [(r["inicio"], r["fin"]) for r in consecutivas] == [(3, 17), (12, 26)]
    assert all(r["secuencia"].startswith("TATAAA") and r["secuencia"].endswith("TATAAA")
               for r in consecutivas)
    assert len(cajas_tata.regiones_promotoras(sec, modo="pares")) == 1


def test_sin_cajas_suficientes():
    assert cajas_tata.regiones_promotoras("ACGTTATAAAGC") == []


def test_reverso_complementario():
    assert cajas_tata.reverso_complementario("TTTATA") == "TATAAA"


def test_lectura_fasta(tmp_path):
    archivo = tmp_path / "x.fasta"
    archivo.write_text(">a desc\nacgt\nTATA\n>b\nGG\n")
    assert cajas_tata.leer_secuencias(archivo) == [("a", "ACGTTATA"), ("b", "GG")]


# ---------- Desafío IV ----------
def test_desafios_del_juego_son_consistentes():
    rng = random.Random(1)
    for _ in range(50):
        sec, opciones, correcta = juego_rpg.desafio_promotor(rng, 9)
        assert sec.count("TATAAA") == 1 and sec.find("TATAAA") + 1 == correcta
        assert correcta in opciones and len(opciones) == 4
        pre, opciones, maduro = juego_rpg.desafio_splicing(rng, 6)
        assert "".join(c for c in pre if c.isupper()) == maduro and maduro in opciones
    assert juego_rpg.transcribir("TACGGA") == "AUGCCU"


class _Bot(juego_rpg.Juego):
    """Jugador automático: responde bien (o siempre mal si `perdedor`)."""

    def __init__(self, perdedor=False, **kw):
        super().__init__("Bot", "normal", random.Random(5), juego_rpg.Colores(False), False, **kw)
        self.perdedor, self._resp = perdedor, []

    def elegir_opcion(self, enunciado, opciones, correcta, pista=None, puntos=10):
        indice = opciones.index(correcta)
        if self.perdedor:
            indice = (indice + 1) % len(opciones)
        self._resp = [str(indice + 1)] * 100
        super().elegir_opcion(enunciado, opciones, correcta, pista, puntos)

    def escribir_respuesta(self, enunciado, correcta, pista=None, puntos=15):
        self._resp = ["???" if self.perdedor else correcta] * 100
        super().escribir_respuesta(enunciado, correcta, pista, puntos)

    def preguntar(self, prompt, pista=None):
        return self._resp.pop(0)


def test_bot_gana_la_partida(capsys):
    bot = _Bot()
    puntos = bot.jugar()
    assert bot.proteina.startswith("M") and puntos > 0
    assert "Lo lograste" in capsys.readouterr().out


def test_bot_pierde_por_apoptosis():
    with pytest.raises(juego_rpg.Apoptosis):
        _Bot(perdedor=True).jugar()
