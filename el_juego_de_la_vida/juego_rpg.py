#!/usr/bin/env python3
"""DESAFÍO IV - "La odisea del gen": un RPG de consola sobre la expresión génica.

Sos la ARN Polimerasa II y tenés que llevar la información de un gen desde el
ADN hasta una proteína funcional. Cada acción consume ATP (tu energía): si te
quedás sin ATP, la célula entra en apoptosis y perdés.

Niveles:
  1. El promotor      -> encontrá la caja TATA para poder unirte al ADN.
  2. Transcripción    -> copiá la hebra molde a ARNm (A->U, T->A, C->G, G->C).
  3. Splicing         -> eliminá el intrón y quedate con los exones.
  4. Exportación      -> sacá el ARNm del núcleo sano y salvo.
  5. Traducción       -> en el ribosoma, traducí los codones a aminoácidos.

Ejemplos:
    python juego_rpg.py
    python juego_rpg.py --nombre Rosalind --dificultad dificil --semilla 7
Comandos durante el juego: 'pista' (cuesta ATP), 'ayuda', 'salir'.
"""

import argparse
import random
import sys
import time

from codigo_genetico import CODON_A_AMINOACIDO, NOMBRES_AMINOACIDOS

CAJA_TATA = "TATAAA"
MOLDE_A_ARN = str.maketrans("ATCG", "UAGC")
DIFICULTADES = {  # ATP inicial, costo de error, largo de los puzzles
    "facil": {"atp": 30, "error": 3, "largo": 6},
    "normal": {"atp": 20, "error": 4, "largo": 9},
    "dificil": {"atp": 14, "error": 5, "largo": 12},
}


class Salir(Exception):
    """Se lanza cuando la persona escribe 'salir'."""


class Colores:
    def __init__(self, activado):
        codigos = {"verde": "92", "rojo": "91", "amarillo": "93", "azul": "94",
                   "magenta": "95", "negrita": "1"}
        for nombre, codigo in codigos.items():
            setattr(self, nombre, (lambda t, c=codigo: f"\033[{c}m{t}\033[0m")
                    if activado else (lambda t: t))


# ---------------------------------------------------------------------------
# Generadores de desafíos (funciones puras, fáciles de testear)
# ---------------------------------------------------------------------------

def adn_al_azar(rng, largo, evitar=CAJA_TATA):
    """ADN aleatorio que no contiene el motivo `evitar`."""
    while True:
        secuencia = "".join(rng.choice("ACGT") for _ in range(largo))
        if evitar not in secuencia:
            return secuencia


def desafio_promotor(rng, largo):
    """ADN con una única caja TATA. Devuelve (secuencia, opciones, posicion_correcta)."""
    total = largo * 3
    posicion = rng.randint(0, total - len(CAJA_TATA))
    relleno = adn_al_azar(rng, total)
    secuencia = relleno[:posicion] + CAJA_TATA + relleno[posicion + len(CAJA_TATA):]
    # Si al pegar la caja se formó otra por accidente, reintentamos.
    if secuencia.count(CAJA_TATA) != 1:
        return desafio_promotor(rng, largo)
    distractores = set()
    while len(distractores) < 3:
        otra = rng.randint(1, total - len(CAJA_TATA) + 1)
        if otra != posicion + 1:
            distractores.add(otra)
    opciones = sorted(distractores | {posicion + 1})
    return secuencia, opciones, posicion + 1


def transcribir(molde):
    """Hebra molde de ADN (3'->5') a ARNm (5'->3')."""
    return molde.translate(MOLDE_A_ARN)


def desafio_splicing(rng, largo):
    """Pre-ARNm con exones en MAYÚSCULA y un intrón en minúscula (GU...AG)."""
    exon1 = transcribir(adn_al_azar(rng, largo // 2 + 2))
    exon2 = transcribir(adn_al_azar(rng, largo // 2 + 2))
    intron = "gu" + "".join(rng.choice("acgu") for _ in range(largo)) + "ag"
    maduro = exon1 + exon2
    incorrectas = {exon1 + intron.upper(), intron.upper() + exon2, exon1 + intron.upper() + exon2}
    opciones = sorted(incorrectas | {maduro}, key=lambda _: rng.random())
    return exon1 + intron + exon2, opciones, maduro


def desafio_traduccion(rng, cantidad):
    """Lista de codones sentido (sin stop) al azar y la proteína correspondiente."""
    sentido = [c for c, aa in CODON_A_AMINOACIDO.items() if aa != "*"]
    codones = ["AUG"] + [rng.choice(sentido) for _ in range(cantidad - 1)]
    proteina = "".join(CODON_A_AMINOACIDO[c] for c in codones)
    return codones, proteina


EVENTOS = [
    ("¡Una nucleasa se acerca a tu ARNm! ¿Qué estructura en el extremo 5' lo protege?",
     ["cola poli-A", "caperuza (cap) de 7-metilguanosina", "caja TATA"], 1),
    ("Un ribosoma perdido pregunta: ¿qué codón indica el INICIO de la traducción?",
     ["UAA", "AUG", "UGA"], 1),
    ("Un virus de ARN te desafía: ¿en qué organela de la célula eucariota ocurre la transcripción?",
     ["núcleo", "mitocondria... ¡ah no! solo sus propios genes", "aparato de Golgi"], 0),
    ("Un factor de transcripción te guiña el ojo: ¿qué proteína se une primero a la caja TATA?",
     ["TBP (proteína de unión a TATA)", "ADN ligasa", "helicasa"], 0),
]


# ---------------------------------------------------------------------------
# Motor del juego
# ---------------------------------------------------------------------------

class Juego:
    def __init__(self, nombre, dificultad, rng, colores, pausa):
        self.nombre = nombre
        self.config = DIFICULTADES[dificultad]
        self.atp = self.config["atp"]
        self.rng = rng
        self.c = colores
        self.pausa = pausa
        self.puntos = 0
        self.proteina = ""

    # --- utilidades de entrada/salida ---
    def decir(self, texto="", demora=0.0):
        print(texto)
        if self.pausa and demora:
            time.sleep(demora)

    def preguntar(self, prompt, pista=None):
        while True:
            try:
                respuesta = input(self.c.azul(f"{prompt} > ")).strip()
            except EOFError:
                raise Salir()
            comando = respuesta.lower()
            if comando == "salir":
                raise Salir()
            if comando == "ayuda":
                self.decir("Comandos: 'pista' (cuesta 2 ATP), 'ayuda', 'salir'. "
                           f"ATP actual: {self.atp}")
                continue
            if comando == "pista":
                if pista:
                    self.gastar(2, "pedir una pista")
                    self.decir(self.c.amarillo(f"💡 {pista}"))
                else:
                    self.decir("No hay pistas para esta pregunta.")
                continue
            return respuesta

    def gastar(self, cantidad, motivo):
        self.atp -= cantidad
        self.decir(self.c.rojo(f"   -{cantidad} ATP por {motivo}. Te quedan {max(self.atp, 0)} ATP."))
        if self.atp <= 0:
            raise Apoptosis()

    def acertar(self, puntos, mensaje="¡Correcto!"):
        self.puntos += puntos
        self.decir(self.c.verde(f"   ✔ {mensaje} (+{puntos} puntos)"), 0.4)

    def elegir_opcion(self, enunciado, opciones, correcta, pista=None, puntos=10):
        """Pregunta de opción múltiple; repite hasta acertar (cada error cuesta ATP)."""
        self.decir(enunciado)
        for i, opcion in enumerate(opciones, 1):
            self.decir(f"   {i}) {opcion}")
        while True:
            respuesta = self.preguntar("Elegí una opción", pista)
            if respuesta.isdigit() and 1 <= int(respuesta) <= len(opciones):
                if opciones[int(respuesta) - 1] == correcta:
                    self.acertar(puntos)
                    return
                self.gastar(self.config["error"], "una respuesta incorrecta")
            else:
                self.decir(f"   Escribí un número entre 1 y {len(opciones)}.")

    def escribir_respuesta(self, enunciado, correcta, pista=None, puntos=15):
        self.decir(enunciado)
        while True:
            respuesta = self.preguntar("Tu respuesta", pista).upper().replace(" ", "")
            if respuesta == correcta:
                self.acertar(puntos)
                return
            self.gastar(self.config["error"], "un error de copia")

    def encabezado(self, numero, titulo):
        self.decir()
        self.decir(self.c.magenta(self.c.negrita(f"═══ NIVEL {numero}: {titulo} ═══")), 0.3)
        self.decir(f"   ⚡ ATP: {self.atp}   ⭐ Puntos: {self.puntos}")

    def evento_aleatorio(self):
        if self.rng.random() < 0.6:
            texto, opciones, correcta = self.rng.choice(EVENTOS)
            self.decir()
            self.decir(self.c.amarillo("⚠ EVENTO ALEATORIO"))
            self.elegir_opcion(texto, opciones, opciones[correcta], puntos=5)

    # --- niveles ---
    def nivel_promotor(self):
        self.encabezado(1, "EL PROMOTOR")
        self.decir("Flotás por el núcleo buscando dónde unirte. Para empezar a transcribir\n"
                   f"necesitás encontrar la caja TATA ({CAJA_TATA}) en la región promotora.", 0.5)
        secuencia, opciones, correcta = desafio_promotor(self.rng, self.config["largo"])
        regla = "".join(str(i % 10) if i % 5 == 0 else "·" for i in range(1, len(secuencia) + 1))
        self.decir(f"\n   5' {secuencia} 3'\n      {regla}")
        self.elegir_opcion("¿En qué posición empieza la caja TATA?",
                           [str(o) for o in opciones], str(correcta),
                           pista=f"La caja empieza con la letra número {correcta}... ¡upa, se me escapó!")

    def nivel_transcripcion(self):
        self.encabezado(2, "TRANSCRIPCIÓN")
        self.decir("¡Te uniste! Ahora tenés que copiar la hebra MOLDE a ARN mensajero.\n"
                   "Recordá: A→U, T→A, C→G, G→C (en el ARN no hay T, hay U).", 0.5)
        molde = adn_al_azar(self.rng, self.config["largo"])
        self.escribir_respuesta(f"\n   Hebra molde  3' {molde} 5'\n   ARNm         5' ??? 3'",
                                transcribir(molde),
                                pista=f"Las primeras 3 bases son {transcribir(molde)[:3]}")

    def nivel_splicing(self):
        self.encabezado(3, "SPLICING")
        self.decir("El pre-ARNm tiene intrones (en minúscula) que no codifican.\n"
                   "El espliceosoma te pide ayuda: ¿cómo queda el ARNm maduro?", 0.5)
        pre, opciones, maduro = desafio_splicing(self.rng, self.config["largo"])
        self.decir(f"\n   pre-ARNm: {pre}")
        self.elegir_opcion("Elegí el ARNm maduro:", opciones, maduro,
                           pista="Los intrones empiezan con GU y terminan con AG. ¡Se van!")

    def nivel_exportacion(self):
        self.encabezado(4, "EXPORTACIÓN")
        self.decir("Tu ARNm maduro tiene que salir del núcleo hacia el citoplasma.", 0.5)
        self.elegir_opcion("¿Por dónde sale?",
                           ["Atraviesa la membrana plasmática",
                            "Por los poros nucleares",
                            "Se queda en el núcleo y los ribosomas van a buscarlo"],
                           "Por los poros nucleares",
                           pista="La envoltura nuclear tiene 'agujeritos' muy regulados.")
        self.evento_aleatorio()

    def nivel_traduccion(self):
        self.encabezado(5, "TRADUCCIÓN")
        cantidad = max(3, self.config["largo"] // 3)
        codones, proteina = desafio_traduccion(self.rng, cantidad)
        self.decir("Llegaste al ribosoma. Traducí cada codón a su aminoácido (código de 1 letra).", 0.5)
        for codon, aa in zip(codones, proteina):
            self.escribir_respuesta(f"\n   Codón {codon} → ?", aa,
                                    pista=f"Es la {NOMBRES_AMINOACIDOS[aa]}", puntos=10)
            self.proteina += aa
            self.decir(f"   Cadena naciente: N-{self.proteina}-C")

    def jugar(self):
        self.decir(self.c.negrita("\n🧬  LA ODISEA DEL GEN  🧬"))
        self.decir(f"Hola, {self.nombre}. Sos la ARN Polimerasa II y tenés {self.atp} ATP.\n"
                   "Escribí 'ayuda' en cualquier momento.", 0.5)
        for nivel in (self.nivel_promotor, self.nivel_transcripcion, self.evento_aleatorio,
                      self.nivel_splicing, self.nivel_exportacion, self.nivel_traduccion):
            nivel()
        self.puntos += self.atp  # el ATP que sobra suma puntos
        self.decir(self.c.verde(self.c.negrita(
            f"\n🎉 ¡Lo lograste! Sintetizaste la proteína N-{self.proteina}-C.")))
        self.decir(f"Puntaje final: {self.puntos} (incluye {self.atp} de ATP restante).")
        return self.puntos


class Apoptosis(Exception):
    """La célula se quedó sin ATP."""


def main(argv=None):
    parser = argparse.ArgumentParser(description="RPG de consola sobre la expresión génica.")
    parser.add_argument("--nombre", default="Polimerasa", help="Tu nombre de jugador/a")
    parser.add_argument("--dificultad", choices=DIFICULTADES, default="normal")
    parser.add_argument("--semilla", type=int, help="Semilla para repetir la misma partida")
    parser.add_argument("--sin-color", action="store_true", help="Desactivar colores ANSI")
    parser.add_argument("--rapido", action="store_true", help="Sin pausas dramáticas")
    args = parser.parse_args(argv)

    colores = Colores(not args.sin_color and sys.stdout.isatty())
    juego = Juego(args.nombre, args.dificultad, random.Random(args.semilla),
                  colores, pausa=not args.rapido)
    try:
        juego.jugar()
    except Apoptosis:
        print(colores.rojo("\n💀 Te quedaste sin ATP. La célula entró en apoptosis. ¡Game over!"))
        print(f"Puntaje: {juego.puntos}")
        return 1
    except (Salir, KeyboardInterrupt):
        print("\n👋 Abandonaste la partida. ¡Hasta la próxima mitosis!")
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
