#!/usr/bin/env python3
"""DESAFÍO III - Árboles filogenéticos con IQ-TREE, enraizado y gráfico.

1. Corre IQ-TREE sobre los alineamientos del Desafío II con distintas
   configuraciones, para poder discutir modelo de sustitución y bootstrap:

     set9_mfp      ModelFinder elige el modelo (-m MFP) + UFBoot 1000 + SH-aLRT 1000
     set9_poisson  modelo más simple posible (Poisson: todas las sustituciones iguales) + UFBoot 1000
     set9_b10      modelo elegido + bootstrap no paramétrico estándar con 10 réplicas
     set9_b100     ídem con 100 réplicas
     set9_b1000    ídem con 1000 réplicas
     set3_mfp      el set de 3 (humano, mosca, bartmosca)

2. Lee los .treefile (Newick), **enraíza** el árbol con un grupo externo
   (Rhizobium, una bacteria: ver README) y lo dibuja con matplotlib, con los
   soportes en cada nodo. FigTree hace lo mismo de forma interactiva.

Requiere IQ-TREE (`brew install iqtree3`, o iqtree2). Las corridas usan una
semilla fija (-seed) para que sean reproducibles.

Uso:
    python arboles.py            # corre IQ-TREE (si no hay resultados) y grafica
    python arboles.py --solo-graficar
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
RESULTADOS = AQUI / "resultados"
IQTREE_DIR = RESULTADOS / "iqtree"
SEMILLA = "20261124"
OUTGROUP = "Rhizobium_sp__CF097"


# ---------- Newick ----------
class Nodo:
    def __init__(self, nombre="", largo=0.0, soporte=""):
        self.nombre = nombre
        self.largo = largo
        self.soporte = soporte
        self.hijos = []
        self.padre = None

    def agregar(self, hijo):
        hijo.padre = self
        self.hijos.append(hijo)
        return hijo

    def es_hoja(self):
        return not self.hijos

    def hojas(self):
        if self.es_hoja():
            return [self.nombre]
        return [h for hijo in self.hijos for h in hijo.hojas()]

    def to_newick(self):
        if self.es_hoja():
            texto = self.nombre
        else:
            texto = "(" + ",".join(h.to_newick() for h in self.hijos) + ")" + self.soporte
        return f"{texto}:{self.largo:.6g}" if self.padre is not None else texto + ";"


def leer_newick(texto):
    """Parser Newick mínimo: nombres, largos de rama y etiquetas de soporte
    en nodos internos (IQ-TREE escribe 'SH-aLRT/UFBoot' o solo 'bootstrap')."""
    tokens = re.findall(r"[(),:;]|[^(),:;]+", texto.strip())
    raiz = Nodo()
    actual = raiz
    esperando_largo = False
    ultimo = raiz
    for tok in tokens:
        if tok == "(":
            actual = actual.agregar(Nodo())
            ultimo = actual
        elif tok == ",":
            actual = actual.padre.agregar(Nodo())
            ultimo = actual
        elif tok == ")":
            actual = actual.padre
            ultimo = actual
        elif tok == ":":
            esperando_largo = True
        elif tok == ";":
            break
        else:
            tok = tok.strip()
            if esperando_largo:
                ultimo.largo = float(tok)
                esperando_largo = False
            elif ultimo.es_hoja() and not ultimo.nombre:
                ultimo.nombre = tok
            else:
                ultimo.soporte = tok
    return raiz


def enraizar(raiz, outgroup):
    """Re-enraíza el árbol en el punto medio de la rama que lleva al outgroup.

    Los árboles de IQ-TREE salen **sin raíz** (dibujados desde un nodo
    cualquiera). Elegir un grupo externo fija la dirección del tiempo:
    el ancestro común de todo lo demás queda del otro lado de esa rama.
    """
    # 1. Pasar el árbol a un grafo no dirigido. Cada arista guarda su largo y
    #    su soporte: en el Newick de IQ-TREE el soporte escrito en un nodo
    #    interno es el de la rama que lo une con su padre (es un dato de la
    #    *rama*, no del nodo), así que al re-enraizar viaja con la arista.
    vecinos = {}
    nodos = {}
    for n in _todos(raiz):
        nodos[id(n)] = n
        vecinos.setdefault(id(n), [])
        if n.padre is not None:
            vecinos[id(n)].append((id(n.padre), n.largo, n.soporte))
            vecinos.setdefault(id(n.padre), []).append((id(n), n.largo, n.soporte))

    hoja = next(n for n in nodos.values() if n.nombre == outgroup)
    (id_vecino, largo, _), = vecinos[id(hoja)]

    # 2. Reconstruir un árbol dirigido desde el vecino del outgroup, sin volver atrás.
    def copiar(id_nodo, desde, largo_rama, soporte):
        original = nodos[id_nodo]
        nuevo = Nodo(original.nombre, largo_rama, soporte if original.nombre == "" else "")
        for id_hijo, l, s in vecinos[id_nodo]:
            if id_hijo != desde:
                nuevo.agregar(copiar(id_hijo, id_nodo, l, s))
        # Un nodo interno que quedó con un solo hijo (la vieja raíz binaria) se colapsa.
        if len(nuevo.hijos) == 1 and not nuevo.nombre:
            unico = nuevo.hijos[0]
            unico.largo += nuevo.largo
            return unico
        return nuevo

    # 3. Raíz nueva en el punto medio de la rama del outgroup.
    nueva = Nodo()
    nueva.agregar(Nodo(outgroup, largo / 2))
    nueva.agregar(copiar(id_vecino, id(hoja), largo / 2, ""))
    return nueva


def _todos(nodo):
    yield nodo
    for h in nodo.hijos:
        yield from _todos(h)


def clados(raiz):
    """Conjunto de clados (frozenset de hojas) de los nodos internos no raíz."""
    return {frozenset(n.hojas()) for n in _todos(raiz) if n.hijos and n.padre is not None}


def soporte_de(raiz, hojas):
    """Etiqueta de soporte del nodo cuyo clado es exactamente `hojas`."""
    objetivo = frozenset(hojas)
    for n in _todos(raiz):
        if n.hijos and frozenset(n.hojas()) == objetivo:
            return n.soporte
    return None


def hermanas(raiz, hoja):
    """Hojas del grupo hermano de `hoja` (lo que comparte el nodo padre con ella)."""
    n = next(x for x in _todos(raiz) if x.nombre == hoja)
    return sorted(set(n.padre.hojas()) - {hoja})


# ---------- Gráfico ----------
def dibujar(raiz, ruta, titulo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y = {}
    hojas = []

    def asignar_y(n):
        if n.es_hoja():
            y[id(n)] = len(hojas)
            hojas.append(n)
        else:
            for h in n.hijos:
                asignar_y(h)
            y[id(n)] = sum(y[id(h)] for h in n.hijos) / len(n.hijos)

    asignar_y(raiz)
    x = {}

    def asignar_x(n, base):
        x[id(n)] = base + (n.largo if n.padre is not None else 0)
        for h in n.hijos:
            asignar_x(h, x[id(n)])

    asignar_x(raiz, 0.0)
    ancho = max(x.values())

    fig, ax = plt.subplots(figsize=(9, 0.5 * len(hojas) + 1.5))
    for n in _todos(raiz):
        if n.padre is not None:
            ax.plot([x[id(n.padre)], x[id(n)]], [y[id(n)], y[id(n)]], color="black", lw=1.4)
        if n.hijos:
            ys = [y[id(h)] for h in n.hijos]
            ax.plot([x[id(n)], x[id(n)]], [min(ys), max(ys)], color="black", lw=1.4)
            if n.soporte and n.padre is not None:
                ax.text(x[id(n)] - 0.005 * ancho, y[id(n)] - 0.15, n.soporte, ha="right", va="bottom",
                        fontsize=7, color="#555555")
    for h in hojas:
        color = "#d62728" if h.nombre == "bartmosca" else "black"
        peso = "bold" if h.nombre == "bartmosca" else "normal"
        ax.text(x[id(h)] + 0.01 * ancho, y[id(h)], h.nombre.replace("__", ". ").replace("_", " "), va="center",
                fontsize=9, color=color, fontweight=peso, style="italic" if h.nombre != "bartmosca" else "normal")
    ax.plot([0, 0.1], [-0.9, -0.9], color="black", lw=1.2)
    ax.text(0.05, -1.05, "0.1 sust./sitio", ha="center", va="bottom", fontsize=7)
    ax.set_xlim(-0.02 * ancho, ancho * 1.35)
    ax.set_ylim(-1.6, len(hojas) - 0.4)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(titulo, fontsize=10)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


# ---------- IQ-TREE ----------
def binario_iqtree():
    for nombre in ("iqtree3", "iqtree2", "iqtree"):
        if shutil.which(nombre):
            return nombre
    raise SystemExit("No se encontró IQ-TREE (brew install iqtree3)")


def alineamiento(nombre_set):
    """Prefiere el alineamiento de Clustal Omega (Desafío II); si no está, MAFFT."""
    for ruta in (RESULTADOS / f"aln_{nombre_set}.fasta", RESULTADOS / f"aln_mafft_{nombre_set}.fasta"):
        if ruta.exists():
            return ruta
    raise SystemExit(f"Falta el alineamiento de {nombre_set}")


def modelo_elegido(prefijo):
    texto = (IQTREE_DIR / f"{prefijo}.iqtree").read_text()
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", texto)
    return m.group(1) if m else None


def correr(prefijo, aln, *opciones):
    IQTREE_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [binario_iqtree(), "-s", str(aln), "--prefix", str(IQTREE_DIR / prefijo),
           "-seed", SEMILLA, "-T", "1", "-redo", "--quiet", *opciones]
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def correr_todo():
    aln9, aln3 = alineamiento("set9"), alineamiento("set3")
    correr("set9_mfp", aln9, "-m", "MFP", "-B", "1000", "-alrt", "1000")
    modelo = modelo_elegido("set9_mfp")
    print(f"Modelo elegido por ModelFinder (BIC): {modelo}")
    correr("set9_poisson", aln9, "-m", "Poisson", "-B", "1000", "-alrt", "1000")
    for n in (10, 100, 1000):
        correr(f"set9_b{n}", aln9, "-m", modelo, "-b", str(n))
    correr("set3_mfp", aln3, "-m", "MFP")
    arboles_por_region(modelo)


def subalineamiento(ruta_aln, referencia, desde, hasta, destino):
    """Guarda las columnas donde `referencia` tiene sus residuos desde..hasta
    (posiciones 1-based, sin contar gaps). Sirve para armar un árbol por región."""
    from secuencias import escribir_fasta, leer_fasta

    aln = leer_fasta(ruta_aln)
    ref = dict(aln)[referencia]
    columnas, pos = [], 0
    for c, residuo in enumerate(ref):
        if residuo != "-":
            pos += 1
            if desde <= pos <= hasta:
                columnas.append(c)
    escribir_fasta([(n, "".join(s[c] for c in columnas)) for n, s in aln], destino)
    return len(columnas)


REGIONES = {"nterm": (1, 25), "resto": (31, 107)}  # posiciones de bartmosca; 26-30 (GGKHK) es igual en ambas fuentes


def arboles_por_region(modelo="DAYHOFF+G4"):
    """Si bartmosca es una quimera, cada región debería contar una historia distinta."""
    sys.path.insert(0, str(AQUI))
    aln4 = alineamiento("set4")
    resultados = {}
    for region, (desde, hasta) in REGIONES.items():
        destino = IQTREE_DIR / f"set4_{region}.fasta"
        IQTREE_DIR.mkdir(parents=True, exist_ok=True)
        n = subalineamiento(aln4, "bartmosca", desde, hasta, destino)
        correr(f"set4_{region}", destino, "-m", modelo, "-B", "1000")
        arbol = leer_newick((IQTREE_DIR / f"set4_{region}.treefile").read_text())
        resultados[region] = (n, arbol)
    return resultados


def resumir_regiones():
    for region in REGIONES:
        ruta = IQTREE_DIR / f"set4_{region}.treefile"
        if not ruta.exists():
            continue
        arbol = leer_newick(ruta.read_text())
        # 4 taxones sin raíz: el único split interno dice quién va con quién.
        internos = [c for c in clados(arbol) if 1 < len(c) < 4]
        hojas = set(arbol.hojas())
        for c in internos:
            lado = c if "bartmosca" in c else frozenset(hojas - c)
            print(f"set4_{region}: bartmosca agrupa con {sorted(lado - {'bartmosca'})} "
                  f"(UFBoot {soporte_de(arbol, c)})")


def graficar_todo():
    salida = []
    for prefijo, titulo in [
        ("set9_mfp", "9 secuencias - modelo ModelFinder, soporte SH-aLRT / UFBoot (1000)"),
        ("set9_poisson", "9 secuencias - modelo Poisson, soporte SH-aLRT / UFBoot (1000)"),
        ("set9_b100", "9 secuencias - modelo ModelFinder, bootstrap estándar (100)"),
    ]:
        arbol = enraizar(leer_newick((IQTREE_DIR / f"{prefijo}.treefile").read_text()), OUTGROUP)
        (IQTREE_DIR / f"{prefijo}_enraizado.nwk").write_text(arbol.to_newick() + "\n")
        dibujar(arbol, RESULTADOS / f"fig_arbol_{prefijo}.png", titulo + f"\nenraizado con {OUTGROUP}")
        salida.append((prefijo, arbol))
    return salida


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--solo-graficar", action="store_true")
    args = parser.parse_args(argv)
    if not args.solo_graficar:
        correr_todo()
    for prefijo, arbol in graficar_todo():
        print(f"{prefijo}: hermana(s) de bartmosca = {hermanas(arbol, 'bartmosca')} "
              f"(soporte {soporte_de(arbol, ['bartmosca', *hermanas(arbol, 'bartmosca')])})")
    resumir_regiones()
    return 0


if __name__ == "__main__":
    sys.exit(main())
