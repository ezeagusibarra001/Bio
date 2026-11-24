#!/usr/bin/env python3
"""Análisis del alineamiento múltiple de citocromo c (TP "¿En qué se parecen una
gallina y una mosca?").

Toma el alineamiento en formato Clustal (el que devuelve Clustal Omega o MAFFT
con --clustalout) y calcula todo lo que se pregunta en los PARA PENSAR:

  - identidad de a pares (matriz 9 x 9) y cuántas posiciones difieren entre
    humano y gallina;
  - la línea de conservación de Clustal (* : .) recalculada a mano, con las
    mismas reglas que usa Clustal, y cuántas columnas hay de cada tipo;
  - dos árboles de distancias implementados acá (UPGMA y Neighbor-Joining),
    para comparar con el árbol que da Clustal;
  - figuras: heatmap de identidad y árboles.

Solo usa la biblioteca estándar para el análisis; matplotlib se usa únicamente
para las figuras.

Ejemplos:
    python analisis_citocromo_c.py
    python analisis_citocromo_c.py --aln datos/citocromo_c_mafft.aln --sin-figuras
"""

import argparse
import sys
from pathlib import Path

GAP = "-"
CARPETA = Path(__file__).resolve().parent
ALN_POR_DEFECTO = CARPETA / "datos" / "citocromo_c_mafft.aln"

# Tabla del enunciado completa: accesión -> (nombre científico, nombre común)
ESPECIES = {
    "NP_061820.1": ("Homo sapiens", "Humano"),
    "NP_001072946.1": ("Gallus gallus", "Gallina"),
    "NP_001065289.1": ("Pan troglodytes", "Chimpancé"),
    "NP_001157486.1": ("Equus caballus", "Caballo"),
    "NP_001183974.1": ("Canis lupus familiaris", "Perro"),
    "AEP27192.1": ("Gorilla gorilla", "Gorila"),
    "XP_024245566.1": ("Oncorhynchus tshawytscha", "Salmón chinook"),
    "NP_001086101.1": ("Xenopus laevis", "Rana africana de uñas"),
    "NP_477164.1": ("Drosophila melanogaster", "Mosca de la fruta"),
}

# Grupos que usa Clustal para la línea de conservación (los mismos de ClustalW,
# que heredó Clustal Omega). Si todos los residuos de una columna caen dentro de
# un mismo grupo "fuerte" se marca ':'; si caen en un grupo "débil", '.'.
GRUPOS_FUERTES = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
GRUPOS_DEBILES = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK",
                  "NDEQHK", "NEQHRK", "FVLIM", "HFY"]

# Esquema de colores de Clustal Omega en el EBI ("Show Colors").
COLORES_CLUSTAL = {
    "AVFPMILW": "rojo (chicos e hidrofóbicos, incluye aromáticos menos Y)",
    "DE": "azul (ácidos)",
    "RK": "magenta (básicos, menos H)",
    "STYHCNGQ": "verde (hidroxilo, sulfhidrilo, amina y G)",
}


# ---------------------------------------------------------------- lectura
def leer_fasta(ruta):
    """Devuelve un dict {identificador: secuencia} respetando el orden del archivo.

    El identificador es la primera palabra del encabezado (la accesión de NCBI).
    """
    secuencias = {}
    actual = None
    for linea in Path(ruta).read_text().splitlines():
        linea = linea.strip()
        if not linea:
            continue
        if linea.startswith(">"):
            actual = linea[1:].split()[0]
            secuencias[actual] = []
        elif actual is not None:
            secuencias[actual].append(linea)
    return {k: "".join(v) for k, v in secuencias.items()}


def leer_clustal(ruta):
    """Lee un alineamiento en formato Clustal.

    Devuelve (alineamiento, conservacion): un dict {id: secuencia alineada} y la
    línea de conservación tal como la escribió el programa (puede ser '' si el
    archivo no la trae).
    """
    lineas = Path(ruta).read_text().splitlines()
    alineamiento = {}
    conservacion = []
    ancho_nombre = None
    for linea in lineas[1:]:  # la primera línea es el encabezado "CLUSTAL ..."
        if not linea.strip():
            continue
        if linea[0] != " ":
            partes = linea.split()
            nombre, bloque = partes[0], partes[1]
            alineamiento.setdefault(nombre, []).append(bloque)
            ancho_nombre = linea.index(bloque)
            largo_bloque = len(bloque)
        elif ancho_nombre is not None:
            # línea de conservación: misma columna de inicio que los bloques
            conservacion.append(linea[ancho_nombre:ancho_nombre + largo_bloque].ljust(largo_bloque))
    alineamiento = {k: "".join(v) for k, v in alineamiento.items()}
    largo = len(next(iter(alineamiento.values())))
    return alineamiento, "".join(conservacion).ljust(largo)[:largo]


# ---------------------------------------------------------------- conservación
def simbolo_conservacion(columna):
    """Símbolo de Clustal para una columna del alineamiento.

    '*' todos idénticos, ':' propiedades muy parecidas (grupo fuerte),
    '.' propiedades algo parecidas (grupo débil), ' ' sin conservación.
    Una columna con algún gap nunca se marca.
    """
    if GAP in columna:
        return " "
    residuos = set(columna)
    if len(residuos) == 1:
        return "*"
    if any(residuos <= set(g) for g in GRUPOS_FUERTES):
        return ":"
    if any(residuos <= set(g) for g in GRUPOS_DEBILES):
        return "."
    return " "


def linea_conservacion(alineamiento):
    seqs = list(alineamiento.values())
    return "".join(simbolo_conservacion([s[i] for s in seqs]) for i in range(len(seqs[0])))


def contar_conservacion(linea):
    return {
        "*": linea.count("*"),
        ":": linea.count(":"),
        ".": linea.count("."),
        "sin marca": linea.count(" "),
    }


def columnas_con_gap(alineamiento):
    seqs = list(alineamiento.values())
    return sum(1 for i in range(len(seqs[0])) if any(s[i] == GAP for s in seqs))


# ---------------------------------------------------------------- identidad
def identidad(a, b):
    """% de identidad entre dos secuencias alineadas, sobre las columnas donde
    ninguna de las dos tiene gap (mismo criterio que alineamientos_secuenciales)."""
    pares = [(x, y) for x, y in zip(a, b) if x != GAP and y != GAP]
    if not pares:
        return 0.0
    return 100 * sum(1 for x, y in pares if x == y) / len(pares)


def matriz_identidad(alineamiento):
    ids = list(alineamiento)
    return ids, [[identidad(alineamiento[i], alineamiento[j]) for j in ids] for i in ids]


def diferencias(a, b):
    """Lista de (posición en a, residuo en a, residuo en b) donde difieren.

    La posición se cuenta sobre la secuencia a sin gaps (empezando en 1, es
    decir contando la metionina inicial tal como aparece en NCBI). Si b tiene un
    gap en esa columna, aparece '-'.
    """
    salida = []
    pos = 0
    for x, y in zip(a, b):
        if x != GAP:
            pos += 1
        if x != y and not (x == GAP and y == GAP):
            salida.append((pos if x != GAP else None, x, y))
    return salida


# ---------------------------------------------------------------- árboles
class Nodo:
    """Nodo de un árbol: hoja (con nombre) o interno (con hijos y largos de rama)."""

    def __init__(self, nombre=None, hijos=None):
        self.nombre = nombre
        self.hijos = hijos or []  # lista de (Nodo, largo_de_rama)

    def es_hoja(self):
        return not self.hijos

    def hojas(self):
        if self.es_hoja():
            return [self.nombre]
        return [h for hijo, _ in self.hijos for h in hijo.hojas()]

    def newick(self, decimales=4):
        if self.es_hoja():
            return self.nombre.replace(" ", "_")
        partes = [f"{h.newick(decimales)}:{max(0.0, l):.{decimales}f}" for h, l in self.hijos]
        return "(" + ",".join(partes) + ")"


def distancias_p(matriz_id):
    """Distancia p (proporción de sitios distintos) a partir de identidades en %."""
    return [[round(1 - v / 100, 6) for v in fila] for fila in matriz_id]


def upgma(nombres, D):
    """Agrupamiento UPGMA (promedio). Devuelve la raíz (árbol enraizado y ultramétrico)."""
    grupos = [(Nodo(n), 1, 0.0) for n in nombres]  # (nodo, tamaño, altura)
    D = [fila[:] for fila in D]
    while len(grupos) > 1:
        n = len(grupos)
        i, j = min(((i, j) for i in range(n) for j in range(i + 1, n)), key=lambda p: D[p[0]][p[1]])
        (ni, si, hi), (nj, sj, hj) = grupos[i], grupos[j]
        altura = D[i][j] / 2
        nuevo = (Nodo(hijos=[(ni, altura - hi), (nj, altura - hj)]), si + sj, altura)
        nuevas = [(D[i][k] * si + D[j][k] * sj) / (si + sj) for k in range(n) if k not in (i, j)]
        resto = [k for k in range(n) if k not in (i, j)]
        D = [[D[a][b] for b in resto] for a in resto]
        for fila, d in zip(D, nuevas):
            fila.append(d)
        D.append(nuevas + [0.0])
        grupos = [grupos[k] for k in resto] + [nuevo]
    return grupos[0][0]


def neighbor_joining(nombres, D):
    """Neighbor-Joining (Saitou y Nei 1987). Devuelve un árbol SIN raíz,
    representado con una trifurcación en el nodo central."""
    nodos = [Nodo(n) for n in nombres]
    D = [fila[:] for fila in D]
    while len(nodos) > 3:
        n = len(nodos)
        r = [sum(fila) for fila in D]
        i, j = min(((i, j) for i in range(n) for j in range(i + 1, n)),
                   key=lambda p: (n - 2) * D[p[0]][p[1]] - r[p[0]] - r[p[1]])
        li = D[i][j] / 2 + (r[i] - r[j]) / (2 * (n - 2))
        lj = D[i][j] - li
        nuevo = Nodo(hijos=[(nodos[i], max(0.0, li)), (nodos[j], max(0.0, lj))])
        resto = [k for k in range(n) if k not in (i, j)]
        nuevas = [(D[i][k] + D[j][k] - D[i][j]) / 2 for k in resto]
        D = [[D[a][b] for b in resto] for a in resto]
        for fila, d in zip(D, nuevas):
            fila.append(d)
        D.append(nuevas + [0.0])
        nodos = [nodos[k] for k in resto] + [nuevo]
    a, b, c = nodos
    la = (D[0][1] + D[0][2] - D[1][2]) / 2
    lb = (D[0][1] + D[1][2] - D[0][2]) / 2
    lc = (D[0][2] + D[1][2] - D[0][1]) / 2
    return Nodo(hijos=[(a, max(0.0, la)), (b, max(0.0, lb)), (c, max(0.0, lc))])


def enraizar(raiz, grupo_externo):
    """Re-enraiza un árbol en el punto medio de la rama que lleva a la hoja
    `grupo_externo`. Sirve para leer el árbol NJ (que no tiene raíz) en el
    sentido del tiempo, usando la especie más lejana como outgroup."""
    vecinos = {}

    def conectar(a, b, largo):
        vecinos.setdefault(id(a), (a, []))[1].append((b, largo))
        vecinos.setdefault(id(b), (b, []))[1].append((a, largo))

    def recorrer(nodo):
        vecinos.setdefault(id(nodo), (nodo, []))
        for hijo, largo in nodo.hijos:
            conectar(nodo, hijo, largo)
            recorrer(hijo)

    recorrer(raiz)
    hoja = next(n for n, _ in vecinos.values() if n.es_hoja() and n.nombre == grupo_externo)
    padre, largo = vecinos[id(hoja)][1][0]

    def construir(nodo, desde):
        if nodo.es_hoja():
            return Nodo(nodo.nombre)
        hijos = [(construir(v, nodo), l) for v, l in vecinos[id(nodo)][1] if v is not desde]
        return Nodo(hijos=hijos)

    return Nodo(hijos=[(Nodo(hoja.nombre), largo / 2), (construir(padre, hoja), largo / 2)])


# ---------------------------------------------------------------- figuras
def dibujar_heatmap(nombres, matriz, ruta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 7))
    im = ax.imshow(matriz, cmap="viridis", vmin=min(min(f) for f in matriz), vmax=100)
    ax.set_xticks(range(len(nombres)), nombres, rotation=45, ha="right")
    ax.set_yticks(range(len(nombres)), nombres)
    for i, fila in enumerate(matriz):
        for j, v in enumerate(fila):
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    color="black" if v > 90 else "white", fontsize=9)
    fig.colorbar(im, ax=ax, label="% identidad")
    ax.set_title("Citocromo c: identidad de a pares (%)")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def _dibujar_arbol(ax, raiz, titulo):
    y_hojas = {}
    contador = [0]

    def posiciones(nodo, x):
        if nodo.es_hoja():
            y = contador[0]
            contador[0] += 1
            ax.text(x + 0.003, y, " " + nodo.nombre, va="center", fontsize=9)
            y_hojas[id(nodo)] = y
            return y
        ys = []
        for hijo, largo in nodo.hijos:
            yh = posiciones(hijo, x + largo)
            ax.plot([x, x + largo], [yh, yh], color="black", lw=1.2)
            ys.append(yh)
        ax.plot([x, x], [min(ys), max(ys)], color="black", lw=1.2)
        return sum(ys) / len(ys)

    posiciones(raiz, 0.0)
    ax.set_title(titulo)
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_xlabel("distancia p (proporción de sitios distintos)")
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.set_xlim(right=ax.get_xlim()[1] * 1.45)


def dibujar_arboles(arbol_upgma, arbol_nj, ruta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    _dibujar_arbol(a1, arbol_upgma, "UPGMA (enraizado, ultramétrico)")
    _dibujar_arbol(a2, arbol_nj, "Neighbor-Joining (raíz en la mosca)")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- reporte
def analizar(ruta_aln):
    aln, cons_programa = leer_clustal(ruta_aln)
    comunes = {acc: ESPECIES.get(acc, (acc, acc))[1] for acc in aln}
    aln_comun = {comunes[k]: v for k, v in aln.items()}
    nombres, mid = matriz_identidad(aln_comun)
    D = distancias_p(mid)
    arbol_upgma = upgma(nombres, D)
    arbol_nj = enraizar(neighbor_joining(nombres, D), "Mosca de la fruta")
    cons = linea_conservacion(aln)
    return {
        "alineamiento": aln_comun,
        "conservacion_programa": cons_programa,
        "conservacion": cons,
        "conteo": contar_conservacion(cons),
        "columnas_con_gap": columnas_con_gap(aln),
        "nombres": nombres,
        "identidad": mid,
        "humano_vs_gallina": diferencias(aln_comun["Humano"], aln_comun["Gallina"]),
        "upgma": arbol_upgma,
        "nj": arbol_nj,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Análisis del alineamiento de citocromo c.")
    parser.add_argument("--aln", default=str(ALN_POR_DEFECTO), help="Alineamiento en formato Clustal")
    parser.add_argument("--salida", default=str(CARPETA / "resultados"), help="Carpeta para las figuras")
    parser.add_argument("--sin-figuras", action="store_true", help="No generar los PNG")
    args = parser.parse_args(argv)

    r = analizar(args.aln)
    aln = r["alineamiento"]
    largo = len(next(iter(aln.values())))

    print(f"Alineamiento: {args.aln}  ({len(aln)} secuencias, {largo} columnas)\n")
    for nombre, seq in aln.items():
        print(f"{nombre:<22}{seq}")
    print(f"{'':<22}{r['conservacion']}\n")

    c = r["conteo"]
    print(f"Columnas: '*' = {c['*']}  ':' = {c[':']}  '.' = {c['.']}  sin marca = {c['sin marca']}"
          f"  (con gap: {r['columnas_con_gap']})")
    if r["conservacion_programa"].strip():
        iguales = r["conservacion_programa"] == r["conservacion"]
        print(f"¿La línea recalculada coincide con la del programa? {'sí' if iguales else 'NO'}")

    print("\nIdentidad de a pares (%):")
    nombres = r["nombres"]
    cortos = [n[:8] for n in nombres]
    print(" " * 22 + "".join(f"{n:>9}" for n in cortos))
    for n, fila in zip(nombres, r["identidad"]):
        print(f"{n:<22}" + "".join(f"{v:9.1f}" for v in fila))

    difs = r["humano_vs_gallina"]
    print(f"\nHumano vs gallina: {len(difs)} posiciones distintas")
    print("  " + ", ".join(f"{a}{p}{b}" for p, a, b in difs))

    print("\nUPGMA:", r["upgma"].newick() + ";")
    print("NJ (raíz en la mosca):", r["nj"].newick() + ";")

    if not args.sin_figuras:
        salida = Path(args.salida)
        salida.mkdir(exist_ok=True)
        # ordeno filas/columnas como las hojas del árbol NJ para que se vean los grupos
        orden = [nombres.index(h) for h in r["nj"].hojas()]
        mid = [[r["identidad"][i][j] for j in orden] for i in orden]
        dibujar_heatmap([nombres[i] for i in orden], mid, salida / "heatmap_identidad.png")
        dibujar_arboles(r["upgma"], r["nj"], salida / "arboles.png")
        print(f"\nFiguras guardadas en {salida}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
