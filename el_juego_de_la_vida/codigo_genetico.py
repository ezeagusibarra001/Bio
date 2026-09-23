"""Código genético universal (ARNm -> aminoácido, código de una letra).

Módulo compartido por los scripts del TP "El juego de la vida".
"""

# Codón de ARNm -> aminoácido. '*' indica codón de stop.
CODON_A_AMINOACIDO = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Tabla inversa: aminoácido -> lista de codones que lo codifican.
AMINOACIDO_A_CODONES = {}
for _codon, _aa in CODON_A_AMINOACIDO.items():
    AMINOACIDO_A_CODONES.setdefault(_aa, []).append(_codon)

CODON_INICIO = "AUG"
CODONES_STOP = AMINOACIDO_A_CODONES["*"]

NOMBRES_AMINOACIDOS = {
    "A": "Alanina", "R": "Arginina", "N": "Asparagina", "D": "Ácido aspártico",
    "C": "Cisteína", "Q": "Glutamina", "E": "Ácido glutámico", "G": "Glicina",
    "H": "Histidina", "I": "Isoleucina", "L": "Leucina", "K": "Lisina",
    "M": "Metionina", "F": "Fenilalanina", "P": "Prolina", "S": "Serina",
    "T": "Treonina", "W": "Triptófano", "Y": "Tirosina", "V": "Valina",
    "*": "STOP",
}


def traducir(arn):
    """Traduce una secuencia de ARNm a proteína (se detiene en el primer stop)."""
    arn = arn.upper().replace("T", "U")
    proteina = []
    for i in range(0, len(arn) - 2, 3):
        aa = CODON_A_AMINOACIDO.get(arn[i:i + 3], "X")
        if aa == "*":
            break
        proteina.append(aa)
    return "".join(proteina)
