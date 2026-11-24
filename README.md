# UNQ · Bio

Monorepo con los trabajos prácticos de **Introducción a la Bioinformática** (UNQ), basados en el material de la [Dra. Ana Julia Velez Rueda](https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica).

| Proyecto | TP | Contenido |
|---|---|---|
| [`el_juego_de_la_vida/`](el_juego_de_la_vida/) | El juego de la vida | Célula procariota vs. eucariota, proteína → ARN, cajas TATA, RPG de expresión génica |
| [`biomoleculas/`](biomoleculas/) | Biomoléculas | Representación de proteínas, Rosalind Franklin, predicción de estructura secundaria, motivos en UniProt, sequence logos |
| [`alineamientos_secuenciales/`](alineamientos_secuenciales/) | Alineamientos secuenciales | Homología, identidad y gaps, Needleman-Wunsch (implementación propia), matrices de sustitución, BLAST |
| [`citocromo_c/`](citocromo_c/) | ¿En qué se parecen una gallina y una mosca? | Citocromo c de 9 organismos (NCBI), alineamiento múltiple, identidad y conservación, árboles UPGMA y NJ propios |
| [`inferencias_evolutivas/`](inferencias_evolutivas/) | Inferencias evolutivas | El caso bartmosca: mosaico/quimera, parálogos, alineamiento múltiple, IQ-TREE (ModelFinder, bootstrap), enraizado con outgroup |
| [`la_vida_en_tres_dimensiones/`](la_vida_en_tres_dimensiones/) | La vida en tres dimensiones | Formato PDB (ubiquitina), cavidades y pKa (EGFR, AKR1C4, albúminas), pLDDT de AlphaFold para E7 de HPV16 |

Cada carpeta es un proyecto independiente con su propio `README.md` (respuestas teóricas + cómo usar los scripts), sus datos de ejemplo y sus tests.

## Puesta en marcha

Requiere Python 3.9 o superior.

```bash
cd UNQ/Bio
python -m venv .venv
source .venv/bin/activate          # En Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Librerías externas: `biomoleculas/logo.py` usa logomaker, matplotlib y pandas; las figuras de los TPs nuevos usan matplotlib; `la_vida_en_tres_dimensiones/` usa pyKVFinder (cavidades) y propka (pKa). Algunos análisis llaman además a programas externos que no son de Python: MAFFT e IQ-TREE (`brew install mafft iqtree3`). Los tests no los necesitan: corren sobre los resultados guardados.

## Correr todos los tests

```bash
python -m pytest            # desde UNQ/Bio
```

## Agregar un TP nuevo

1. Crear una carpeta `nombre_del_tp/` con su `README.md`, sus scripts y una carpeta `tests/`.
2. Agregar la carpeta de tests en `testpaths` dentro de `pytest.ini`.
3. Sumar las dependencias nuevas a `requirements.txt` y el proyecto a la tabla de arriba.
