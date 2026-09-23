# UNQ · Bio

Monorepo con los trabajos prácticos de **Introducción a la Bioinformática** (UNQ), basados en el material de la [Dra. Ana Julia Velez Rueda](https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica).

| Proyecto | TP | Contenido |
|---|---|---|
| [`el_juego_de_la_vida/`](el_juego_de_la_vida/) | El juego de la vida | Célula procariota vs. eucariota, proteína → ARN, cajas TATA, RPG de expresión génica |
| [`biomoleculas/`](biomoleculas/) | Biomoléculas | Representación de proteínas, Rosalind Franklin, predicción de estructura secundaria, motivos en UniProt, sequence logos |

Cada carpeta es un proyecto independiente con su propio `README.md` (respuestas teóricas + cómo usar los scripts), sus datos de ejemplo y sus tests.

## Puesta en marcha

Requiere Python 3.9 o superior.

```bash
cd UNQ/Bio
python -m venv .venv
source .venv/bin/activate          # En Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Solo `biomoleculas/logo.py` necesita librerías externas (logomaker, matplotlib, pandas). El resto usa únicamente la biblioteca estándar.

## Correr todos los tests

```bash
python -m pytest            # desde UNQ/Bio
```

## Agregar un TP nuevo

1. Crear una carpeta `nombre_del_tp/` con su `README.md`, sus scripts y una carpeta `tests/`.
2. Agregar la carpeta de tests en `testpaths` dentro de `pytest.ini`.
3. Sumar las dependencias nuevas a `requirements.txt` y el proyecto a la tabla de arriba.
