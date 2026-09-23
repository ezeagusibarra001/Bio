# TP: El juego de la vida

Resolución del TP [El juego de la vida](https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/blob/master/Teorico_Practicos/Intro_a_la_Biolog%C3%ADa/El_juego_de_la_vida.md) (Dra. Ana Julia Velez Rueda).

## Estructura

| Archivo | Qué resuelve |
|---|---|
| `codigo_genetico.py` | Tabla del código genético compartida por los scripts |
| `proteina_a_arn.py` | Desafío II: proteína → ARNm |
| `cajas_tata.py` | Desafío III: regiones promotoras delimitadas por cajas TATA |
| `juego_rpg.py` | Desafío IV: RPG de consola sobre la expresión génica |
| `datos/ejemplo_adn.fasta` | Secuencias de prueba para el Desafío III |
| `tests/` | Tests con pytest |

Todos los comandos se corren **desde esta carpeta** (`cd el_juego_de_la_vida`).

---

## 1. Por la célula

### 🤔 Para pensar: genes "informacionales" vs. "operacionales"

La frase se refiere a la **hipótesis de la quimera** sobre el origen de la célula eucariota. Al comparar genomas se vio que los genes eucariotas no provienen de un único ancestro:

- **Genes informacionales**: son los que manejan la información genética, es decir los involucrados en la **replicación del ADN, la transcripción y la traducción** (ARN polimerasas, proteínas ribosomales, factores de iniciación, histonas, etc.). En eucariotas estos genes se parecen mucho más a los de **arqueas**.
- **Genes operacionales**: son los del "mantenimiento" diario de la célula, es decir **metabolismo energético, biosíntesis de aminoácidos, lípidos y cofactores, transporte de membrana**, etc. Estos se parecen más a los de **bacterias**.

La interpretación más aceptada es que la célula eucariota surgió de la fusión/simbiosis entre una arquea hospedadora (que aportó la maquinaria informacional) y una o más bacterias (entre ellas la alfaproteobacteria que dio origen a la mitocondria), muchos de cuyos genes pasaron al núcleo por **transferencia génica endosimbiótica**. Esto conecta con la teoría endosimbiótica de Lynn Margulis que menciona el texto. Hoy esta idea se refuerza con el descubrimiento de las arqueas de Asgard, que son las parientes procariotas más cercanas a los eucariotas conocidas.

Además, los genes informacionales suelen participar en complejos con muchas interacciones y por eso se transfieren con menos facilidad entre linajes (la "hipótesis de la complejidad"), lo que hace que conserven mejor la señal de su origen arqueano.

### 🧗 Desafío I: diferencias entre célula procariota y eucariota

| Característica | Procariota | Eucariota |
|---|---|---|
| Núcleo | No tiene; el ADN está en el **nucleoide**, sin membrana | Núcleo delimitado por **envoltura nuclear** doble |
| Organelas membranosas | No tiene | Mitocondrias, RE, Golgi, lisosomas, cloroplastos (en plantas), etc. |
| Tamaño típico | 1 a 10 µm | 10 a 100 µm |
| Cromosomas | Generalmente uno, **circular** | Varios, **lineales**, con telómeros |
| Asociación del ADN | Proteínas tipo histona (en arqueas hay histonas) | ADN enrollado en **histonas** formando cromatina |
| Plásmidos | Frecuentes | Raros (algunas levaduras) |
| Intrones | Muy raros | Frecuentes; requieren **splicing** |
| Transcripción y traducción | **Acopladas**, ambas en el citoplasma | Separadas: transcripción en el núcleo, traducción en el citoplasma |
| ARNm | A menudo **policistrónico** (operones) | **Monocistrónico**, con caperuza 5' y cola poli-A |
| Ribosomas | **70S** (subunidades 50S + 30S) | **80S** (60S + 40S) en el citoplasma; 70S en mitocondrias y cloroplastos |
| Citoesqueleto | Proteínas homólogas simples (FtsZ, MreB) | Complejo: actina, microtúbulos, filamentos intermedios |
| División | **Fisión binaria** | **Mitosis** y **meiosis** |
| Pared celular | En bacterias, de **peptidoglicano** (las arqueas tienen otras) | Solo en plantas (celulosa) y hongos (quitina); las animales no tienen |
| Organización | Unicelulares (a veces colonias) | Uni o pluricelulares, con diferenciación celular |
| Reproducción sexual | No; intercambian genes por conjugación, transformación y transducción | Frecuente |

---

## 2. ¡Abracadabra una PROTEÍNA!

### 🤔 Para pensar: ¿qué cambiaría en procariotas?

- **No hay núcleo**: la transcripción ocurre en el citoplasma, por lo que el ARNm **no tiene que salir del núcleo**.
- **Transcripción y traducción acopladas**: los ribosomas empiezan a traducir el ARNm mientras todavía se está transcribiendo.
- **No hay splicing**: los genes casi no tienen intrones, así que no hay que eliminarlos.
- **Sin procesamiento del ARNm**: no se agrega caperuza 5' ni cola poli-A.
- **Una sola ARN polimerasa** (con distintos factores sigma) en lugar de las ARN pol I, II y III.
- **Otro promotor**: en bacterias no es la caja TATA con TBP sino las cajas **-10 (caja de Pribnow, `TATAAT`)** y **-35**, reconocidas por el factor sigma. (Las arqueas sí usan una caja TATA y TBP, ¡otra pista de su parentesco con eucariotas!)
- **Inicio de la traducción**: el ribosoma se ubica gracias a la secuencia **Shine-Dalgarno** y el primer aminoácido es **formil-metionina**.
- **ARNm policistrónicos**: un mismo ARNm puede codificar varias proteínas (operones).

### 🧗 Desafío II: de proteína a ARN (`proteina_a_arn.py`)

El código genético es **degenerado**: la mayoría de los aminoácidos tiene más de un codón. Por lo tanto, para una proteína existen muchísimos ARNm posibles (para Sec1 son ~1.4 × 10²⁰). El script invierte la tabla del código genético y permite elegir el primer codón o uno al azar, y verifica el resultado traduciéndolo de vuelta.

```bash
python proteina_a_arn.py ATVEKGGKHKTGPNEKGKKIFVQKCSQCHTVLHGLFGRKTGQA
python proteina_a_arn.py ATVEKGGKHKTGPNEKGKKIFVQKCSQCHTVLHGLFGRKTGQA --modo aleatorio --semilla 42
python proteina_a_arn.py ATVEK --inicio --stop   # agrega AUG y un codón stop
python proteina_a_arn.py MW --todas               # lista todas las opciones (péptidos cortos)
```

Salida para Sec1 (modo `primero`):
```
GCUACUGUUGAAAAAGGUGGUAAACAUAAAACUGGUCCUAAUGAAAAAGGUAAAAAAAUUUUUGUUCAAAAAUGUUCUCAAUGUCAUACUGUUUUACAUGGUUUAUUUGGUCGUAAAACUGGUCAAGCU
```

> Nota: Sec1 no empieza con M, así que es un fragmento interno de una proteína; por eso `--inicio` y `--stop` son opcionales.

### 🧗 Desafío III: cajas TATA (`cajas_tata.py`)

Toma un archivo FASTA (o de texto plano) y reporta como región promotora cada tramo que **empieza y termina con `TATAAA`**, con posiciones en base 1.

```bash
python cajas_tata.py datos/ejemplo_adn.fasta
python cajas_tata.py datos/ejemplo_adn.fasta --modo pares        # cajas de a pares sin compartir
python cajas_tata.py datos/ejemplo_adn.fasta --ambas-hebras      # busca también en la hebra reversa
python cajas_tata.py datos/ejemplo_adn.fasta --min-largo 30
```

Decisiones de diseño:
- Si hay 3 cajas, en modo `consecutivas` se reportan las regiones 1→2 y 2→3; en modo `pares`, solo 1→2.
- Con `--ambas-hebras` se busca en el reverso complementario, porque un gen puede estar en cualquiera de las dos hebras.
- Biológicamente, la caja TATA está ~25–30 pb antes del inicio de la transcripción; delimitar la región "entre dos cajas" es la simplificación que propone la consigna.

### 🧗 Desafío IV: "La odisea del gen" (`juego_rpg.py`)

Sos la **ARN Polimerasa II** y tenés que llevar un gen hasta proteína gastando lo menos posible de tu **ATP**. Si te quedás sin ATP, la célula entra en apoptosis. 💀

Niveles: **1. Promotor** (encontrar la caja TATA) → **2. Transcripción** (copiar la hebra molde) → **3. Splicing** (sacar el intrón) → **4. Exportación** (salir por el poro nuclear) → **5. Traducción** (codones → aminoácidos). En el medio aparecen eventos aleatorios con preguntas.

```bash
python juego_rpg.py
python juego_rpg.py --nombre Rosalind --dificultad dificil
python juego_rpg.py --semilla 7          # misma partida para toda la clase
python juego_rpg.py --help
```

Durante la partida: `pista` (cuesta 2 ATP), `ayuda`, `salir`. Los puzzles se generan al azar, así que cada partida es distinta; con `--semilla` se puede compartir la misma partida con la clase para competir por puntaje.

## Tests

```bash
python -m pytest tests -v
```
