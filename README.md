# Atlas of Persistence

Exploratory mathematical and computational framework for studying physical
identity persistence, characteristic frequency, stability, gravitational
compactness, exergy, and accessible futures.

Marco matematico y computacional exploratorio para estudiar persistencia de
identidad fisica, frecuencia caracteristica, estabilidad, compactacion
gravitatoria, exergia y futuros accesibles.

## Documentation / Documentacion

- [English](README.en.md)
- [Espanol](README.es.md)

## Quick Start / Inicio Rapido

Requires Python 3.10 or later. Install the analysis dependencies first.

Requiere Python 3.10 o posterior. Instale primero las dependencias de analisis.

```powershell
python -m pip install -r requirements.txt
python -m src.atlas_persistence
python -m src.persistence
python -m src.complex_phase
python -m unittest discover -s tests -v
```

The first two commands regenerate the CSV reports, summaries, and persistence
landscape SVG plots in `output/`.

Los dos primeros comandos regeneran los informes CSV, resumenes y graficos SVG
del paisaje de persistencia en `output/`.

Focused framework documentation / Documentacion del marco enfocado:
[docs/persistence_framework.md](docs/persistence_framework.md).

Complex phase-gap analysis / Analisis complejo de huecos de fase:
[docs/complex_phase_gap_analysis.md](docs/complex_phase_gap_analysis.md).

## Persistence Landscape / Paisaje de Persistencia

The focused framework investigates whether the dimensionless index
`N = f * tau` reveals useful descriptive structures across measured physical
entities. For the particle sample, `f` is the Compton frequency and `tau` is a
measured lifetime or experimental lower bound.

El marco enfocado investiga si el indice adimensional `N = f * tau` revela
estructuras descriptivas utiles entre entidades fisicas medidas. Para la
muestra de particulas, `f` es la frecuencia de Compton y `tau` es una vida
media medida o una cota inferior experimental.

The current table includes electron, muon, tau, neutron, proton, W, Z, and
Higgs. Electron and proton persistence values are lower bounds. The plots are
exploratory visualizations, not evidence of a discovered physical law or
stability island.

La tabla actual incluye electron, muon, tau, neutron, proton, W, Z y Higgs. Los
valores de persistencia del electron y del proton son cotas inferiores. Los
graficos son visualizaciones exploratorias, no evidencia de una ley fisica o
una isla de estabilidad descubierta.

Generated numerical results / Resultados numericos generados:

- [Persistence table / Tabla de persistencia](data/persistence_table.csv)
- [Landscape comparison / Comparacion de paisajes](output/persistence/landscape_metrics.csv)
- [Largest gaps / Mayores huecos](output/persistence/persistence_gaps.csv)
- [Exploratory clusters / Clusters exploratorios](output/persistence/persistence_clusters.csv)
- [Generated summary / Resumen generado](output/persistence/summary.md)

### Sorted Persistence Index / Indice de Persistencia Ordenado

![Sorted log10 persistence index](output/persistence/a_log10_n_sorted.svg)

### Frequency vs Lifetime / Frecuencia vs Vida Media

![Frequency versus lifetime on logarithmic scales](output/persistence/b_frequency_vs_lifetime.svg)

### Frequency vs Persistence / Frecuencia vs Persistencia

![Frequency versus persistence index on logarithmic scales](output/persistence/c_frequency_vs_persistence.svg)

### Adjacent-Gap Histogram / Histograma de Huecos Adyacentes

![Histogram of adjacent persistence gaps](output/persistence/d_gap_histogram.svg)

### Exploratory Clusters / Clusters Exploratorios

![Exploratory adjacent-gap clusters](output/persistence/e_cluster_landscape.svg)
