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
python -m src.gap_geometry
python -m src.gap_pair_structure
python -m src.dynamic_phase_structure
python -m unittest discover -s tests -v
```

The six analysis commands regenerate the CSV reports, summaries, plots, and
dynamic base-scan animations in `output/`, `outputs/complex_phase/`,
`outputs/gap_geometry/`, and `outputs/base_scan/`.

Los cinco comandos de analisis regeneran los informes CSV, resumenes y graficos
en `output/`, `outputs/complex_phase/` y `outputs/gap_geometry/`.

Focused framework documentation / Documentacion del marco enfocado:
[docs/persistence_framework.md](docs/persistence_framework.md).

Complex phase-gap analysis / Analisis complejo de huecos de fase:
[docs/complex_phase_gap_analysis.md](docs/complex_phase_gap_analysis.md).

Gap-center polygon analysis / Analisis poligonal de centros de huecos:
[docs/gap_geometry_analysis.md](docs/gap_geometry_analysis.md).

Gap-pair structure analysis / Analisis estructural de pares de huecos:
[docs/gap_pair_structure_analysis.md](docs/gap_pair_structure_analysis.md).

Dynamic base and phase-structure scan / Barrido dinamico de base y estructura
de fase:
[docs/dynamic_phase_structure_report.md](docs/dynamic_phase_structure_report.md).

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

## Complex Phase Landscape / Paisaje de Fase Compleja

The complex phase branch maps finite particle persistence cycles onto the unit
circle for bases `2`, `e`, `pi`, `phi`, and `10`. This is an exploratory
numerical representation. Angular alignment does not imply a physical law,
preferred base, or resonance mechanism.

La rama de fase compleja proyecta ciclos de persistencia finitos sobre el
circulo unidad para las bases `2`, `e`, `pi`, `phi` y `10`. Es una
representacion numerica exploratoria. La alineacion angular no implica una ley
fisica, una base preferida ni un mecanismo de resonancia.

Generated numerical results / Resultados numericos generados:

- [Complex phase table / Tabla de fase compleja](data/complex_phase_table.csv)
- [Circular gaps / Huecos circulares](data/complex_phase_gaps.csv)
- [Circular metrics / Metricas circulares](data/complex_phase_summary.csv)
- [Analysis notes / Notas del analisis](docs/complex_phase_gap_analysis.md)
- [Gap-center geometry / Geometria de centros de huecos](data/gap_geometry_summary.csv)
- [Gap-center notes / Notas de centros de huecos](docs/gap_geometry_analysis.md)
- [Gap-pair metrics / Metricas de pares de huecos](data/gap_pair_structure.csv)
- [Gap-pair scale summary / Resumen de escalas de pares](data/gap_pair_structure_summary.csv)
- [Gap-pair frequency structure / Estructura frecuencial de pares](data/gap_pair_frequency_structure.csv)
- [Gap-pair notes / Notas de pares de huecos](docs/gap_pair_structure_analysis.md)

### Unit Circles / Circulos Unidad

| base 2 | base e |
| --- | --- |
| <img src="outputs/complex_phase/unit_circle_base_2.png" alt="Complex persistence phase on the unit circle for base 2" width="420"> | <img src="outputs/complex_phase/unit_circle_base_e.png" alt="Complex persistence phase on the unit circle for base e" width="420"> |

| base pi | base phi | base 10 |
| --- | --- | --- |
| <img src="outputs/complex_phase/unit_circle_base_pi.png" alt="Complex persistence phase on the unit circle for base pi" width="300"> | <img src="outputs/complex_phase/unit_circle_base_phi.png" alt="Complex persistence phase on the unit circle for base phi" width="300"> | <img src="outputs/complex_phase/unit_circle_base_10.png" alt="Complex persistence phase on the unit circle for base 10" width="300"> |

### Polar Landscapes / Paisajes Polares

| base 2 | base e |
| --- | --- |
| <img src="outputs/complex_phase/polar_base_2.png" alt="Normalized persistence radius versus phase for base 2" width="420"> | <img src="outputs/complex_phase/polar_base_e.png" alt="Normalized persistence radius versus phase for base e" width="420"> |

| base pi | base phi | base 10 |
| --- | --- | --- |
| <img src="outputs/complex_phase/polar_base_pi.png" alt="Normalized persistence radius versus phase for base pi" width="300"> | <img src="outputs/complex_phase/polar_base_phi.png" alt="Normalized persistence radius versus phase for base phi" width="300"> | <img src="outputs/complex_phase/polar_base_10.png" alt="Normalized persistence radius versus phase for base 10" width="300"> |

### Phase Comparisons / Comparaciones de Fase

| Frequency vs phase / Frecuencia vs fase | Persistence vs phase / Persistencia vs fase |
| --- | --- |
| <img src="outputs/complex_phase/frequency_vs_phase.png" alt="Compton frequency versus complex persistence phase" width="420"> | <img src="outputs/complex_phase/persistence_vs_phase.png" alt="Persistence cycles versus complex persistence phase" width="420"> |

| Largest circular gaps / Mayores huecos circulares | Phase-fraction heatmap / Mapa de calor de fraccion de fase |
| --- | --- |
| <img src="outputs/complex_phase/gap_ranking_top10.png" alt="Largest circular gaps by logarithmic base" width="420"> | <img src="outputs/complex_phase/phase_fraction_heatmap.png" alt="Particle phase fractions across logarithmic bases" width="420"> |

## Gap-Center Geometry / Geometria de Centros de Huecos

The gap-center geometry analysis connects the centers of the largest circular
phase gaps in angular order. These polygons are exploratory visual summaries.
Recognizable geometry does not imply a physical law or a preferred base.

El analisis de geometria de centros de huecos conecta en orden angular los
centros de los mayores huecos circulares de fase. Estos poligonos son resumenes
visuales exploratorios. Una geometria reconocible no implica una ley fisica ni
una base preferida.

### Polygon Gallery / Galeria de Poligonos

#### Base 2

| top 3 | top 4 | top 5 | top 6 |
| --- | --- | --- | --- |
| <img src="outputs/gap_geometry/gap_center_polygon_base_2_top_3.png" alt="Gap-center polygon for base 2 and top 3 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_2_top_4.png" alt="Gap-center polygon for base 2 and top 4 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_2_top_5.png" alt="Gap-center polygon for base 2 and top 5 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_2_top_6.png" alt="Gap-center polygon for base 2 and top 6 gaps" width="220"> |

#### Base e

| top 3 | top 4 | top 5 | top 6 |
| --- | --- | --- | --- |
| <img src="outputs/gap_geometry/gap_center_polygon_base_e_top_3.png" alt="Gap-center polygon for base e and top 3 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_e_top_4.png" alt="Gap-center polygon for base e and top 4 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_e_top_5.png" alt="Gap-center polygon for base e and top 5 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_e_top_6.png" alt="Gap-center polygon for base e and top 6 gaps" width="220"> |

#### Base pi

| top 3 | top 4 | top 5 | top 6 |
| --- | --- | --- | --- |
| <img src="outputs/gap_geometry/gap_center_polygon_base_pi_top_3.png" alt="Gap-center polygon for base pi and top 3 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_pi_top_4.png" alt="Gap-center polygon for base pi and top 4 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_pi_top_5.png" alt="Gap-center polygon for base pi and top 5 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_pi_top_6.png" alt="Gap-center polygon for base pi and top 6 gaps" width="220"> |

#### Base phi

| top 3 | top 4 | top 5 | top 6 |
| --- | --- | --- | --- |
| <img src="outputs/gap_geometry/gap_center_polygon_base_phi_top_3.png" alt="Gap-center polygon for base phi and top 3 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_phi_top_4.png" alt="Gap-center polygon for base phi and top 4 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_phi_top_5.png" alt="Gap-center polygon for base phi and top 5 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_phi_top_6.png" alt="Gap-center polygon for base phi and top 6 gaps" width="220"> |

#### Base 10

| top 3 | top 4 | top 5 | top 6 |
| --- | --- | --- | --- |
| <img src="outputs/gap_geometry/gap_center_polygon_base_10_top_3.png" alt="Gap-center polygon for base 10 and top 3 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_10_top_4.png" alt="Gap-center polygon for base 10 and top 4 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_10_top_5.png" alt="Gap-center polygon for base 10 and top 5 gaps" width="220"> | <img src="outputs/gap_geometry/gap_center_polygon_base_10_top_6.png" alt="Gap-center polygon for base 10 and top 6 gaps" width="220"> |

### Geometry Comparisons / Comparaciones de Geometria

| Symmetry score / Simetria | Polygon area / Area poligonal | Best symmetry / Mejor simetria |
| --- | --- | --- |
| <img src="outputs/gap_geometry/symmetry_score_by_base_and_k.png" alt="Gap-center symmetry score by base and top K" width="300"> | <img src="outputs/gap_geometry/polygon_area_by_base_and_k.png" alt="Gap-center polygon area by base and top K" width="300"> | <img src="outputs/gap_geometry/best_symmetry_score_radar.png" alt="Best gap-center symmetry score by base" width="300"> |
