# Atlas of Persistence

Exploratory mathematical and computational framework for studying physical
identity persistence, characteristic frequency, stability, gravitational
compactness, exergy, and accessible futures.

Marco matematico y computacional exploratorio para estudiar persistencia de
identidad fisica, frecuencia caracteristica, estabilidad, compactacion
gravitatoria, exergia y futuros accesibles.

## Archipelago Explorer

Archipelago Explorer is an exploratory visualization tool. It does not claim
to discover new particles, new laws, or hidden physical symmetries. It provides
an interface for testing whether frequency, lifetime, persistence cycles and
phase mappings reveal robust structures in known physical datasets.

Run the Dockerized browser app at `http://localhost:8000`:

```powershell
cd archipelago-explorer
docker compose up --build -d
```

See [the app README](archipelago-explorer/README.md) and
[scientific scope](archipelago-explorer/docs/scientific_scope.md).

## Documentation / Documentacion

- [English](README.en.md)
- [Espanol](README.es.md)

## License / Licencia

License: MIT. See [LICENSE](LICENSE).

Licencia: MIT. Consulte [LICENSE](LICENSE).

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
python -m src.harmonic_base_geometry
python -m src.stable_particle_sensitivity
python -m src.exploratory_3d_geometry
python -m src.cylindrical_helicoid_bases
python -m src.invariant_graph_validation
python -m unittest discover -s tests -v
```

The eleven analysis commands regenerate the CSV reports, summaries, plots, and
dynamic base-scan animations in `output/`, `outputs/complex_phase/`,
`outputs/gap_geometry/`, `outputs/base_scan/`, and
`outputs/invariant_graph_validation/`.

Los once comandos de analisis regeneran los informes CSV, resumenes y graficos
en `output/`, `outputs/complex_phase/`, `outputs/gap_geometry/` y
`outputs/invariant_graph_validation/`.

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

Invariant graph validation with explicit tolerances / Validacion de grafos
invariantes con tolerancias explicitas:
[docs/invariant_graph_validation.md](docs/invariant_graph_validation.md).

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

## Dynamic Base Scan / Barrido Dinamico de Base

The dynamic scan evaluates 5,000 evenly spaced logarithmic bases from `1.2`
to `50`. It tracks circular statistics, gaps, gap-center polygons, family and
interaction centroids, and persistent particle neighbors. The resulting
patterns are exploratory and do not establish a physically preferred base.

El barrido dinamico evalua 5.000 bases logaritmicas equiespaciadas entre `1.2`
y `50`. Sigue estadisticas circulares, huecos, poligonos de centros de huecos,
centroides por familia e interaccion y vecinos persistentes entre particulas.
Los patrones resultantes son exploratorios y no establecen una base fisicamente
preferida.

Generated numerical results / Resultados numericos generados:

- [Dynamic scan report / Informe del barrido dinamico](docs/dynamic_phase_structure_report.md)
- [Base scan metrics / Metricas del barrido de base](data/base_scan_metrics.csv)
- [Gap geometry metrics / Metricas geometricas de huecos](data/gap_geometry_metrics.csv)
- [Family and interaction alignment / Alineacion por familia e interaccion](data/family_alignment_metrics.csv)
- [Persistent neighbors / Vecinos persistentes](data/persistent_neighbors.csv)
- [Particle phase stability / Estabilidad de fase por particula](data/particle_phase_stability.csv)
- [Special-base metrics / Metricas de bases especiales](data/special_base_metrics.csv)

### Continuous Metrics / Metricas Continuas

| Largest gap / Mayor hueco | Entropy / Entropia |
| --- | --- |
| <img src="outputs/base_scan/base_vs_largest_gap.png" alt="Largest circular phase gap across logarithmic bases" width="420"> | <img src="outputs/base_scan/base_vs_entropy.png" alt="Circular phase entropy across logarithmic bases" width="420"> |

| Symmetry score / Simetria | Resultant length / Longitud resultante |
| --- | --- |
| <img src="outputs/base_scan/base_vs_symmetry_score.png" alt="Gap-center polygon symmetry score across logarithmic bases" width="420"> | <img src="outputs/base_scan/base_vs_resultant_length.png" alt="Circular resultant length across logarithmic bases" width="420"> |

### Persistent Relationships / Relaciones Persistentes

| Persistent neighbor network / Red de vecinos persistentes | Family centroid evolution / Evolucion de centroides por familia | Interaction centroid evolution / Evolucion de centroides por interaccion |
| --- | --- | --- |
| <img src="outputs/base_scan/persistent_neighbor_network.png" alt="Network of particle pairs that remain neighbors across logarithmic bases" width="300"> | <img src="outputs/base_scan/family_centroid_evolution.png" alt="Family centroid evolution across logarithmic bases" width="300"> | <img src="outputs/base_scan/interaction_centroid_evolution.png" alt="Interaction centroid evolution across logarithmic bases" width="300"> |

### Structured-Base Rankings / Clasificacion de Bases Estructuradas

| Most structured / Mas estructuradas | Least structured / Menos estructuradas |
| --- | --- |
| <img src="outputs/base_scan/top_10_most_structured_bases.png" alt="Top ten sampled bases by exploratory symmetry score" width="420"> | <img src="outputs/base_scan/top_10_least_structured_bases.png" alt="Bottom ten sampled bases by exploratory symmetry score" width="420"> |

### Dynamic Animation / Animacion Dinamica

![Dynamic persistence phase scan](outputs/base_scan/dynamic_phase_scan.gif)

[MP4 animation / Animacion MP4](outputs/base_scan/dynamic_phase_scan.mp4)

## Harmonic Base Geometry / Geometria Armonica de Bases

The harmonic geometry analysis detects local symmetry-score peaks, compares
them with simple mathematical expressions, computes an FFT spectrum, and
exports 2D and 3D geometry for the top 20 distinct peaks. Numerical proximity
does not establish a preferred constant or physical resonance.

El analisis de geometria armonica detecta maximos locales de simetria, los
compara con expresiones matematicas sencillas, calcula un espectro FFT y
exporta geometria 2D y 3D para los 20 maximos distintos principales. La
proximidad numerica no establece una constante preferida ni una resonancia
fisica.

- [Harmonic geometry report / Informe de geometria armonica](docs/harmonic_base_geometry.md)
- [Detected peaks / Maximos detectados](data/base_peaks.csv)
- [Constant matches / Coincidencias con constantes](data/peak_constant_matches.csv)
- [Gap geometry / Geometria de huecos](data/harmonic_gap_geometry.csv)
- [Persistent neighbors / Vecinos persistentes](data/harmonic_persistent_neighbors.csv)

### FFT Spectrum / Espectro FFT

![FFT spectrum of symmetry score](outputs/harmonic_base_geometry/fft_spectrum.png)

### 3D Animations / Animaciones 3D

![Rotating radial persistence spiral](outputs/harmonic_base_geometry/rotating_3d_view.gif)

![Base-sweep radial persistence spiral](outputs/harmonic_base_geometry/base_sweep_3d.gif)

[Rotating MP4](outputs/harmonic_base_geometry/rotating_3d_view.mp4) |
[Base-sweep MP4](outputs/harmonic_base_geometry/base_sweep_3d.mp4)

## Stable Particle Sensitivity / Sensibilidad de Particulas Estables

This sensitivity test compares finite-lifetime-only phases with controlled
stable-particle truncations. The truncation values are numerical assumptions,
not measured lifetimes.

Este test de sensibilidad compara fases con vidas medias finitas con
truncamientos controlados para particulas estables. Los valores de truncamiento
son supuestos numericos, no vidas medias medidas.

- [Sensitivity report / Informe de sensibilidad](docs/stable_particle_sensitivity.md)
- [Sensitivity summary / Resumen de sensibilidad](data/stable_sensitivity_summary.csv)
- [Centroid table / Tabla de centroides](data/stable_sensitivity_centroids.csv)
- [Particle phases / Fases de particulas](data/stable_sensitivity_phases.csv)
- [Continuous comparison / Comparacion continua](data/stable_sensitivity_base_scan.csv)

| Family centroids: finite only / Centroides: solo finitas | Family centroids: stable truncation / Centroides: truncamiento estable |
| --- | --- |
| <img src="outputs/stable_sensitivity/family_centroids_finite_only.png" alt="Family centroids using finite-lifetime particles only" width="420"> | <img src="outputs/stable_sensitivity/family_centroids_stable_truncated.png" alt="Family centroids with stable particles truncated at 1e35 seconds" width="420"> |

| Centroid difference / Diferencia de centroides | Base-scan comparison / Comparacion del barrido |
| --- | --- |
| <img src="outputs/stable_sensitivity/family_centroid_difference.png" alt="Family centroid difference between stable-truncated and finite-only variants" width="420"> | <img src="outputs/stable_sensitivity/base_scan_comparison.png" alt="Top-4 symmetry comparison for stable-truncated and finite-only variants" width="420"> |

| Unit circles / Circulos unidad | Lepton trajectories / Trayectorias de leptones |
| --- | --- |
| <img src="outputs/stable_sensitivity/unit_circles_variants.png" alt="Unit-circle comparison for finite-only and stable-truncated variants" width="420"> | <img src="outputs/stable_sensitivity/lepton_phase_trajectories.png" alt="Electron, muon, and tau phase trajectories under stable truncation" width="420"> |

## Exploratory 3D Geometry / Geometria 3D Exploratoria

This study investigates why helix-like views appear when wrapped logarithmic
phases are plotted against persistence. It searches for unexpected geometry
without treating visual patterns as physical evidence.

Este estudio investiga por que aparecen vistas helicoidales cuando las fases
logaritmicas envueltas se representan frente a la persistencia. Busca geometria
inesperada sin tratar los patrones visuales como evidencia fisica.

- [Exploratory geometry report / Informe de geometria exploratoria](docs/exploratory_geometry.md)
- [Reference coordinates / Coordenadas de referencia](data/exploratory_geometry_coordinates.csv)
- [Density cloud / Nube de densidad](data/exploratory_density_cloud.csv)
- [Gap constellations / Constelaciones de huecos](data/exploratory_gap_constellations.csv)
- [Interactive radial viewer / Visor radial interactivo](outputs/exploratory_3d/interactive_radial_geometry.html)

| Cylindrical helix / Helice cilindrica | Radial helix / Helice radial | Conical helix / Helice conica |
| --- | --- | --- |
| <img src="outputs/exploratory_3d/cylindrical.png" alt="Cylindrical persistence helix" width="300"> | <img src="outputs/exploratory_3d/radial.png" alt="Radial persistence helix" width="300"> | <img src="outputs/exploratory_3d/conical.png" alt="Conical persistence helix" width="300"> |

| Toroidal projection / Proyeccion toroidal | Spherical projection / Proyeccion esferica | Randomized controls / Controles aleatorios |
| --- | --- | --- |
| <img src="outputs/exploratory_3d/toroidal.png" alt="Toroidal persistence projection" width="300"> | <img src="outputs/exploratory_3d/spherical.png" alt="Spherical persistence projection" width="300"> | <img src="outputs/exploratory_3d/randomized_controls.png" alt="Cylindrical helix randomized controls" width="300"> |

![Connected gap-center constellations](outputs/exploratory_3d/gap_center_polygon_panel.png)

![Rotating radial helix](outputs/exploratory_3d/rotating_radial_helix.gif)

[Rotating MP4](outputs/exploratory_3d/rotating_radial_helix.mp4) |
[Base-sweep MP4](outputs/exploratory_3d/base_sweep_radial_helix.mp4)

## Cylindrical Helicoids Across Bases / Helicoides Cilindricos por Base

This exploratory comparison plots the same finite-lifetime entities on a
cylindrical helicoid for nine logarithmic bases. Height remains `log10(N)`;
base changes rotate particles around the reference cylinder.

Esta comparacion exploratoria representa las mismas entidades con vida media
finita sobre un helicoide cilindrico para nueve bases logaritmicas. La altura
permanece como `log10(N)`; cambiar la base rota las particulas alrededor del
cilindro de referencia.

- [Helicoid report / Informe de helicoides](docs/cylindrical_helicoid_bases.md)
- [Grouping metrics / Metricas de agrupacion](data/cylindrical_helicoid_metrics.csv)
- [Persistent neighbors / Vecinos persistentes](data/helicoid_persistent_neighbors.csv)
- [MATLAB/Wolfram coordinates / Coordenadas MATLAB/Wolfram](data/helicoid_coordinates_by_base.csv)

| All bases / Todas las bases | Grouping metrics / Metricas de agrupacion |
| --- | --- |
| <img src="outputs/cylindrical_helicoid/all_bases_grid.png" alt="Cylindrical helicoid comparison across logarithmic bases" width="420"> | <img src="outputs/cylindrical_helicoid/base_comparison_metrics.png" alt="Cylindrical helicoid grouping metrics across logarithmic bases" width="420"> |

| base 2 | base pi | base 10 |
| --- | --- | --- |
| <img src="outputs/cylindrical_helicoid/helicoid_base_2.png" alt="Cylindrical helicoid for base 2" width="300"> | <img src="outputs/cylindrical_helicoid/helicoid_base_pi.png" alt="Cylindrical helicoid for base pi" width="300"> | <img src="outputs/cylindrical_helicoid/helicoid_base_10.png" alt="Cylindrical helicoid for base 10" width="300"> |

## Golden Ratio Spiral Robustness / Robustez de la Espiral Aurea

This falsification-oriented analysis tests whether base `phi` produces unusually
coherent radial or helical geometry compared with nearby and unrelated bases.
The golden ratio is not assumed physically relevant. A visual spiral, a fit, or
the exploratory combined score does not imply physical significance.

Este analisis orientado a la falsacion comprueba si la base `phi` produce una
geometria radial o helicoidal inusualmente coherente frente a bases cercanas y
no relacionadas. No se supone que la proporcion aurea tenga relevancia fisica.
Una espiral visual, un ajuste o la puntuacion exploratoria combinada no implican
significado fisico.

- [Golden ratio spiral report / Informe de la espiral aurea](docs/golden_ratio_spiral_analysis.md)
- [Base metrics / Metricas por base](data/golden_ratio_base_metrics.csv)
- [Randomized controls / Controles aleatorios](data/golden_ratio_control_comparison.csv)
- [Ordering sensitivity / Sensibilidad al orden](data/golden_ratio_phi_ordering_fits.csv)
- [Stable-particle modes / Modos para particulas estables](data/golden_ratio_stable_mode_metrics.csv)

### Base Scan / Barrido de Bases

| Golden-spiral error / Error de espiral aurea | Spiral fit R2 / R2 del ajuste | Heuristic coherence / Coherencia heuristica |
| --- | --- | --- |
| <img src="outputs/golden_ratio_spiral/base_vs_golden_spiral_error.png" alt="Golden spiral error across logarithmic bases" width="300"> | <img src="outputs/golden_ratio_spiral/base_vs_spiral_fit_R2.png" alt="Radial spiral fit R2 across logarithmic bases" width="300"> | <img src="outputs/golden_ratio_spiral/base_vs_geometry_coherence_score.png" alt="Exploratory geometry coherence score across logarithmic bases" width="300"> |

| Phi neighborhood / Entorno de phi | Real versus randomized controls / Datos reales frente a controles aleatorios |
| --- | --- |
| <img src="outputs/golden_ratio_spiral/phi_neighborhood_zoom.png" alt="Exploratory score near phi" width="420"> | <img src="outputs/golden_ratio_spiral/real_vs_randomized_controls.png" alt="Real phi score compared with randomized controls" width="420"> |

### Helicoid Views / Vistas Helicoidales

| Radial phi | Radial pi | Radial e | Radial 10 |
| --- | --- | --- | --- |
| <img src="outputs/golden_ratio_spiral/radial_helicoid_phi.png" alt="Radial helicoid for base phi" width="220"> | <img src="outputs/golden_ratio_spiral/radial_helicoid_pi.png" alt="Radial helicoid for base pi" width="220"> | <img src="outputs/golden_ratio_spiral/radial_helicoid_e.png" alt="Radial helicoid for base e" width="220"> | <img src="outputs/golden_ratio_spiral/radial_helicoid_10.png" alt="Radial helicoid for base 10" width="220"> |

| Cylindrical phi / Cilindrica phi | Conical phi / Conica phi | Top-down phi / Vista superior phi |
| --- | --- | --- |
| <img src="outputs/golden_ratio_spiral/cylindrical_helicoid_phi.png" alt="Cylindrical helicoid for base phi" width="300"> | <img src="outputs/golden_ratio_spiral/conical_helicoid_phi.png" alt="Conical helicoid for base phi" width="300"> | <img src="outputs/golden_ratio_spiral/top_down_polar_projection_phi.png" alt="Top-down radial projection for base phi" width="300"> |

### Fit Diagnostics / Diagnostico del Ajuste

| Spiral fit overlay / Superposicion del ajuste | Residuals / Residuos |
| --- | --- |
| <img src="outputs/golden_ratio_spiral/spiral_fit_overlay_phi.png" alt="Exploratory logarithmic spiral fit overlay for base phi" width="420"> | <img src="outputs/golden_ratio_spiral/spiral_fit_residuals_phi.png" alt="Residuals for the exploratory phi spiral fit" width="420"> |

### Animations / Animaciones

![Animated base sweep around phi](outputs/golden_ratio_spiral/animated_base_sweep_around_phi.gif)

![Rotating radial helicoid for phi](outputs/golden_ratio_spiral/rotating_3d_phi.gif)

## Invariant Graph Validation / Validacion de Grafos Invariantes

This validation layer tests whether particle/entity relationships remain close
across logarithmic bases, coordinate representations, tolerance thresholds, and
randomized controls. It is explicitly exploratory: persistent graph structure
does not imply physical interaction.

Esta capa de validacion comprueba si las relaciones entre particulas o
entidades permanecen cercanas al cambiar bases logaritmicas, representaciones
de coordenadas, tolerancias y controles aleatorios. Es explicitamente
exploratoria: una estructura persistente de grafo no implica interaccion
fisica.

- [Validation report / Informe de validacion](docs/invariant_graph_validation.md)
- [Invariant edges / Aristas invariantes](data/invariant_edges.csv)
- [Graph validation summary / Resumen de validacion](data/graph_validation_summary.csv)
- [Robust clusters / Clusters robustos](data/robust_clusters.csv)
- [Tolerance configuration / Configuracion de tolerancias](config/tolerances.json)

Current conclusion so far: a small set of edges is persistent under the primary
10 degree angular tolerance across 2009 tested bases, including B0 / Bs0,
W / Z, B_plus / Bs0, and omega_782 / top. These are candidates for further
robustness testing, not discoveries. Several distance and nearest-neighbor
clusters look visually stable but do not pass randomized-control checks, so
they should be treated as weak candidates or likely artifacts.

Conclusion actual: un pequeno conjunto de aristas es persistente con la
tolerancia angular primaria de 10 grados en 2009 bases probadas, incluyendo
B0 / Bs0, W / Z, B_plus / Bs0 y omega_782 / top. Son candidatos para pruebas
adicionales de robustez, no descubrimientos. Varios clusters por distancia y
vecinos cercanos parecen estables visualmente, pero no superan los controles
aleatorios, por lo que deben tratarse como candidatos debiles o artefactos
probables.

| Persistent graph / Grafo persistente | Tolerance sweep / Barrido de tolerancia |
| --- | --- |
| <img src="outputs/invariant_graph_validation/01_persistent_neighbor_graph.png" alt="Persistent neighbor graph with explicit tolerance labels" width="420"> | <img src="outputs/invariant_graph_validation/02_tolerance_sweep_heatmap.png" alt="Tolerance sweep heatmap of persistence fractions" width="420"> |

| Base robustness / Robustez por base | Randomized controls / Controles aleatorios |
| --- | --- |
| <img src="outputs/invariant_graph_validation/03_base_robustness_heatmap.png" alt="Base robustness heatmap for graph edges" width="420"> | <img src="outputs/invariant_graph_validation/04_real_vs_randomized_graph_metric_plots.png" alt="Observed graph metrics compared with randomized controls" width="420"> |

| Family modularity / Modularidad por familia | Interaction modularity / Modularidad por interaccion |
| --- | --- |
| <img src="outputs/invariant_graph_validation/05_modularity_by_family_vs_controls.png" alt="Family modularity versus randomized controls" width="420"> | <img src="outputs/invariant_graph_validation/06_modularity_by_interaction_vs_controls.png" alt="Interaction modularity versus randomized controls" width="420"> |

| Spectral comparison / Comparacion espectral | Robust clusters / Clusters robustos |
| --- | --- |
| <img src="outputs/invariant_graph_validation/07_spectral_graph_comparison.png" alt="Spectral graph comparison against controls" width="420"> | <img src="outputs/invariant_graph_validation/08_robust_cluster_diagram.png" alt="Robust cluster diagram with persistence fractions" width="420"> |

| Edge survival / Supervivencia de aristas | 3D helicoid overlay / Superposicion helicoidal 3D |
| --- | --- |
| <img src="outputs/invariant_graph_validation/09_edge_survival_curve.png" alt="Edge survival curve by persistence fraction" width="420"> | <img src="outputs/invariant_graph_validation/10_3d_helicoid_with_edges.png" alt="3D helicoid with persistent graph edges overlaid" width="420"> |
