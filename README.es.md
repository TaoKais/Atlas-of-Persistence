# Atlas de Persistencia

[English](README.en.md) | **Espanol**

Marco matematico y computacional exploratorio para estudiar persistencia de
identidad, frecuencia caracteristica, estabilidad, compactacion gravitatoria y
exergia. El proyecto genera hipotesis contrastables y calculos reproducibles;
no propone una teoria fisica demostrada.

## Resultado conceptual minimo

La frecuencia de Compton

```text
f_C = m c^2 / h
```

es proporcional a la masa. Por tanto, ordenar particulas por `f_C` o por masa
produce exactamente el mismo orden y los mismos huecos logaritmicos. `f_C` es
una coordenada fisica valida, pero no puede clasificar la persistencia mejor que
la masa si se usa de forma aislada.

El atlas explora parametros que si incorporan informacion adicional:

| Parametro | Definicion | Interpretacion exploratoria |
| --- | --- | --- |
| `N_C` | `f_C * tau` | ciclos de Compton durante la vida media |
| `Q_C` | `2 pi f_C tau` | inverso de la anchura relativa `Gamma/(mc^2)` |
| `rho_E` | `E_total/(mc^2)` | energia total relativista respecto al reposo |
| `Phi` | `GM/(Rc^2)` | compactacion gravitatoria |
| `delta_X` | `X_destruida/X_entrada` | fraccion de exergia destruida |

`N_C` y `Q_C` no son explicaciones causales de la estabilidad: son
descriptores adimensionales derivados de masa y vida media. El numero de
canales de decaimiento incluido en los datos es una variable exploratoria y no
un volumen riguroso del espacio de fases.

## Ejecucion

Requiere Python 3.10 o posterior y no instala dependencias:

```powershell
python -m src.atlas_persistence
python -m unittest discover -s tests -v
```

El primer comando regenera los CSV y el resumen en `output/`.

## Validacion predictiva exploratoria

La fase actual incluye una tabla de ablacion en `output/model_validation.csv`.
Compara modelos ridge mediante validacion leave-one-out sobre las identidades
inestables. La variable objetivo es `log10(tau)`.

Los predictores se incorporan progresivamente:

1. Intercepto.
2. Masa.
3. Masa y numero de modos de decaimiento representativos anotados.
4. Masa, modos y familia.
5. Masa, modos, familia e interaccion dominante.

`N_C = f_C tau` queda excluido de la prediccion porque contiene la propia vida
media. Su correlacion con `tau` es descriptiva y circular.

En la muestra actual, el numero de modos anotados no mejora el modelo de masa.
La interaccion dominante si reduce el error fuera de muestra, pero puede
codificar informacion cercana al mecanismo de decaimiento. Es una pista para
ampliar el estudio, no evidencia de una nueva ley de persistencia.

## Estructura

```text
data/       datos de entrada trazables y editables
src/        calculos y CLI
tests/      comprobaciones matematicas y fisicas basicas
output/     informes regenerables
```

## Datos y fuentes

Los valores de particulas son una seleccion exploratoria de valores centrales
publicados por Particle Data Group (PDG). Las constantes se basan en CODATA
2022/NIST. Cada CSV incluye una columna `source`.

- [PDG particle listings, 2025 update](https://pdg.lbl.gov/2025/listings/contents_listings.html)
- [PDG summary tables, 2025 update](https://pdg.lbl.gov/2025/tables/contents_tables.html)
- [NIST CODATA values](https://physics.nist.gov/cuu/Constants/)

Los escenarios de exergia y algunos objetos compactos son ejemplos de referencia
para validar formulas, no observaciones ajustadas.

`representative_decay_mode_count` es una anotacion curada de modos
representativos. No pretende contar todos los canales exclusivos ni sustituir
una integral de espacio de fases. `dominant_interaction` es una clasificacion
fisica amplia y tambien debe tratarse como variable exploratoria.

## Preguntas falsables

1. Tras controlar por masa, carga, espin y leyes de conservacion, comprobar si
   un descriptor de canales accesibles mejora la prediccion de vida media fuera
   de muestra.
2. Comparar huecos logaritmicos observados con catalogos simulados que respeten
   sesgos experimentales. Un hueco visual por si solo no demuestra una frecuencia
   privilegiada.
3. Evaluar si `delta_X` predice perdida de capacidad de reconstruccion en
   procesos termodinamicos concretos mejor que eficiencia energetica aislada.
4. Tratar `Phi` como coordenada de regimen y contrastar donde una aproximacion
   newtoniana deja de cumplir una tolerancia prefijada.
5. Sustituir el proxy de modos por observables definidos de forma uniforme:
   anchuras parciales, branching fractions, espacio de fases y reglas de
   seleccion.
