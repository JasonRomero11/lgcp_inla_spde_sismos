# Modelado espacio-temporal de eventos sísmicos en Colombia entre 2005 – 2020 mediante procesos de Cox Log-Gaussian y aproximación Bayesiana INLA-SPDE 

**Tesis de Maestría en Ciencias de la Información y Comunicaciones — Énfasis en Geomática**
Universidad Distrital Francisco José de Caldas
**Autor:** Jason Mauricio Romero Ríos
**Periodo de estudio:** 2005 – 2020

---

## Descripción del Proyecto

Este repositorio contiene los scripts de R desarrollados para el ajuste de un modelo de **Proceso Cox Log-Gaussiano (LGCP)** aplicado a la sismicidad continental de Colombia para el período 2005–2020, usando la metodología Bayesiana **INLA-SPDE**.

El trabajo aborda tres etapas:

1. **Análisis Exploratorio de Datos Espaciales (ESDA):** caracterización descriptiva y visual del catálogo sísmico.
2. **Estimación de distribuciones a priori informativas** mediante simulación de procesos LGCP y entrenamiento de una red neuronal convolucional (CNN-1D).
3. **Ajuste del modelo LGCP espacio-temporal** con INLA-SPDE, incorporando covariables geológicas y el campo espacial latente.

### Hipótesis central

> La ocurrencia de eventos sísmicos en Colombia puede modelarse como un Proceso Cox Log-Gaussiano espacio-temporal, cuya intensidad es explicada parcialmente por factores geológicos (fallas, volcanes, isostasia, topografía) y por un campo Gaussiano latente que captura la dependencia espacial residual.

---

## Estructura del Repositorio

```
RepositorioTesisMCIC/
│
├── scripts/                         # Scripts principales del análisis
│   ├── ESDA/
│   │   └── ESDA.R                   # Análisis exploratorio de datos espaciales
│   ├── INLA - SPDE/
│   │   ├── utils.R                  # Funciones auxiliares reutilizables
│   │   ├── spatial_model_INLA_SPDE_2020.R       # Modelo espacial (año 2020)
│   │   └── spatio_temporal_INLA_SPDE_2005_2020.R # Modelo espacio-temporal
│   ├── SIMULACIONES_LGCP/
│   │   └── simulations_rglcp.R      # Simulaciones LGCP para entrenamiento CNN
│   └── ANOMALIAS_GRAVIMETRICAS/
│       └── anomaliasGravimetricas.ipynb  # Cálculo de anomalías isostáticas (Python)
│
├── covariables_rds/                 # Covariables preprocesadas (objetos im de spatstat)
│   ├── script_rds.R                 # Script de preprocesamiento de covariables
│   ├── topografia_im_scaled.rds     # Elevación (MSNM) escalada a [-1,1]
│   ├── isostasia_im_scaled.rds      # Anomalía isostática escalada
│   ├── volcanes_im_scaled.rds       # Distancia a volcanes escalada
│   ├── inversa_im_scaled.rds        # Distancia a fallas inversas escalada
│   ├── normal_im_scaled.rds         # Distancia a fallas normales escalada
│   ├── dextral_im_scaled.rds        # Distancia a fallas dextrales escalada
│   ├── sinestral_im_scaled.rds      # Distancia a fallas sinestrales escalada
│   └── shapeZona_sp                 # Shapefile de la zona de estudio (continental)
│
├── R/                               # Scripts de referencia (INLA-SPDE book)
│   ├── spde-book-functions.R        # Funciones auxiliares del libro SPDE (Lindgren)
│   ├── discrete_gradient.R          # Gradiente discreto para visualización
│   └── ...                          # Otros ejemplos de referencia INLA
│
├── DiagramasArquitectura/           # Diagramas de la arquitectura CNN
│   ├── arquitectura.png             # Arquitectura CNN-1D propuesta
│   ├── arquitecturaReferencia.png   # Arquitectura CNN de referencia (VIHRS 2022)
│   ├── simulacionlgcp.png           # Flujo de simulación LGCP (versión base)
│   └── simulacionlgcp_improve.png   # Flujo de simulación LGCP (versión mejorada)
│
└── Imagenes/                        # Figuras organizadas por etapa del análisis
    ├── ESDA/                        # Mapas y gráficos exploratorios
    ├── Metodologia_SPDE/            # Malla triangular y teselación de Voronoi
    ├── CNN_Priors/                  # Entrenamiento CNN y estimación de priors
    ├── Modelos_Espaciales/          # Resultados del modelo espacial (2020)
    └── Modelo_Espacio_Temporal/     # Resultados del modelo espacio-temporal
```

---

## Metodología

### 1. Análisis Exploratorio (ESDA)

Caracterización del catálogo sísmico 2005–2020 del **Servicio Geológico Colombiano (SGC)**, proyectado en **EPSG:3116** (Magna Sirgas Colombia Bogotá). Se analizan distribuciones temporales, espaciales, magnitud y profundidad, junto con las covariables geológicas (fallas, topografía, isostasia, volcanes).

### 2. Simulación LGCP y Estimación de Priors mediante CNN

Para construir distribuciones a priori informativas sobre los hiperparámetros del campo Matérn (rango `r` y desviación estándar `σ`), se adoptó el enfoque de **Verönneau-Iphigenie et al. (2022)**, extendido con features de primer orden:

- Se generan **15,000 realizaciones** de procesos LGCP para entrenamiento y **1,000 para test**, sobre la ventana real de Colombia (EPSG:3116), bajo parámetros aleatorios: `μ ∈ [−19.5, −18.5]`, `σ² ∈ [1.0, 3.0]`, `scale ∈ [80 km, 200 km]`.
- Cada realización se simula con `rLGCP(model = "matern", nu = 1)` en una grilla de 128×128, filtrando patrones con entre 30 y 150,000 puntos.
- Para cada realización se extrae la **función L de Besag** `L̂(r) − r` (corrección de borde, `rmax = 200 km`, 128 valores) y **8 features de primer orden** que caracterizan la heterogeneidad espacial de la intensidad:
  - **Quadrat-based (3):** varianza de conteos, índice de dispersión (VMR), ratio max/min.
  - **Kernel density (5):** varianza, asimetría, curtosis, entropía normalizada y coeficiente de variación del campo suavizado (bandwidth = 50 km, grilla 64×64).
- La simulación se ejecuta en paralelo (`pbmcapply`, 10 cores) por chunks de 50 realizaciones, guardando resultados incrementales en archivos `.rds`.
- Una **red neuronal CNN-1D** (entrenada externamente en Python con Keras/TensorFlow) aprende la relación entre L̂(r) + features de primer orden y los parámetros generadores `(μ, σ², scale)`.
- La CNN se aplica al catálogo real para obtener estimaciones iniciales de los parámetros a priori.

**Parámetros estimados para Colombia:**

| Parámetro | Valor estimado |
|-----------|---------------|
| μ (intensidad media log) | −17.77 |
| σ² (varianza del campo latente) | 2.40 |
| scale / rango espacial | 110,342 m (~110 km) |

Estos valores alimentan los **PC-priors** (Penalised Complexity priors) del modelo INLA-SPDE.

### 3. Modelo LGCP Espacial (año 2020)

Se ajustan cuatro especificaciones del modelo LGCP para evaluar el efecto de los priors informativos y las covariables:

| Modelo | PC-priors | Covariables |
|--------|-----------|-------------|
| M0     | No        | No          |
| M1     | Sí        | No          |
| M2     | No        | Sí          |
| M3     | Sí        | Sí          |

**Covariables finales** (tras análisis de multicolinealidad – VIF < 10):
- Distancia a volcanes
- Distancia a fallas inversas
- Distancia a fallas normales
- Anomalía isostática
- Elevación (MSNM) *(no significativa, excluida en modelo final)*

### 4. Modelo LGCP Espacio-Temporal (2005–2020)

Extensión del modelo espacial incorporando una dinámica temporal **AR(1)** sobre el campo Gaussiano latente, con períodos bienales (8 grupos: 2005-2006, ..., 2019-2020). La estructura de covarianza espacio-temporal usa el producto de Kronecker:

```
Q = Q_T ⊗ Q_S
```

donde `Q_S` es la matriz de precisión espacial SPDE y `Q_T` es la matriz de precisión AR(1).

---

## Resultados Principales

### Efectos fijos del modelo M3 (espacial)

| Covariable      | Media   | IC 95%              |
|-----------------|---------|---------------------|
| Intercepto      | −16.264 | [−17.462, −15.067]  |
| Volcanes        | −3.561  | [−5.084, −2.039]    |
| Falla inversa   | −8.589  | [−10.999, −6.178]   |
| Falla normal    | −3.214  | [−5.324, −1.103]    |
| Isostasia       | −0.187  | [−0.299, −0.074]    |

> Los coeficientes negativos en las variables de distancia indican **mayor intensidad sísmica en zonas próximas** a estas estructuras geológicas (efecto de proximidad).

### Hiperparámetros del modelo espacio-temporal

| Parámetro              | Media     | IC 95%                       |
|------------------------|-----------|------------------------------|
| Rango espacial (ρ)     | 142,063 m | [130,087 m, 154,311 m]       |
| Desviación (σ)         | 2.204     | [2.025, 2.387]               |
| Coeficiente AR(1) (a)  | 0.943     | [0.937, 0.948]               |

> El coeficiente AR(1) de 0.943 indica **fuerte persistencia temporal**: los patrones de sismicidad evolucionan gradualmente entre períodos consecutivos.

### Validación

- El modelo **M3** supera consistentemente a M2 en Log-Score y LCPO.
- Bootstrap con B=10,000: `P(M3 > M2) = 100%`, IC 95% completamente positivo `[0.35, 0.78]`.
- Residuos de Pearson sin estructura espacial sistemática (ausencia de sesgo).

---

## Figuras Clave

| Figura | Descripción |
|--------|-------------|
| ![Malla](Imagenes/Metodologia_SPDE/Triangulacion.png) | Triangulación de Delaunay (K=4,309 vértices) |
| ![Voronoi](Imagenes/Metodologia_SPDE/TriangulacionVoronoi.png) | Teselación de Voronoi para pesos de integración |
| ![CNN](Imagenes/CNN_Priors/arquitectura.png) | Arquitectura CNN-1D propuesta |
| ![Spatial M0-M1](Imagenes/Modelos_Espaciales/spatial_effects_M0_M1.png) | Campo latente: M0 vs M1 |
| ![Spatial M2-M3](Imagenes/Modelos_Espaciales/spatial_effects_M2_M3.png) | Campo latente: M2 vs M3 |
| ![Log-Score](Imagenes/Modelos_Espaciales/log_score_lcpo.png) | Comparación Log-Score y LCPO |
| ![Bootstrap](Imagenes/Modelos_Espaciales/boostrap.png) | Distribución bootstrap M3 vs M2 |
| ![Intensidad](Imagenes/Modelos_Espaciales/result_intensidades.png) | Intensidad estimada λ(s) – Modelo M3 |
| ![Temporal](Imagenes/Modelo_Espacio_Temporal/espacio_temporal.png) | Campo latente espacio-temporal |

---

## Dependencias de Python

### Paquetes requeridos (anomalías gravimétricas)

```python
numpy       # Cálculo numérico
rasterio    # Lectura/escritura de rasters geoespaciales
```

Instalación:
```bash
pip install numpy rasterio
```

---

## Dependencias de R

### Paquetes requeridos

```r
# Inferencia Bayesiana
library(INLA)       # Integrated Nested Laplace Approximation (>= 24.11.25)
library(fmesher)    # Construcción de mallas FEM para INLA-SPDE
library(inlabru)    # Interfaz de alto nivel para modelos INLA

# Datos espaciales
library(sf)         # Simple Features para vectoriales
library(sp)         # Clases espaciales legacy (requerido por INLA)
library(terra)      # Manejo de rasters
library(raster)     # Rasters (compatibilidad legacy)

# Procesos puntuales
library(spatstat)          # Análisis de procesos puntuales
library(spatstat.geom)     # Geometría de patrones puntuales
library(spatstat.explore)  # Funciones de segundo orden (K, L, G, F)
library(spatstat.random)   # Simulación de procesos (rLGCP)

# Manipulación de datos
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

# Visualización
library(ggplot2)
library(viridis)
library(patchwork)
library(cowplot)
library(ggspatial)
library(tidyterra)

# Dependencias adicionales
library(spdep)
library(pbmcapply)   # Paralelización con barra de progreso
library(xtable)      # Tablas LaTeX
library(scales)
```

### Instalación de INLA

INLA no está disponible en CRAN. Instalar desde el repositorio oficial:

```r
install.packages("INLA",
  repos = c(INLA = "https://inla.r-inla-download.org/R/stable"),
  dependencies = TRUE
)
```

---

## Entorno Computacional

| Componente | Especificación |
|------------|---------------|
| Procesador | AMD Ryzen 9 7950X (16 núcleos / 32 hilos) |
| RAM | 64 GB DDR5 |
| GPU | NVIDIA GeForce RTX 4090 (24 GB) |
| SO | Ubuntu 22.04.5 LTS |
| R | Versión 4.5.2 |
| INLA | Versión 24.11.25 |

> Los modelos espaciales tardan ~50 segundos. El modelo espacio-temporal (~8 períodos × 4,309 vértices) tarda ~33 minutos sin priors y ~30 minutos con PC-priors.

---

## Orden de Ejecución

### Paso 0: Cálculo de anomalías gravimétricas isostáticas (Python)
```
scripts/ANOMALIAS_GRAVIMETRICAS/anomaliasGravimetricas.ipynb
```
*Calcula las anomalías isostáticas de gravedad (AIG) a partir de rasters DEM y gravedad observada. Incluye: reducción de Aire Libre, Placa de Bouguer, correcciones topográficas e isostáticas (modelo Airy-Heiskanen), y generación del mosaico final. Requiere: `numpy`, `rasterio`.*

### Paso 1: Preprocesamiento de covariables
```
covariables_rds/script_rds.R
```
*Requiere los datos originales (rasters, shapefiles). Genera los archivos `.rds` en `covariables_rds/`.*

### Paso 2: Análisis exploratorio
```
scripts/ESDA/ESDA.R
```
*Requiere: catálogo sísmico `.gpkg`, shapefiles de fallas, rasters de topografía e isostasia.*

### Paso 3: Simulaciones LGCP y entrenamiento CNN
```
scripts/SIMULACIONES_LGCP/simulations_rglcp.R
```
*Genera 15,000 realizaciones LGCP (train) + 1,000 (test) sobre la ventana continental de Colombia. Para cada realización extrae la función L̂(r) (128 valores) y 8 features de primer orden (quadrat-based + kernel density). La simulación se paraleliza en 10 cores por chunks de 50 y guarda `.rds` incrementales. La CNN se entrena externamente en Python (Keras/TensorFlow).*

### Paso 4: Modelo espacial (año 2020)
```
scripts/INLA - SPDE/spatial_model_INLA_SPDE_2020.R
```
*Fuente: `scripts/INLA - SPDE/utils.R`. Requiere covariables `.rds`.*

### Paso 5: Modelo espacio-temporal (2005–2020)
```
scripts/INLA - SPDE/spatio_temporal_INLA_SPDE_2005_2020.R
```
*Fuente: `scripts/INLA - SPDE/utils.R`. Requiere covariables `.rds`.*

---

## Datos de Entrada

| Dato | Fuente | Formato |
|------|--------|---------|
| Catálogo sísmico 2005–2020 | Servicio Geológico Colombiano (SGC) | `.gpkg` (EPSG:3116) |
| Límite continental Colombia | SGC / IGAC | `.geojson` (EPSG:3116) |
| Topografía (DEM 2000m) | SRTM / NASA | `.tif` (raster) |
| Anomalía isostática | SGC | `.tif` (raster) |
| Fallas geológicas (4 tipos) | SGC | `.shp` |
| Volcanes | SGC | `.shp` |

> **Nota:** Los datos de entrada **no están incluidos** en este repositorio por restricciones de tamaño y licencias. Los archivos `.rds` preprocesados sí están incluidos en `covariables_rds/`.

---

## Referencias

- Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *JRSS-B*, 73(4), 423–498.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations. *JRSS-B*, 71(2), 319–392.
- Fuglstad, G.A., Simpson, D., Lindgren, F., & Rue, H. (2019). Constructing priors that penalize the complexity of Gaussian random fields. *JASA*, 114(525), 445–452.
- Simpson, D., Rue, H., Riebler, A., Martins, T.K., & Sørbye, S.H. (2017). Penalising model component complexity: a principled practical approach to constructive priors. *Statistical Science*, 32(1), 1–28.
- Verönneau-Iphigenie, H.R.S. et al. (2022). Simulation-based inference for spatial point processes using summary statistics. *Spatial Statistics*, 100668.
- Cameletti, M., Lindgren, F., Simpson, D., & Rue, H. (2013). Spatio-temporal modeling of particulate matter concentration through the SPDE approach. *AStA*, 97(2), 109–131.
- Gómez-Rubio, V. (2020). *Bayesian inference with INLA*. CRC Press.

---

*Repositorio preparado como parte de la Tesis de Maestría, Universidad Distrital Francisco José de Caldas, 2025–2026.*
