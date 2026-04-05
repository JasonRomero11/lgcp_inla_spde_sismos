# =============================================================================
# SCRIPT: spatial_model_INLA_SPDE_2020.R
# =============================================================================
# PROPÓSITO:
#   Ajusta cuatro variantes del modelo LGCP (Log-Gaussian Cox Process) sobre
#   el catálogo de eventos sísmicos de Colombia para el año 2020, usando la
#   aproximación bayesiana INLA con el enfoque SPDE de Lindgren et al. (2011).
#
# DESCRIPCIÓN:
#   Se evalúan cuatro especificaciones del modelo para comparar el efecto de
#   incorporar priors informativos (estimados mediante CNN sobre simulaciones
#   LGCP) y covariables geológicas:
#
#     M0: sin priors informativos, sin covariables
#     M1: con priors informativos, sin covariables
#     M2: sin priors informativos, con covariables
#     M3: con priors informativos, con covariables  ← modelo final seleccionado
#
#   El script construye:
#     1. Un patrón puntual (ppp) y la malla de Delaunay (SPDE mesh)
#     2. La malla dual de Voronoi como pesos de integración
#     3. El stack INLA con datos de proceso puntual (y = 0/1, E = areas Voronoi)
#     4. Covariables evaluadas en vértices del mesh + ubicaciones de eventos
#     5. Los cuatro modelos mediante la función modelo_lgcp_col()
#     6. Métricas de comparación: Log-Score, LCPO, residuos de Pearson
#     7. Prueba Bootstrap (B=5000) para comparar M2 vs M3
#
# FUENTE (source):
#   - scripts/INLA - SPDE/utils.R   (funciones de graficación y utilidades)
#   - R/spde-book-functions.R        (book.mesh.dual, funciones del libro SPDE)
#   - R/discrete_gradient.R          (helper para gradientes discretos)
#
# INPUTS:
#   - covariables_rds/shapeZona_sp               : polígono continental Colombia (EPSG:3116)
#   - Data/EventosColPointsPlanas31162005_2020_continental.gpkg : catálogo sísmico
#   - covariables_rds/*_im_scaled.rds            : covariables como objetos im de spatstat
#
# OUTPUTS:
#   - pp.resM0, pp.resM1, pp.resM2, pp.resM3  : listas con resultado INLA + grid + stack
#   - spdf_sf                                  : SpatialPixelsDataFrame con predicciones e intensidades
#   - Tablas de efectos fijos en LaTeX (via results_model_to_tableLatex)
#   - Gráficos: campo latente, intensidad predicha, Log-Score, LCPO, bootstrap
#
# PARÁMETROS A PRIORI (estimados por CNN sobre 10,000 simulaciones LGCP):
#   mu_simulated    = -17.76713   (intercepto log-lineal esperado)
#   scale_simulated = 2.400884    (desviación estándar del campo Matérn)
#   range_simulated = 110342.3 m  (rango de correlación espacial, ~110 km)
#
# AUTOR: Jason Mauricio Romero Ríos
# UNIVERSIDAD: Universidad Distrital Francisco José de Caldas
# TESIS: Maestría en Ciencias de la Información y Comunicaciones – Geomática
# =============================================================================


# =============================================================================
# SECCIÓN 1: CARGA DE PAQUETES
# =============================================================================

library(INLA)        # Inferencia bayesiana INLA
library(sp)          # Clases espaciales legacy (requeridas por INLA)
library(ggplot2)     # Visualización
library(sf)          # Simple Features: lectura y manejo de vectoriales
library(spatial)     # Utilidades espaciales básicas
library(spData)      # Datos de ejemplo para spatial
library(spdep)       # Dependencia espacial
library(maps)        # Mapas base
library(gridExtra)   # Composición de gráficos en grilla
library(spatstat)    # Análisis de procesos puntuales (ppp, owin, etc.)
library(deldir)      # Teselación de Voronoi (usado indirectamente por book.mesh.dual)
library(raster)      # Manejo de rasters (legacy, para compatibilidad con sp)
library(viridis)     # Paletas de color
library(terra)       # Rasters modernos (SpatRaster)
library(lubridate)   # Manejo de fechas
library(patchwork)   # Composición de plots ggplot2
library(fmesher)     # Construcción de mallas FEM para INLA-SPDE
library(tidyr)       # Transformación de datos

# Instala fmesher si no está disponible
if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable")
}


# =============================================================================
# SECCIÓN 2: CONFIGURACIÓN DE RUTAS Y PARÁMETROS GLOBALES
# =============================================================================

setwd("~/Documents/Personal/lgcp_inla_spde_sismos")

# Directorio de covariables preprocesadas (.rds)
files_rds <- "covariables_rds"

# Directorio de salida de imágenes
path_image_results <- 'imagenes_doc'

# Cargar funciones del libro SPDE (book.mesh.dual, etc.) y utilitarias
source("R/spde-book-functions.R")
source("R/discrete_gradient.R")
source("R/utils.R")   # → create_ppp(), create_mesh(), plot_spatial_effects_*(), etc.

# Ruta al catálogo sísmico (GeoPackage, EPSG:3116)
path_file_seismic <- "Data/EventosColPointsPlanas31162005_2020_continental.gpkg"


# =============================================================================
# SECCIÓN 3: PARÁMETROS A PRIORI (estimados mediante CNN – ver simulations_rglcp.R)
# =============================================================================

# Intercepto: media del proceso log-Gaussiano (en log-escala)
# Corresponde a log(intensidad media) ≈ exp(-17.77) ≈ 2.1e-8 eventos/m²
mu_simulated <- -20.6119

# Desviación estándar del campo latente Matérn
# Prior PC: P(sigma > scale_simulated) = 0.5
scale_simulated <- 3.8939

# Rango espacial de correlación (metros, ~110 km)
# Prior PC: P(rango < range_simulated) = 0.5
range_simulated <- 59900

# Banderas de control del script
save_graphics         <- F     # Si TRUE, guarda los gráficos en disco
use_aprior_information <- T    # Si TRUE, usa PC-priors informativos
use_pts_pred_wAreas   <- T     # Si TRUE, usa cuadrícula de predicción con áreas
size_grid_pred        <- 10000 # Tamaño de celda para la cuadrícula de predicción (m)

# Lista con prior para el intercepto (usado en control.fixed de INLA)
# mean.intercept: media a priori de mu
# prec.intercept: precisión a priori (= 1/varianza), 0.1 → varianza = 10 (poco informativo)
control_fixed_list <- list(
  mean.intercept = mu_simulated,
  prec.intercept = 0.1
)


# =============================================================================
# SECCIÓN 4: CARGA Y FILTRADO DEL CATÁLOGO SÍSMICO
# =============================================================================

# Cargar el límite del área de estudio (Colombia continental, EPSG:3116)
shapeZona_sp <- readRDS(paste0(files_rds, "/shapeZona_sp"))

# Simplificar la geometría para acelerar intersecciones (tolerancia = 5 km)
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = TRUE)

# Leer el catálogo sísmico completo (2005-2020)
sismosSp <- st_read(path_file_seismic)
sismosSp$year <- sismosSp$YEAR

# Filtrar solo el año 2020 (modelo espacial puro, sin estructura temporal)
sismosSp <- subset(sismosSp, year >= 2020)
sismosSp <- subset(sismosSp, year <= 2020)

# Extraer coordenadas X e Y para uso posterior
sismosSp$X <- st_coordinates(sismosSp)[, 1]
sismosSp$Y <- st_coordinates(sismosSp)[, 2]

# Verificar integridad de coordenadas
if (any(is.na(sismosSp$X) | is.na(sismosSp$Y))) {
  warning("Hay valores NA en las coordenadas. Eliminando...")
  sismosSp <- sismosSp[!is.na(sismosSp$X) & !is.na(sismosSp$Y), ]
}


# =============================================================================
# SECCIÓN 5: CREACIÓN DEL PATRÓN PUNTUAL (ppp) Y MALLA SPDE
# =============================================================================

# Crear el patrón puntual de procesos puntuales (spatstat::ppp) dentro de la
# ventana definida por shapeZona_sp. Devuelve:
#   $point_pattern      : objeto ppp con los eventos sísmicos
#   $coordinates_eventos: matriz nx2 con coordenadas de los eventos
creating_ppp <- create_ppp(shapeZona_sp, sismosSp)

p  <- creating_ppp$point_pattern       # Patrón puntual (ppp)
xy <- creating_ppp$coordinates_eventos # Coordenadas de los eventos

# ----- Construcción de la malla de Delaunay (SPDE mesh) -----
param_cutoff  <- 5000           # Distancia mínima entre vértices (m)
offset_param  <- c(100, 20000)  # Extensión interior y exterior de la malla (m)

# Asegurar geometría válida del polígono de estudio
shapeZona <- st_make_valid(shapeZona_sp)
shapeZona <- as_Spatial(shapeZona)
shapeZona <- st_as_sf(shapeZona)

# Calcular el tamaño máximo de arista en función de la extensión espacial.
# División por 40 produce una malla más fina que por 30 (más vértices, mayor precisión).
max.edge_params <- max(c(diff(range(xy[, 1])), diff(range(xy[, 2])))) / 40

# Construcción de la malla 2D de Delaunay para INLA-SPDE
# max.edge = c(interior, exterior): longitud máxima de aristas dentro y fuera de la región
# El resultado (meshSismos) tiene ~4219 vértices
meshSismos <- inla.mesh.2d(
  boundary = shapeZona,
  loc      = cbind(xy),
  max.edge = c(1, 3) * max.edge_params,
  cutoff   = param_cutoff,
  offset   = offset_param,
  crs      = st_crs(shapeZona)
)

cat("Número de vértices del mesh:", meshSismos$n)  # ~4219
plot(meshSismos, main = "Malla de Delaunay – Colombia 2020")


# =============================================================================
# SECCIÓN 6: CONSTRUCCIÓN DEL STACK INLA (datos del proceso puntual)
# =============================================================================

# --- Puntos del mesh y de los eventos ---
ptsSismos     <- as.matrix(xy)                     # Coordenadas de los eventos (n × 2)
mesh.ptsSismos <- as.matrix(meshSismos$loc[, 1:2]) # Vértices del mesh (nv × 2)
allptsSismos  <- rbind(mesh.ptsSismos, ptsSismos)  # Todos los puntos (nv + n) × 2

nvSismos <- meshSismos$n    # Número de vértices del mesh (nv)
nSismos  <- nrow(ptsSismos) # Número de eventos (n)

# --- Malla dual de Voronoi para los pesos de integración ---
# book.mesh.dual() construye los polígonos de Voronoi duales a cada vértice del mesh.
# Estos polígonos se intersectan con shapeZona_sp para obtener las áreas de integración.
dmesh <- st_as_sf(book.mesh.dual(meshSismos))
st_crs(dmesh) <- st_crs(shapeZona_sp)

# Identificar los polígonos de Voronoi que intersectan el área de estudio
intersect_idx <- st_intersects(dmesh, shapeZona_sp, sparse = FALSE)[, 1]
dmesh_in      <- dmesh[intersect_idx, ]
intersections <- st_intersection(dmesh_in, shapeZona_sp)

# Calcular el área de cada polígono de Voronoi dentro del estudio
areas <- st_area(intersections)

# Vector de pesos (áreas) de integración por vértice
# w[i] = área del polígono de Voronoi del vértice i dentro de shapeZona_sp
# w[i] = 0 para vértices fuera del dominio de estudio
w <- numeric(nrow(dmesh))
w[intersect_idx] <- as.numeric(areas)

# Verificación: la suma de pesos debe aproximar el área total de la región
total_area  <- as.numeric(st_area(shapeZona_sp))
sum_wSismos <- sum(w)
cat("Área total de Colombia continental:", total_area, "m²\n")
cat("Suma de áreas de Voronoi:", sum_wSismos, "m²\n")

# --- Variable respuesta y exposición del proceso puntual ---
# La formulación LGCP en INLA usa una Poisson con:
#   y = 0 (en vértices del mesh) y y = 1 (en ubicaciones de eventos)
#   E = área de Voronoi (en vértices) y E = 0 (en ubicaciones de eventos)
# Esto aproxima la log-verosimilitud del proceso puntual vía integración numérica.
y.ppSismos       <- rep(0:1, c(nvSismos, nSismos))
wSismos_numeric  <- as.numeric(unlist(w))
e.ppSismos       <- as.numeric(c(w, rep(0, nSismos)))

# --- Matrices de proyección A ---
# lmatSismos: proyecta los eventos al mesh (matriz de pesos bilineales)
lmatSismos <- inla.spde.make.A(meshSismos, ptsSismos)
# imaSismos: identidad para los vértices del mesh
imaSismos  <- Diagonal(nvSismos, rep(1, nvSismos))
# A.ppSismos: matriz de proyección combinada (nv + n) × nv
A.ppSismos <- rbind(imaSismos, lmatSismos)


# =============================================================================
# SECCIÓN 7: EXTRACCIÓN DE COVARIABLES EN PUNTOS DEL MESH Y EVENTOS
# =============================================================================

# Crear patrón puntual expandido con los vértices del mesh + los eventos
expanded_window  <- grow.rectangle(as.rectangle(p$window), 50000)
allpts.pppSismos <- ppp(allptsSismos[, 1], allptsSismos[, 2], expanded_window)

# Cargar covariables preprocesadas como objetos im (image) de spatstat.
# Cada objeto im es una grilla raster escalada a [0,1] con escala_im_to_unit_range()
topografia_im_scaled     <- readRDS(paste0(files_rds, "/topografia_im_scaled.rds"))
isostasia_im_scaled      <- readRDS(paste0(files_rds, "/isostasia_im_scaled.rds"))
volcanes_im_scaled       <- readRDS(paste0(files_rds, "/volcanes_im_scaled.rds"))
falla_sinestral_im_scaled <- readRDS(paste0(files_rds, "/sinestral_im_scaled.rds"))
falla_dextral_im_scaled   <- readRDS(paste0(files_rds, "/dextral_im_scaled.rds"))
falla_normal_im_scaled    <- readRDS(paste0(files_rds, "/normal_im_scaled.rds"))
falla_inversa_im_scaled   <- readRDS(paste0(files_rds, "/inversa_im_scaled.rds"))

# Lista de covariables para el modelo
covar <- list(
  msnm           = topografia_im_scaled,      # Elevación (m s.n.m.), escalada
  isostasia      = isostasia_im_scaled,        # Anomalía isostática (mGal), escalada
  volcanes       = volcanes_im_scaled,         # Distancia a volcanes, escalada
  falla_sinestral = falla_sinestral_im_scaled, # Distancia a fallas sinestrales, escalada
  falla_dextral  = falla_dextral_im_scaled,   # Distancia a fallas dextrales, escalada
  falla_inversa  = falla_inversa_im_scaled,   # Distancia a fallas inversas, escalada
  falla_normal   = falla_normal_im_scaled     # Distancia a fallas normales, escalada
)

# Extraer el valor de cada covariable en todos los puntos (vértices del mesh + eventos).
# nearest.pixel() encuentra la celda de la imagen más cercana a cada punto.
covs100 <- lapply(covar, function(X) {
  pixels <- nearest.pixel(allpts.pppSismos$x, allpts.pppSismos$y, X)
  sapply(1:npoints(allpts.pppSismos), function(i) {
    X[pixels$row[i], pixels$col[i]]
  })
})

# Añadir intercepto manual (columna de unos)
covs100$b0 <- rep(1, nvSismos + nSismos)

# Verificar consistencia de longitud de covariables
stopifnot(
  length(covs100$isostasia) == length(covs100$volcanes),
  length(covs100$volcanes)  == length(covs100$b0)
)


# =============================================================================
# SECCIÓN 8: FUNCIÓN PRINCIPAL – modelo_lgcp_col()
# =============================================================================

#' Ajustar modelo LGCP espacial con INLA-SPDE
#'
#' Construye el stack INLA completo (estimación + predicción) y ajusta
#' un modelo de proceso de Cox Log-Gaussiano (LGCP) usando la aproximación
#' SPDE de Lindgren et al. (2011).
#'
#' @param interval1 numeric. Tamaño de celda (m) de la cuadrícula de predicción.
#'   Define la resolución de la estimación de intensidad λ(s).
#' @param use_covariables_model character vector. Nombres de las covariables a
#'   incluir en el predictor lineal. Subconjunto de names(covs100).
#'   Usar c() para modelo sin covariables.
#' @param use_aprior_information logical. Si TRUE, usa PC-priors informativos
#'   (inla.spde2.pcmatern) con range_simulated y scale_simulated.
#'   Si FALSE, usa la SPDE de Matérn estándar sin priors informativos.
#'
#' @return Lista con tres elementos:
#'   $result_inla : objeto inla con el resultado del ajuste
#'   $grid        : objeto sf con los puntos de predicción (cuadrícula)
#'   $full_stack  : stack INLA completo (estimación + predicción combinado)
#'
#' @details
#'   El modelo se formula como:
#'     y_i | λ(s_i) ~ Poisson(λ(s_i) × E_i)
#'     log(λ(s)) = β₀ + β₁X₁(s) + … + u(s)
#'     u(s) ~ GF_Matérn(rango, σ²)  ← campo latente SPDE
#'
#'   Donde y_i = 0 en vértices del mesh (con E_i = área de Voronoi)
#'         y_i = 1 en ubicaciones de eventos (con E_i = 0)
#'
#'   int.strategy = "eb" (Empirical Bayes) para reducir tiempo de cómputo.
#'
modelo_lgcp_col <- function(interval1, use_covariables_model, use_aprior_information) {

  # --- Especificación del campo latente Matérn (SPDE) ---
  if (use_aprior_information) {
    cat("-------------------------------------\n")
    cat("Modelo con PC-priors informativos\n")
    # PC-prior para el rango: P(rango < range_simulated) = 0.5
    # PC-prior para sigma:    P(sigma > scale_simulated) = 0.5
    spdesismos <- inla.spde2.pcmatern(
      mesh        = meshSismos,
      prior.range = c(range_simulated, 0.5),
      prior.sigma = c(scale_simulated, 0.5)
    )
  } else {
    # SPDE de Matérn estándar: sin priors informativos
    # alpha = 2 → ν = 1 (suavidad del campo), proceso de Whittle-Matérn
    spdesismos <- inla.spde2.matern(meshSismos, alpha = 2)
  }
  n_spde <- spdesismos$n.spde

  # --- Índice SPDE para identificar el efecto aleatorio espacial ---
  spde.indexSismos <- inla.spde.make.index(
    name  = "spatial.field",
    n.spde = n_spde
  )

  # --- Stack de estimación: combina vértices del mesh + eventos sísmicos ---
  # data$y = vector de respuesta (0 en vértices, 1 en eventos)
  # data$e = vector de exposición (área de Voronoi en vértices, 0 en eventos)
  # A[[1]] = A.ppSismos: proyecta el campo SPDE
  # A[[2]] = 1: proyecta los efectos fijos (covariables)
  spde.stack <- inla.stack(
    data    = list(y = y.ppSismos, e = e.ppSismos),
    A       = list(A.ppSismos, 1),
    effects = list(spde.indexSismos, covs100),
    tag     = "pp"
  )

  # --- Cuadrícula de predicción (puntos interiores al área de estudio) ---
  # Se usa un buffer negativo para evitar predicciones en el borde
  shapeZona_sp_interno <- st_buffer(shapeZona_sp, -1000)
  grid     <- st_make_grid(shapeZona_sp_interno, cellsize = interval1, what = "centers", square = TRUE)
  grid_in  <- grid[st_within(grid, shapeZona_sp_interno, sparse = FALSE)]
  pts.pred <- as.data.frame(st_coordinates(grid_in))

  # Patrón puntual de predicción para extracción de covariables
  ppp.pred <- ppp(pts.pred[, 1], pts.pred[, 2], window = p$window)

  # Matriz de proyección para los puntos de predicción
  A.pred <- inla.spde.make.A(mesh = meshSismos, loc = as.matrix(pts.pred))

  # Covariables evaluadas en los puntos de predicción (interpolación bilineal)
  covs100.pred <- lapply(covar, function(X) { X[ppp.pred] })
  covs100.pred$b0 <- rep(1, nrow(pts.pred))

  # --- Stack de predicción: y = NA indica que INLA debe predecir ---
  spde.stack.pred <- inla.stack(
    data    = list(y = NA),
    A       = list(A.pred, 1),
    effects = list(spde.indexSismos, covs100.pred),
    tag     = "pred"
  )

  # --- Stack combinado: estimación + predicción ---
  join.stack <- inla.stack(spde.stack, spde.stack.pred)

  # --- Fórmula del modelo ---
  # Intercepto + covariables seleccionadas + efecto aleatorio espacial SPDE
  formula <- as.formula(
    paste("y ~ 1 +",
          paste(use_covariables_model, collapse = " + "),
          "+ f(spatial.field, model = spdesismos)")
  )

  # --- Argumentos para la llamada a INLA ---
  args_inla <- list(
    formula  = formula,
    family   = "poisson",           # Proceso de Poisson (aproxima LGCP)
    data     = inla.stack.data(join.stack),
    control.inla      = list(int.strategy = "eb"),    # Empirical Bayes para rapidez
    control.predictor = list(
      A       = inla.stack.A(join.stack),
      compute = TRUE,
      link    = 1                   # log-link (Poisson)
    ),
    control.compute = list(
      config = TRUE,   # Guarda configuración para simulaciones posteriores
      cpo    = TRUE,   # Calcula el CPO (Conditional Predictive Ordinate) para LCPO
      waic   = TRUE,   # Calcula WAIC
      dic    = TRUE    # Calcula DIC
    ),
    E = inla.stack.data(join.stack)$e  # Exposición (áreas de Voronoi)
  )

  # Añadir prior en efectos fijos solo si se usan priors informativos
  if (use_aprior_information) args_inla$control.fixed <- control_fixed_list

  # --- Ajuste del modelo ---
  pp.res <- do.call(inla, args_inla)

  return(list(
    result_inla = pp.res,    # Objeto inla completo
    grid        = grid_in,   # sf con los puntos de predicción
    full_stack  = join.stack # Stack combinado (estimación + predicción)
  ))
}


# =============================================================================
# SECCIÓN 9: AJUSTE DE LOS CUATRO MODELOS (M0 – M3)
# =============================================================================

# Tamaño de celda de predicción: raíz cuadrada del área media de Voronoi
# Se redefine explícitamente a 15 km para mayor resolución
interval1 <- sqrt(mean(w))
print(interval1)
interval1 <- 15000  # 15 km (resolución de la cuadrícula de predicción)

# M0: sin prior informativo, sin covariables (modelo nulo de referencia)
time_M0 <- system.time({
  pp.resM0 <- modelo_lgcp_col(interval1, c(), FALSE)
})

# M1: con prior informativo, sin covariables (efecto del prior sobre el campo latente)
time_M1 <- system.time({
  pp.resM1 <- modelo_lgcp_col(interval1, c(), TRUE)
})

# M2: sin prior informativo, con covariables geológicas
time_M2 <- system.time({
  pp.resM2 <- modelo_lgcp_col(
    interval1,
    c("volcanes", "falla_inversa", "falla_normal", "isostasia"),
    FALSE
  )
})

# M3: con prior informativo + covariables geológicas (modelo final seleccionado)
time_M3 <- system.time({
  pp.resM3 <- modelo_lgcp_col(
    interval1,
    c("volcanes", "falla_inversa", "falla_normal", "isostasia"),
    TRUE
  )
})

# Tabla resumen de tiempos de ajuste
times <- rbind(M0 = time_M0, M1 = time_M1, M2 = time_M2, M3 = time_M3)
print(times)

# Tablas de efectos fijos en formato LaTeX (via utils.R)
results_model_to_tableLatex(pp.resM0$result_inla)
results_model_to_tableLatex(pp.resM1$result_inla)
results_model_to_tableLatex(pp.resM2$result_inla)
results_model_to_tableLatex(pp.resM3$result_inla)


# =============================================================================
# SECCIÓN 10: VISUALIZACIÓN DEL CAMPO LATENTE ESPACIAL
# =============================================================================

if (save_graphics) {
  # Panel 2×2: campo latente para M0, M1, M2, M3
  out_spatialeffect <- plot_spatial_effects_2x2(
    resM0 = pp.resM0$result_inla,
    resM1 = pp.resM1$result_inla,
    resM2 = pp.resM2$result_inla,
    resM3 = pp.resM3$result_inla,
    output_path = file.path(path_image_results, "modelos_INLA", "spatial_effects_M0_M3_2x2.png")
  )

  # Panel 1×2: comparación M2 vs M3 con estadísticos superpuestos
  out_M1_M2 <- plot_spatial_effects_1x2(
    resA          = pp.resM2$result_inla,
    resB          = pp.resM3$result_inla,
    titles        = c("", ""),
    output_path   = file.path(path_image_results, "modelos_INLA", "spatial_effects_M2_M3.png"),
    show_stats    = TRUE,
    stats_position    = "topright",
    use_quantile_limits = TRUE,
    q             = c(0.02, 0.98)
  )
}


# =============================================================================
# SECCIÓN 11: EXTRACCIÓN DE INTENSIDADES PREDICHAS (λ̂(s))
# =============================================================================

# Usar el grid de M1 (todos los modelos comparten la misma cuadrícula de predicción)
grid_in    <- pp.resM1$grid
join.stack <- pp.resM1$full_stack

# Convertir los puntos de predicción a SpatialPixelsDataFrame para visualización
pts.pred_ <- as.data.frame(st_coordinates(grid_in))
names(pts.pred_) <- c("x", "y")
pts.pred_$dummy <- 1
coordinates(pts.pred_) <- ~x + y
proj4string(pts.pred_) <- st_crs(shapeZona_sp)$proj4string
gridded(pts.pred_)     <- TRUE

# Convertir a sf para operaciones espaciales
spdf       <- as(pts.pred_, "SpatialPixelsDataFrame")
spdf_poly  <- as(spdf, "SpatialPolygonsDataFrame")
spdf_sf    <- st_as_sf(spdf_poly)
spdf_sf    <- st_transform(spdf_sf, crs = st_crs(sismosSp))

# Índice de los datos de predicción en el stack combinado
idx <- inla.stack.index(join.stack, "pred")$data

# Extraer la media posterior del predictor lineal (en escala log) y exponenciar
# para obtener la intensidad λ̂(s) = exp(η̂(s))
spdf$M0 <- exp(pp.resM0$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M1 <- exp(pp.resM1$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M2 <- exp(pp.resM2$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M3 <- exp(pp.resM3$result_inla$summary.linear.predictor[idx, "mean"])

# Transferir intensidades al objeto sf para análisis espacial
spdf_sf$M0 <- spdf$M0
spdf_sf$M1 <- spdf$M1
spdf_sf$M2 <- spdf$M2
spdf_sf$M3 <- spdf$M3

# Verificar consistencia de CRS entre predicciones y sismos
stopifnot(st_crs(spdf_sf) == st_crs(sismosSp))

# Contar eventos sísmicos observados en cada celda de predicción
spdf_sf$ID       <- seq_len(nrow(spdf_sf))
spdf_sf$observed <- lengths(st_intersects(spdf_sf, sismosSp))
spdf_sf$area     <- as.numeric(st_area(spdf_sf))

# Número esperado de eventos por celda = intensidad × área
spdf_sf$expect_M0 <- spdf_sf$M0 * spdf_sf$area
spdf_sf$expect_M1 <- spdf_sf$M1 * spdf_sf$area
spdf_sf$expect_M2 <- spdf_sf$M2 * spdf_sf$area
spdf_sf$expect_M3 <- spdf_sf$M3 * spdf_sf$area


# =============================================================================
# SECCIÓN 12: RESIDUOS DE PEARSON
# =============================================================================

#' Calcular residuos de Pearson para procesos de Poisson
#'
#' r_i = (O_i - E_i) / sqrt(E_i)
#'
#' @param O numeric. Vector de conteos observados.
#' @param E numeric. Vector de conteos esperados (debe ser > 0).
#' @param eps numeric. Umbral mínimo de E para evitar divisiones por cero (default 1e-12).
#'
#' @return numeric. Vector de residuos de Pearson (NA donde E ≤ 0 o valores no finitos).
pearson_resid <- function(O, E, eps = 1e-12) {
  ok  <- is.finite(O) & is.finite(E) & E > 0
  r   <- rep(NA_real_, length(O))
  E2  <- pmax(E, eps)
  r[ok] <- (O[ok] - E2[ok]) / sqrt(E2[ok])
  r
}

# Calcular residuos para cada modelo
spdf_sf$pearson_M0 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M0)
spdf_sf$pearson_M1 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M1)
spdf_sf$pearson_M2 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M2)
spdf_sf$pearson_M3 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M3)

# Resumen estadístico de los residuos (no debería haber sesgo sistemático)
summary(spdf_sf$pearson_M0)
summary(spdf_sf$pearson_M1)
summary(spdf_sf$pearson_M2)
summary(spdf_sf$pearson_M3)

# Transferir residuos al SpatialPixelsDataFrame para visualización con spplot
spdf$observed   <- spdf_sf$observed
spdf$area       <- spdf_sf$area
spdf$expect_M0  <- spdf_sf$expect_M0
spdf$expect_M1  <- spdf_sf$expect_M1
spdf$expect_M2  <- spdf_sf$expect_M2
spdf$expect_M3  <- spdf_sf$expect_M3
spdf$pearson_M0 <- spdf_sf$pearson_M0
spdf$pearson_M1 <- spdf_sf$pearson_M1
spdf$pearson_M2 <- spdf_sf$pearson_M2
spdf$pearson_M3 <- spdf_sf$pearson_M3


# =============================================================================
# SECCIÓN 13: LOG-SCORE DE POISSON (métrica de ajuste)
# =============================================================================

#' Calcular el Log-Score de Poisson
#'
#' Suma del logaritmo de la función de masa de probabilidad Poisson:
#'   LS = Σᵢ log P(Oᵢ | λ = Eᵢ)
#' Un Log-Score menos negativo (mayor) indica mejor ajuste.
#'
#' @param observed integer. Conteos observados.
#' @param expected numeric. Conteos esperados (λ_i).
#' @param eps numeric. Umbral mínimo para expected (default 1e-12).
#'
#' @return numeric. Escalar: suma del log-score de Poisson.
poisson_log_score <- function(observed, expected, eps = 1e-12) {
  ok <- is.finite(observed) & is.finite(expected) & expected > 0
  O  <- observed[ok]
  E  <- pmax(expected[ok], eps)
  sum(dpois(O, lambda = E, log = TRUE))
}

# Calcular Log-Score para cada modelo
LS_M0 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M0)
LS_M1 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M1)
LS_M2 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M2)
LS_M3 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M3)

cat("Log-Scores (mayor es mejor):\n")
print(c(logscore_M0 = LS_M0, logscore_M1 = LS_M1, logscore_M2 = LS_M2, logscore_M3 = LS_M3))


# =============================================================================
# SECCIÓN 14: LCPO – LOGARITMO DEL CONDITIONAL PREDICTIVE ORDINATE
# =============================================================================

#' Calcular el LCPO (negativo, menor es mejor)
#'
#' LCPO = -Σ log(CPO_i)
#' El CPO_i es la densidad predictiva marginal leave-one-out de la observación i.
#' INLA calcula el CPO directamente cuando control.compute$cpo = TRUE.
#' Un LCPO menor indica mejor poder predictivo del modelo.
#'
#' @param result_inla Objeto inla. Resultado del ajuste de INLA.
#'
#' @return numeric. Escalar: LCPO (NA si CPO no está disponible o es inválido).
lcpo_score <- function(result_inla) {
  cpo <- result_inla$cpo$cpo
  if (is.null(cpo) || length(cpo) == 0) return(NA_real_)
  ok <- is.finite(cpo) & cpo > 0
  -sum(log(cpo[ok]))
}

# Calcular LCPO para cada modelo
LCPO_M0 <- lcpo_score(pp.resM0$result_inla)
LCPO_M1 <- lcpo_score(pp.resM1$result_inla)
LCPO_M2 <- lcpo_score(pp.resM2$result_inla)
LCPO_M3 <- lcpo_score(pp.resM3$result_inla)

cat("LCPO (menor es mejor):\n")
print(c(resLCPO_M0 = LCPO_M0, resLCPO_M1 = LCPO_M1, resLCPO_M2 = LCPO_M2, resLCPO_M3 = LCPO_M3))

# Gráfico comparativo de Log-Score y LCPO para los cuatro modelos (via utils.R)
p1 <- grafico_comparacion_metricas(LS_M0, LS_M1, LS_M2, LS_M3,
                                   LCPO_M0, LCPO_M1, LCPO_M2, LCPO_M3)


# =============================================================================
# SECCIÓN 15: PRUEBA BOOTSTRAP – COMPARACIÓN M2 vs M3
# =============================================================================
# Procedimiento: bootstrap no paramétrico con B = 5000 remuestras.
# En cada remuestra b:
#   1. Seleccionar n celdas con reemplazo de las n celdas válidas
#   2. Calcular el Log-Score de Poisson para M2 y M3 sobre esa remuestra
#   3. δ[b] = LS_M3[b] - LS_M2[b]
# Si δ > 0 → M3 es mejor que M2 en esa remuestra.
# Resultado: P(M3 > M2) = proporción de remuestras donde δ > 0.

B <- 5000
delta <- numeric(B)

# Filtrar celdas válidas para el bootstrap (ambas esperanzas deben ser positivas y finitas)
ok <- is.finite(spdf_sf$observed) & is.finite(spdf_sf$expect_M2) &
  is.finite(spdf_sf$expect_M3) &
  spdf_sf$expect_M2 > 0 & spdf_sf$expect_M3 > 0

O  <- spdf_sf$observed[ok]
E1 <- spdf_sf$expect_M2[ok]   # Esperanza bajo M2
E2 <- spdf_sf$expect_M3[ok]   # Esperanza bajo M3
n  <- length(O)

for (b in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  ls1     <- sum(dpois(O[idx], lambda = E1[idx], log = TRUE))  # LS de M2 en remuestra b
  ls2     <- sum(dpois(O[idx], lambda = E2[idx], log = TRUE))  # LS de M3 en remuestra b
  delta[b] <- ls2 - ls1  # δ > 0 significa M3 mejor que M2
}

# Resumen del bootstrap
cat("Bootstrap (M3 vs M2):\n")
print(c(
  mean_delta   = mean(delta),
  p_M3_better  = mean(delta > 0),       # Proporción de veces que M3 supera a M2
  q025         = quantile(delta, 0.025),
  q975         = quantile(delta, 0.975)
))

# Guardar el vector delta para uso en grafico_bootstrap()
delta_bootstrap <- delta

# Gráfico de la distribución bootstrap de δ (via utils.R)
p3 <- grafico_bootstrap(delta_bootstrap)


# =============================================================================
# SECCIÓN 16: VISUALIZACIÓN DE RESIDUOS DE PEARSON
# =============================================================================

# Panel 2×2 de residuos de Pearson para los cuatro modelos
spplot(
  spdf,
  c("pearson_M0", "pearson_M1", "pearson_M2", "pearson_M3"),
  names.attr  = c("M0", "M1", "M2", "M3"),
  col.regions = viridis::plasma(32),
  main        = "Residuos de Pearson por modelo"
)


# =============================================================================
# SECCIÓN 17: MAPA DE INTENSIDAD ESTIMADA λ̂(s) – MODELO M3
# =============================================================================

# Visualizar la intensidad predicha λ̂(s) del modelo M3 (el mejor modelo)
# La función plot_intensidades_modelos() está definida en utils.R
resultado_int <- plot_intensidades_modelos(
  spdf_sf            = spdf_sf,
  cols_intensidad    = c("M3"),      # Columna de intensidad a graficar
  modelos            = c("M3"),      # Etiqueta del modelo
  shapeZona          = shapeZona_sp, # Polígono de la región de estudio
  output_path        = "intensidades_modelos.png",
  ncol               = 1,
  show_stats         = FALSE,
  use_quantile_limits = TRUE,
  q                  = c(0.02, 0.98) # Recortar outliers extremos
)

