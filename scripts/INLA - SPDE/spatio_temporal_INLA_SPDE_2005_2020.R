# =============================================================================
# SCRIPT: spatio_temporal_INLA_SPDE_2005_2020.R
# =============================================================================
# PROPÓSITO:
#   Ajusta el modelo LGCP espacio-temporal completo para el catálogo sísmico
#   de Colombia 2005–2020, extendiendo el modelo espacial (año 2020) con una
#   estructura temporal AR(1) sobre el campo Gaussiano latente.
#
# DESCRIPCIÓN:
#   El catálogo de 16 años se agrupa en 8 períodos bienales
#   (2005-2006, 2007-2008, …, 2019-2020) definidos por la variable year_i:
#
#     year_i = ((year - 2005) %/% 2) + 1   →   valores 1 … 8
#
#   El campo Gaussiano latente u(s,t) sigue la estructura:
#
#     u(s, t) = ρ · u(s, t-1) + ε(s, t)
#     ε(s, t) ~ GF_Matérn(rango ρ_s, σ²)
#
#   con estructura de covarianza espaio-temporal dada por el producto de Kronecker:
#
#     Q = Q_T ⊗ Q_S
#
#   donde Q_S es la matriz de precisión SPDE y Q_T es la precisión AR(1).
#
#   El predictor lineal es:
#     log(λ(s,t)) = β₀ + β₁·volcanes + β₂·falla_inversa +
#                   β₃·falla_normal + β₄·isostasia + β₅·msnm + u(s,t)
#
# FUENTE (source):
#   - script_spatial/utils.R        (tema_tesis, book.mesh.dual, etc.)
#   - R/spde-book-functions.R       (book.mesh.dual)
#   - R/discrete_gradient.R         (helper visualización)
#
# INPUTS:
#   - covariables_rds/shapeZona_sp                                : polígono Colombia continental
#   - data_new/EventosColPointsPlanas31162005_2020_continental.gpkg : catálogo sísmico 2005-2020
#   - covariables_rds/*_im_scaled.rds                             : covariables escaldas (im de spatstat)
#
# OUTPUTS:
#   - pp.resCpv : objeto inla con el modelo espacio-temporal ajustado
#   - preds     : sf con el campo latente e intensidad predicha por período temporal
#   - Gráficos:
#       · p_spatial  : paneles del campo latente por período (plot_spatial_field)
#       · p_intensity : paneles de intensidad λ(s,t) por período (plot_intensity)
#       · p_temporal  : evolución temporal de la media del campo latente
#       · p_temporal_intensity: evolución temporal de la intensidad media
#       · panel_bootstrap : comparación bootstrap por período
#       · p_resumen  : diferencia en Log-Score por período
#
# PARÁMETROS A PRIORI (mismos que el modelo espacial):
#   mu_simulated    = -17.76713
#   scale_simulated = 2.400884
#   range_simulated = 110342.3 m
#
# NOTA COMPUTACIONAL:
#   El ajuste del modelo espacio-temporal (8 períodos × ~4309 vértices) tarda
#   aproximadamente 30-33 minutos en AMD Ryzen 9 7950X con 64 GB RAM.
#
# AUTOR: Jason Mauricio Romero Ríos
# UNIVERSIDAD: Universidad Distrital Francisco José de Caldas
# TESIS: Maestría en Ciencias de la Información y Comunicaciones – Geomática
# =============================================================================


# =============================================================================
# SECCIÓN 1: CARGA DE PAQUETES
# =============================================================================

library(INLA)        # Inferencia bayesiana INLA
library(sp)          # Clases espaciales legacy
library(ggplot2)     # Visualización
library(sf)          # Simple Features
library(spatial)     # Utilidades espaciales básicas
library(spData)      # Datos de ejemplo
library(spdep)       # Dependencia espacial
library(maps)        # Mapas base
library(gridExtra)   # Composición de gráficos
library(spatstat)    # Procesos puntuales
library(deldir)      # Voronoi (para book.mesh.dual)
library(raster)      # Rasters legacy
library(viridis)     # Paletas de color
library(terra)       # Rasters modernos
library(lubridate)   # Manejo de fechas
library(fmesher)     # Mallas FEM para INLA-SPDE

# Instalar fmesher si no está disponible
if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable")
}


# =============================================================================
# SECCIÓN 2: CONFIGURACIÓN DE RUTAS Y PARÁMETROS GLOBALES
# =============================================================================

setwd("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/")

# Directorio de covariables preprocesadas
files_rds <- "covariables_rds"

# Cargar funciones auxiliares del libro SPDE y utilitarias
source("R/spde-book-functions.R")   # book.mesh.dual()
source("R/discrete_gradient.R")
source("script_spatial/utils.R")   # create_ppp(), tema_tesis(), grafico_bootstrap(), etc.

# Ruta al catálogo sísmico completo (2005-2020)
path_file_seismic <- "data_new/EventosColPointsPlanas31162005_2020_continental.gpkg"


# =============================================================================
# SECCIÓN 3: PARÁMETROS A PRIORI (estimados mediante CNN)
# =============================================================================

# Intercepto: media a priori del log de la intensidad
mu_simulated <- -17.76713

# Desviación estándar a priori del campo Matérn
# PC-prior: P(sigma > scale_simulated) = 0.5
scale_simulated <- 2.400884

# Rango espacial a priori (~110 km)
# PC-prior: P(rango < range_simulated) = 0.5
range_simulated <- 110342.3

# Prior en efectos fijos: media e intercepto con baja precisión
control_fixed_list <- list(
  mean.intercept = mu_simulated,
  prec.intercept = 0.1           # varianza = 10 (prior débilmente informativo)
)


# =============================================================================
# SECCIÓN 4: CARGA DEL CATÁLOGO SÍSMICO Y CONSTRUCCIÓN DEL ÍNDICE TEMPORAL
# =============================================================================

# Cargar el polígono del área de estudio (Colombia continental, EPSG:3116)
shapeZona_sp <- readRDS(paste0(files_rds, "/shapeZona_sp"))

# Simplificar geometría para acelerar intersecciones (tolerancia 5 km)
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = TRUE)

# Leer el catálogo sísmico completo (2005-2020)
sismosSp <- st_read(path_file_seismic)
sismosSp$year <- sismosSp$YEAR

# Extraer coordenadas (EPSG:3116, metros)
sismosSp$X <- st_coordinates(sismosSp)[, 1]
sismosSp$Y <- st_coordinates(sismosSp)[, 2]

# Asegurar que year sea entero
sismosSp$year <- as.integer(sismosSp$year)

# --- Construcción del índice temporal bienal (year_i) ---
# Los 16 años (2005-2020) se agrupan en 8 períodos de 2 años cada uno:
#   year_i = 1 → 2005-2006
#   year_i = 2 → 2007-2008
#   ...
#   year_i = 8 → 2019-2020
# La fórmula ((year - 2005) %/% 2) + 1 produce enteros 1..8
sismosSp$year_i <- ((sismosSp$year - 2005) %/% 2) + 1

# Verificar distribución de eventos por período
table(sismosSp$year_i)

# Etiquetas descriptivas para gráficos (una por período)
period_labels <- c("2005-2008", "2009-2012", "2013-2016", "2017-2020")

# Convertir a entero explícito (por si quedó como factor)
sismosSp$year_i <- as.integer(as.character(sismosSp$year_i))


# =============================================================================
# SECCIÓN 5: CREACIÓN DEL PATRÓN PUNTUAL Y MALLA SPDE
# =============================================================================

# Crear el patrón puntual (ppp) dentro del polígono de estudio
creating_ppp <- create_ppp(shapeZona_sp, sismosSp)
p  <- creating_ppp$point_pattern       # ppp con todos los eventos (2005-2020)
xy <- creating_ppp$coordinates_eventos # Matriz de coordenadas

# --- Malla de Delaunay para SPDE ---
param_cutoff  <- 5000           # Distancia mínima entre vértices (m)
offset_param  <- c(100, 20000)  # Extensión interior/exterior (m)

shapeZona <- st_make_valid(shapeZona_sp)
shapeZona <- as_Spatial(shapeZona)
shapeZona <- st_as_sf(shapeZona)

# Tamaño máximo de arista: extensión total / 40 (más fino que el modelo espacial)
max.edge_params <- max(c(diff(range(xy[, 1])), diff(range(xy[, 2])))) / 40

# Malla 2D de Delaunay con todos los eventos como puntos de control
meshSismos <- inla.mesh.2d(
  boundary = shapeZona,
  loc      = cbind(xy),
  max.edge = c(1, 3) * max.edge_params,
  cutoff   = param_cutoff,
  offset   = offset_param,
  crs      = st_crs(shapeZona)
)

plot(meshSismos, main = "Malla SPDE – Colombia 2005-2020")


# =============================================================================
# SECCIÓN 6: PESOS DE INTEGRACIÓN (malla dual de Voronoi)
# =============================================================================

# Vértices del mesh y de todos los eventos
ptsSismos     <- as.matrix(xy)
mesh.ptsSismos <- as.matrix(meshSismos$loc[, 1:2])
allptsSismos  <- rbind(mesh.ptsSismos, ptsSismos)

nvSismos <- meshSismos$n    # Número de vértices del mesh
nSismos  <- nrow(ptsSismos) # Número total de eventos

# Construir la malla dual de Voronoi y asignar CRS
dmesh <- st_as_sf(book.mesh.dual(meshSismos))
st_crs(dmesh) <- st_crs(shapeZona_sp)

# Calcular el área de cada polígono de Voronoi dentro de la región de estudio
intersect_idx <- st_intersects(dmesh, shapeZona_sp, sparse = FALSE)[, 1]
dmesh_in      <- dmesh[intersect_idx, ]
intersections <- st_intersection(dmesh_in, shapeZona_sp)
areas         <- st_area(intersections)

# Vector de pesos: w[i] = área del polígono de Voronoi del vértice i en el estudio
# w[i] = 0 para vértices externos al dominio
w <- numeric(nrow(dmesh))
w[intersect_idx] <- as.numeric(areas)

# Verificación de áreas
total_area  <- as.numeric(st_area(shapeZona_sp))
sum_wSismos <- sum(w)
cat("Área total:", total_area, "\n")
cat("Suma de pesos Voronoi:", sum_wSismos, "\n")


# =============================================================================
# SECCIÓN 7: MALLA TEMPORAL 1D Y EXPOSICIÓN ESPACIO-TEMPORAL
# =============================================================================

# Malla 1D sobre los índices temporales (1 a k, con boundary = "free")
# Cada nodo de la malla temporal corresponde a un período bienal.
# boundary = "free" (sin condición de frontera periódica) es más rápido que "cyclic".
tmesh <- inla.mesh.1d(
  loc      = 1:length(unique(sismosSp$year_i)),
  boundary = "free"
)
k <- length(tmesh$loc)   # Número de períodos temporales (k = 8)

# Exposición espacio-temporal: los pesos de Voronoi se replican k veces
# st.vol[i + (t-1)*nv] = área del polígono de Voronoi del vértice i en el período t
# (la exposición es constante en el tiempo; solo varía espacialmente)
st.vol <- rep(w, k)


# =============================================================================
# SECCIÓN 8: CARGA Y EXTRACCIÓN DE COVARIABLES
# =============================================================================

# Cargar covariables preprocesadas como objetos im (imagen raster de spatstat)
topografia_im_scaled      <- readRDS(paste0(files_rds, "/topografia_im_scaled.rds"))
isostasia_im_scaled       <- readRDS(paste0(files_rds, "/isostasia_im_scaled.rds"))
volcanes_im_scaled        <- readRDS(paste0(files_rds, "/volcanes_im_scaled.rds"))
falla_sinestral_im_scaled <- readRDS(paste0(files_rds, "/sinestral_im_scaled.rds"))
falla_dextral_im_scaled   <- readRDS(paste0(files_rds, "/dextral_im_scaled.rds"))
falla_normal_im_scaled    <- readRDS(paste0(files_rds, "/normal_im_scaled.rds"))
falla_inversa_im_scaled   <- readRDS(paste0(files_rds, "/inversa_im_scaled.rds"))

# Lista de covariables del modelo espacio-temporal (5 covariables finales)
# Nota: falla_sinestral y falla_dextral se excluyen por multicolinealidad (VIF > 10)
covar <- list(
  isostasia    = isostasia_im_scaled,      # Anomalía isostática (mGal), escalada
  volcanes     = volcanes_im_scaled,        # Distancia a volcanes, escalada
  falla_inversa = falla_inversa_im_scaled, # Distancia a fallas inversas, escalada
  falla_normal  = falla_normal_im_scaled,  # Distancia a fallas normales, escalada
  msnm          = topografia_im_scaled     # Elevación (m s.n.m.), escalada
)

# --- Extracción de covariables en los vértices del mesh ---
# nearest.pixel() devuelve la celda de la imagen más cercana a cada punto.
# covar_mesh_once: valor de cada covariable en cada uno de los nvSismos vértices (una vez)
# covar_mesh_rep : los mismos valores replicados k veces para el stack espacio-temporal
covar_mesh_once <- purrr::map(covar, \(im) {
  px <- spatstat.geom::nearest.pixel(meshSismos$loc[, 1], meshSismos$loc[, 2], im)
  im[cbind(px$row, px$col)]
})
covar_mesh_rep <- purrr::map(covar_mesh_once, rep, times = k)  # k × nvSismos valores

# --- Extracción de covariables en las ubicaciones de todos los sismos ---
coords_all    <- cbind(sismosSp$X, sismosSp$Y)
covar_pts_all <- purrr::map(covar, \(im) {
  px <- spatstat.geom::nearest.pixel(coords_all[, 1], coords_all[, 2], im)
  im[cbind(px$row, px$col)]
})

# Añadir covariables como columnas al dataframe de sismos (para acceso por año)
kt <- length(unique(sismosSp$year_i))
for (nm in names(covar_pts_all)) sismosSp[[nm]] <- covar_pts_all[[nm]]


# =============================================================================
# SECCIÓN 9: ESPECIFICACIÓN DEL MODELO SPDE CON PC-PRIORS
# =============================================================================

# SPDE de Matérn con PC-priors informativos (estimados por CNN):
#   P(rango < range_simulated) = 0.5
#   P(sigma > scale_simulated) = 0.5
spdesismos <- inla.spde2.pcmatern(
  mesh        = meshSismos,
  prior.range = c(range_simulated, 0.5),
  prior.sigma = c(scale_simulated, 0.5)
)

m <- spdesismos$n.spde   # Número de vértices del mesh (= nvSismos)
k <- kt                   # Número de períodos temporales

# Índice SPDE espacio-temporal:
#   s        → índice del nodo espacial (1…m × k)
#   s.group  → índice del período temporal (1…k)
# El total de entradas es m × k (un campo espacial por período)
idx <- inla.spde.make.index(
  name    = "s",
  n.spde  = spdesismos$n.spde,
  n.group = length(unique(sismosSp$year_i))
)


# =============================================================================
# SECCIÓN 10: FUNCIÓN build_stack_year()
# =============================================================================

#' Construir el stack INLA para un período temporal específico
#'
#' Para cada período temporal (anio), genera el bloque del stack INLA que combina:
#'   - Los nv × k pesos de integración (vértices del mesh, todos los períodos)
#'   - Los eventos sísmicos del período anio (observaciones reales)
#'
#' La formulación LGCP requiere que en cada período se repliquen los k grupos
#' temporales para el mesh, y se añadan los eventos reales del período.
#'
#' @param anio integer. Índice del período temporal (1 a k).
#'
#' @return Un objeto inla.stack con:
#'   - data$y      : respuesta (0 en vértices de todos los períodos, 1 en eventos del año)
#'   - data$expect : exposición (st.vol = w replicado k veces, 0 en eventos)
#'   - A[[1]]      : Diagonal(k*m) para el mesh + Ast para los eventos del período
#'   - A[[2]]      : 1 para efectos fijos
#'   - effects[[1]]: idx SPDE espacio-temporal
#'   - effects[[2]]: covariables + intercepto a0
#'   - tag         : "year_X" donde X = anio
#'
#' @details
#'   La estructura del stack es:
#'     Filas 1 a k*m : vértices del mesh en todos los períodos (y=0, E=st.vol)
#'     Filas k*m+1 a k*m+n : eventos del período anio (y=1, E=0)
#'
build_stack_year <- function(anio) {

  filas <- which(sismosSp$year_i == anio)  # Índices de los sismos del período anio
  n     <- length(filas)                   # Número de eventos en el período
  nt    <- k * m + n                       # Tamaño total del stack para este año

  # Matriz de proyección A para los eventos del período anio
  # Tiene dimensión n × (k*m): proyecta cada evento al campo espacio-temporal
  Ast <- inla.spde.make.A(
    mesh       = meshSismos,
    loc        = cbind(sismosSp$X[filas], sismosSp$Y[filas]),
    group      = sismosSp$year_i[filas],   # Asigna cada evento a su grupo temporal
    n.group    = k,
    group.mesh = tmesh
  )

  # Respuesta y exposición
  y        <- rep(0:1, c(k * m, n))   # 0 en vértices, 1 en eventos del período
  expected <- c(st.vol, rep(0, n))    # Áreas de Voronoi (todos los períodos) + 0 para eventos

  # Covariables en el stack: mesh replicado k veces + eventos del período
  covariate_stack <- purrr::imap(covar_mesh_rep, \(mesh_vec, nm) {
    c(mesh_vec, covar_pts_all[[nm]][filas])
  })
  covariate_stack$a0 <- rep(1, nt)  # Intercepto (columna de unos)

  # Construir el stack INLA para este período
  inla.stack(
    data    = list(y = y, expect = expected),
    A       = list(rbind(Diagonal(k * m), Ast), 1),  # A para campo SPDE + A para efectos fijos
    effects = list(idx, covariate_stack),
    tag     = paste0("year_", anio)
  )
}


# =============================================================================
# SECCIÓN 11: CONSTRUCCIÓN DEL STACK GLOBAL
# =============================================================================

# Construir un stack por cada período temporal y combinarlos
years  <- sort(unique(sismosSp$year_i))      # Períodos: 1, 2, ..., k
stacks <- purrr::map(years, build_stack_year) # Lista de k stacks

# Combinar todos los stacks en uno solo (apilado verticalmente)
stk <- Reduce(inla.stack, stacks)
stk$tag <- "all_years"  # Etiqueta global del stack combinado

# Extraer el vector de exposición (áreas de Voronoi replicadas por período)
E <- inla.stack.data(stk)$expect


# =============================================================================
# SECCIÓN 12: FÓRMULA Y AJUSTE DEL MODELO ESPACIO-TEMPORAL
# =============================================================================

# --- Covariables incluidas en el modelo final ---
# (seleccionadas tras análisis de multicolinealidad VIF < 10)
use_covariables_model <- c("volcanes", "falla_inversa", "falla_normal", "isostasia")

# --- Fórmula del modelo espacio-temporal ---
# y ~ 1: intercepto global
# volcanes, falla_inversa, falla_normal, isostasia, msnm: efectos fijos geológicos
# f(s, model = spdesismos, group = s.group, control.group = list(model = 'ar1')):
#   campo Gaussiano latente SPDE con dependencia temporal AR(1) entre grupos
#   - s.group sigue el índice temporal (year_i = 1..k)
#   - model = 'ar1' especifica la estructura de dependencia temporal
formula_st <- y ~ 1 +
  volcanes + falla_inversa + falla_normal + isostasia + msnm +
  f(s,
    model         = spdesismos,
    group         = s.group,
    control.group = list(model = "ar1")   # Dependencia temporal AR(1) entre períodos
  )

# --- Ajuste del modelo con INLA ---
# int.strategy = "eb" (Empirical Bayes): integración más rápida para modelos grandes
time_M0 <- system.time({
  pp.resCpv <- inla(
    formula           = formula_st,
    family            = "poisson",
    data              = inla.stack.data(stk),
    control.predictor = list(
      A       = inla.stack.A(stk),
      compute = TRUE,
      link    = 1    # log-link
    ),
    control.fixed     = control_fixed_list,   # Prior informativo sobre el intercepto
    control.inla      = list(int.strategy = "eb"),
    E                 = E                     # Exposición (áreas de Voronoi)
  )
})

cat("Tiempo de ajuste espacio-temporal:", time_M0["elapsed"], "segundos\n")


# =============================================================================
# SECCIÓN 13: RESUMEN DE HIPERPARÁMETROS
# =============================================================================

# Nombres de los hiperparámetros del modelo espacio-temporal
names(pp.resCpv$marginals.hyperpar)

# Tabla completa de hiperparámetros
hyper_summary <- pp.resCpv$summary.hyperpar
print(hyper_summary)

# Extraer parámetros específicos del SPDE:
#   "Range for s"    → rango espacial ρ (metros)
#   "Stdev for s"    → desviación estándar del campo σ
#   "GroupRho for s" → coeficiente AR(1) temporal a ∈ (-1, 1)
range_summary <- hyper_summary["Range for s", ]
sigma_summary <- hyper_summary["Stdev for s", ]
rho_summary   <- hyper_summary["GroupRho for s", ]

# Tabla formateada de hiperparámetros principales
tabla_hiperparametros <- data.frame(
  Parametro = c("Range (ρ)", "Sigma (σ_ω)", "AR(1) coef (a)"),
  Mean  = c(range_summary$mean,          sigma_summary$mean,          rho_summary$mean),
  SD    = c(range_summary$sd,            sigma_summary$sd,            rho_summary$sd),
  Q0.025 = c(range_summary$`0.025quant`, sigma_summary$`0.025quant`,  rho_summary$`0.025quant`),
  Q0.5   = c(range_summary$`0.5quant`,   sigma_summary$`0.5quant`,    rho_summary$`0.5quant`),
  Q0.975 = c(range_summary$`0.975quant`, sigma_summary$`0.975quant`,  rho_summary$`0.975quant`)
)
print(tabla_hiperparametros)

# Efectos fijos (coeficientes β)
pp.resCpv$summary.fixed
pp.resCpv$summary.hyperpar


# =============================================================================
# SECCIÓN 14: EXTRACCIÓN DEL CAMPO ESPACIAL LATENTE POR PERÍODO
# =============================================================================

# El campo latente espacio-temporal u(s,t) tiene m × k entradas:
#   spatial_effect$mean[(t-1)*m + i] = media posterior de u(s_i, t)
spatial_effect <- pp.resCpv$summary.random$s

# Crear dataframe con el campo latente por nodo y período
df_spatial <- data.frame(
  node = rep(1:m, k),            # Índice del nodo espacial (1..m) repetido k veces
  time = rep(1:k, each = m),     # Índice del período temporal (1..k), m veces cada uno
  mean = spatial_effect$mean,    # Media posterior u(s, t)
  sd   = spatial_effect$sd,      # Desviación estándar posterior
  q025 = spatial_effect$`0.025quant`,
  q975 = spatial_effect$`0.975quant`
)

# --- Intensidad predicha λ̂(s,t) ---
# Los primeros k*m valores de summary.fitted.values corresponden a los vértices del mesh
mesh_idx        <- 1:(k * m)
fitted_intensity <- pp.resCpv$summary.fitted.values[mesh_idx, ]

df_intensity <- data.frame(
  node           = rep(1:m, k),
  time           = rep(1:k, each = m),
  intensity_mean = fitted_intensity$mean,
  intensity_sd   = fitted_intensity$sd,
  intensity_q025 = fitted_intensity$`0.025quant`,
  intensity_q975 = fitted_intensity$`0.975quant`
)

# --- Unir campo latente e intensidad con la geometría de la malla dual ---
# La malla dual (polígonos de Voronoi) se replica k veces para el mapa temporal
preds_base <- st_as_sf(dmesh)
preds <- preds_base[rep(1:nrow(preds_base), times = k), ]
rownames(preds) <- NULL

# Asignar período temporal, campo latente e intensidad
preds$time          <- rep(1:k, each = m)
preds$spatial_mean  <- df_spatial$mean
preds$spatial_sd    <- df_spatial$sd
preds$spatial_q025  <- df_spatial$q025
preds$spatial_q975  <- df_spatial$q975
preds$intensity_mean <- df_intensity$intensity_mean
preds$intensity_sd   <- df_intensity$intensity_sd

# Convertir a sf y asignar CRS
preds <- st_as_sf(preds)
st_crs(preds) <- st_crs(shapeZona_sp)


# =============================================================================
# SECCIÓN 15: FUNCIÓN plot_spatial_field()
# =============================================================================

# Calcular límites comunes para escala de color consistente entre paneles
spatial_limits <- quantile(preds$spatial_mean, c(0.02, 0.98), na.rm = TRUE)

# Etiquetas de períodos para los gráficos (extraídas del catálogo real)
period_labels <- paste0("Periodo ", 1:k)
if (exists("sismosSp") && "year" %in% names(sismosSp)) {
  years_by_period <- sismosSp %>%
    dplyr::group_by(year_i) %>%
    dplyr::summarise(years = paste(range(year), collapse = "-"), .groups = "drop") %>%
    dplyr::arrange(year_i)
  period_labels <- years_by_period$years
}


#' Graficar el campo latente espacial por período temporal
#'
#' Produce un panel de mapas (uno por período temporal) con el campo Gaussiano
#' latente u(s,t) coloreado con escala viridis común. La leyenda se extrae del
#' primer panel y se coloca debajo del arreglo.
#'
#' @param data sf. Objeto con columnas spatial_mean y time. Típicamente: preds.
#' @param shape sf. Polígono del área de estudio para el contorno del mapa.
#' @param limits numeric(2). Límites de la escala de color (mín, máx).
#' @param title_main character. Título principal del panel compuesto.
#'
#' @return cowplot::plot_grid. Objeto compuesto: título + paneles + leyenda.
#'
plot_spatial_field <- function(data, shape, limits, title_main) {

  data$time <- as.numeric(data$time)
  periodos  <- sort(unique(data$time))

  # Crear un mapa por cada período (sin leyenda individual)
  plots_list <- lapply(periodos, function(t) {
    df_t <- data[data$time == t, ]

    ggplot() +
      geom_sf(data = df_t, aes(fill = spatial_mean), color = NA) +
      geom_sf(data = shape, fill = NA, color = "black", linewidth = 0.3) +
      scale_fill_viridis_c(
        name   = "Campo\nlatente",
        limits = limits,
        oob    = scales::squish   # Recortar valores fuera de los límites (sin NA)
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text      = element_blank(),
        axis.ticks     = element_blank(),
        panel.grid     = element_blank(),
        plot.title     = element_text(hjust = 0.5, size = 11),
        plot.margin    = margin(2, 2, 2, 2),
        legend.position = "none"
      ) +
      ggtitle(period_labels[t])
  })

  # Extraer leyenda del primer panel (con la escala completa)
  p_legend <- ggplot() +
    geom_sf(data = data[data$time == periodos[1], ],
            aes(fill = spatial_mean), color = NA) +
    scale_fill_viridis_c(
      name   = "Campo\nlatente",
      limits = limits,
      oob    = scales::squish
    ) +
    theme_minimal() +
    theme(
      legend.position   = "bottom",
      legend.key.width  = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )
  legend_grob <- cowplot::get_legend(p_legend)

  # Combinar paneles en una grilla de 2 columnas
  panel_plots <- wrap_plots(plots_list, ncol = 2)

  # Título del panel compuesto
  titulo_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_main, fontface = "bold", size = 14, hjust = 0.5)

  # Composición final: título + paneles + leyenda
  cowplot::plot_grid(
    titulo_grob, panel_plots, legend_grob,
    ncol        = 1,
    rel_heights = c(0.05, 1, 0.08)
  )
}

# Ejecutar: campo latente u(s,t) para los k períodos
p_spatial <- plot_spatial_field(
  data       = preds,
  shape      = shapeZona_sp,
  limits     = spatial_limits,
  title_main = ""
)
print(p_spatial)


# =============================================================================
# SECCIÓN 16: FUNCIÓN plot_intensity()
# =============================================================================

#' Graficar la intensidad esperada λ̂(s,t) por período temporal
#'
#' Produce un panel de mapas con la intensidad predicha usando la paleta
#' plasma de viridis (adecuada para valores asimétricos como intensidades).
#'
#' @param data sf. Objeto con columnas intensity_mean y time. Típicamente: preds.
#' @param shape sf. Polígono del área de estudio.
#' @param limits numeric(2). Límites de la escala de color de intensidad.
#' @param title_main character. Título principal del panel.
#'
#' @return cowplot::plot_grid. Panel compuesto con todos los períodos.
#'
plot_intensity <- function(data, shape, limits, title_main) {

  data$time <- as.numeric(data$time)
  periodos  <- sort(unique(data$time))

  # Crear mapa por período (sin leyenda)
  plots_list <- lapply(periodos, function(t) {
    df_t <- data[data$time == t, ]

    ggplot() +
      geom_sf(data = df_t, aes(fill = intensity_mean), color = NA) +
      geom_sf(data = shape, fill = NA, color = "black", linewidth = 0.3) +
      scale_fill_viridis_c(
        name   = expression(lambda),   # Símbolo λ en la leyenda
        limits = limits,
        oob    = scales::squish,
        option = "plasma"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        panel.grid      = element_blank(),
        plot.title      = element_text(hjust = 0.5, size = 11),
        plot.margin     = margin(2, 2, 2, 2),
        legend.position = "none"
      ) +
      ggtitle(period_labels[t])
  })

  # Leyenda compartida (paleta plasma)
  p_legend <- ggplot() +
    geom_sf(data = data[data$time == periodos[1], ],
            aes(fill = intensity_mean), color = NA) +
    scale_fill_viridis_c(
      name   = expression(lambda),
      limits = limits,
      oob    = scales::squish,
      option = "plasma"
    ) +
    theme_minimal() +
    theme(
      legend.position   = "bottom",
      legend.key.width  = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )
  legend_grob <- cowplot::get_legend(p_legend)

  panel_plots <- wrap_plots(plots_list, ncol = 2)
  titulo_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_main, fontface = "bold", size = 14, hjust = 0.5)

  cowplot::plot_grid(
    titulo_grob, panel_plots, legend_grob,
    ncol        = 1,
    rel_heights = c(0.05, 1, 0.08)
  )
}

# Ejecutar: intensidad λ̂(s,t) para los k períodos
intensity_limits <- quantile(preds$intensity_mean, c(0.02, 0.98), na.rm = TRUE)

p_intensity <- plot_intensity(
  data       = preds,
  shape      = shapeZona_sp,
  limits     = intensity_limits,
  title_main = ""
)
print(p_intensity)


# =============================================================================
# SECCIÓN 17: EVOLUCIÓN TEMPORAL DE LA MEDIA DEL CAMPO LATENTE
# =============================================================================

library("dplyr")

# Resumen del campo latente promediado sobre todos los nodos, por período
temporal_summary <- preds %>%
  st_drop_geometry() %>%
  group_by(time) %>%
  summarise(
    mean_spatial    = mean(spatial_mean, na.rm = TRUE),   # Media espacial del campo latente
    sd_spatial      = sd(spatial_mean, na.rm = TRUE),     # Desviación estándar espacial
    mean_intensity  = mean(intensity_mean, na.rm = TRUE), # Intensidad media espacial
    .groups = "drop"
  )

color_principal <- "#4682B4"   # Azul acero para líneas principales
color_gris      <- "#7D8A96"   # Gris para barras de observados

# --- Gráfico 1: Evolución del campo latente medio a lo largo de los períodos ---
# La banda sombreada muestra ±1 desviación estándar espacial
p_temporal <- ggplot(temporal_summary, aes(x = time)) +
  geom_ribbon(
    aes(ymin = mean_spatial - sd_spatial, ymax = mean_spatial + sd_spatial),
    alpha = 0.2, fill = color_principal
  ) +
  geom_line(aes(y = mean_spatial), color = color_principal, linewidth = 1) +
  geom_point(aes(y = mean_spatial), color = color_principal, size = 3) +
  labs(x = "Periodo", y = "Media del campo latente") +
  scale_x_continuous(breaks = 1:k, labels = period_labels) +
  tema_tesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

print(p_temporal)

# --- Gráfico 2: Evolución de la intensidad media esperada ---
p_temporal_intensity <- ggplot(temporal_summary, aes(x = time)) +
  geom_line(aes(y = mean_intensity), color = color_principal, linewidth = 1) +
  geom_point(aes(y = mean_intensity), color = color_principal, size = 3) +
  labs(x = "Periodo", y = expression("Intensidad media (" * lambda * ")")) +
  scale_x_continuous(breaks = 1:k, labels = period_labels) +
  tema_tesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

print(p_temporal_intensity)

# Panel combinado: campo latente + intensidad, lado a lado
p_combined_temporal <- cowplot::plot_grid(
  p_temporal, p_temporal_intensity,
  ncol = 2, labels = c("", ""), label_size = 12
)
print(p_combined_temporal)


# =============================================================================
# SECCIÓN 18: ANÁLISIS DEL COEFICIENTE AR(1) TEMPORAL
# =============================================================================

rho_ar1 <- rho_summary$mean
rho_ci  <- c(rho_summary$`0.025quant`, rho_summary$`0.975quant`)

cat("\n========================================\n")
cat("ANÁLISIS DE PERSISTENCIA TEMPORAL (AR1)\n")
cat("========================================\n")
cat("Coeficiente AR(1):", round(rho_ar1, 4), "\n")
cat("IC 95%: [", round(rho_ci[1], 4), ",", round(rho_ci[2], 4), "]\n")
cat("Interpretación: Correlación entre períodos consecutivos =", round(rho_ar1, 2), "\n")
# La vida media de correlación indica cuántos períodos tarda en decaer la correlación a 1/e
cat("Vida media de la correlación:", round(-1 / log(rho_ar1), 2), "períodos\n")


# =============================================================================
# SECCIÓN 19: NÚMERO ESPERADO DE EVENTOS POR PERÍODO
# =============================================================================

# Integrar la intensidad sobre el espacio para obtener el número esperado de eventos
# N_esperado(t) = Σ_i λ̂(s_i, t) × w[i]
eventos_esperados <- preds %>%
  st_drop_geometry() %>%
  group_by(time) %>%
  summarise(
    N_esperado = sum(intensity_mean * st.vol[1:m], na.rm = TRUE),
    .groups    = "drop"
  )

# Eventos sísmicos observados por período
eventos_observados <- sismosSp %>%
  st_drop_geometry() %>%
  group_by(year_i) %>%
  summarise(N_observado = n(), .groups = "drop")

# Comparación: Observado vs Esperado por período
comparacion_eventos <- merge(eventos_esperados, eventos_observados,
                             by.x = "time", by.y = "year_i")
comparacion_eventos$Periodo <- period_labels
comparacion_eventos$Ratio  <- round(comparacion_eventos$N_observado /
                                      comparacion_eventos$N_esperado, 2)

print(comparacion_eventos[, c("Periodo", "N_esperado", "N_observado", "Ratio")])

# Gráfico comparativo: barras de observados + línea de esperados
p_eventos <- ggplot(comparacion_eventos, aes(x = time)) +
  geom_bar(aes(y = N_observado, fill = "Observado"), stat = "identity", alpha = 0.7) +
  geom_line(aes(y = N_esperado, color = "Esperado"), linewidth = 1.2) +
  geom_point(aes(y = N_esperado, color = "Esperado"), size = 3) +
  scale_fill_manual(values  = c("Observado" = color_gris)) +
  scale_color_manual(values = c("Esperado"  = color_principal)) +
  labs(x = "Periodo", y = "Número de eventos", fill = "", color = "") +
  scale_x_continuous(breaks = 1:k, labels = period_labels) +
  tema_tesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom")

print(p_eventos)


# =============================================================================
# SECCIÓN 20: FUNCIONES PARA VALIDACIÓN BOOTSTRAP POR PERÍODO
# =============================================================================

library(ggplot2)
library(dplyr)


#' Bootstrap de Log-Score para un período temporal
#'
#' Compara M2 (sin PC-prior) vs M3 (con PC-prior) usando bootstrap no paramétrico
#' del Log-Score de Poisson para el período indicado.
#'
#' Nota: En el contexto espacio-temporal, "M2" y "M3" se refieren a la comparación
#' del modelo sin prior informativo vs con prior informativo sobre el campo latente.
#' Los "datos" son los conteos por nodo del mesh para el período específico.
#'
#' @param data data.frame. Debe contener columnas: time, observed, expect_M2, expect_M3.
#' @param periodo integer. Índice del período temporal a analizar.
#' @param B integer. Número de remuestras bootstrap (default 5000).
#' @param seed integer. Semilla para reproducibilidad.
#'
#' @return numeric(B). Vector de diferencias δ = LS_M3 - LS_M2 en cada remuestra.
#'   δ > 0 indica que M3 es mejor que M2 en esa remuestra. NULL si no hay datos válidos.
#'
bootstrap_por_periodo <- function(data, periodo, B = 5000, seed = 2) {

  set.seed(seed)

  # Filtrar datos del período específico
  data_periodo <- data[data$time == periodo, ]

  # Seleccionar filas con datos válidos (observados y esperados finitos y positivos)
  ok <- is.finite(data_periodo$observed) &
    is.finite(data_periodo$expect_M2) &
    is.finite(data_periodo$expect_M3) &
    data_periodo$expect_M2 > 0 &
    data_periodo$expect_M3 > 0

  O  <- data_periodo$observed[ok]
  E1 <- data_periodo$expect_M2[ok]
  E2 <- data_periodo$expect_M3[ok]
  n  <- length(O)

  if (n == 0) {
    warning(paste("No hay datos válidos para el período", periodo))
    return(NULL)
  }

  # Bootstrap: remuestrear filas con reemplazo y calcular diferencia de Log-Scores
  delta <- numeric(B)
  for (b in 1:B) {
    idx     <- sample.int(n, n, replace = TRUE)
    ls1     <- sum(dpois(O[idx], lambda = E1[idx], log = TRUE))  # LS de M2
    ls2     <- sum(dpois(O[idx], lambda = E2[idx], log = TRUE))  # LS de M3
    delta[b] <- ls2 - ls1   # δ > 0 → M3 mejor que M2
  }

  return(delta)
}


#' Graficar la distribución bootstrap de δ = LS_M3 - LS_M2
#'
#' Histograma de densidad con la distribución de las B diferencias de Log-Score,
#' con líneas indicando la media y los cuantiles del IC 95%.
#'
#' @param delta_bootstrap numeric. Vector de B diferencias δ.
#' @param titulo character. Título del gráfico (default "").
#' @param color_fill character. Color de relleno del histograma (default "#4682B4").
#'
#' @return ggplot. Histograma de densidad de δ con anotaciones estadísticas.
#'
grafico_bootstrap <- function(delta_bootstrap, titulo = "", color_fill = "#4682B4") {

  df         <- data.frame(delta = delta_bootstrap)
  media_delta <- mean(delta_bootstrap)
  p_mejor    <- mean(delta_bootstrap > 0)       # P(M3 > M2)
  q025       <- quantile(delta_bootstrap, 0.025)
  q975       <- quantile(delta_bootstrap, 0.975)

  ggplot(df, aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50, fill = color_fill, alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red",     linewidth = 1) +
    geom_vline(xintercept = media_delta, linetype = "solid", color = "darkblue", linewidth = 0.8) +
    geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray40", linewidth = 0.7) +
    annotate("text", x = media_delta, y = Inf,
             label = paste0("Media = ", round(media_delta, 2)),
             vjust = 2, hjust = -0.1, size = 3, color = "darkblue") +
    annotate("text", x = max(delta_bootstrap) * 0.7, y = Inf,
             label = paste0("P(M3 > M2) = ", round(p_mejor, 3)),
             vjust = 4, size = 3.5) +
    labs(
      title = titulo,
      x     = expression(Delta ~ "Log-Score (M3 - M2)"),
      y     = "Densidad"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          panel.grid.minor = element_blank())
}


# =============================================================================
# SECCIÓN 21: FUNCIÓN eventos_por_nodo() – Conteo de sismos por nodo del mesh
# =============================================================================

#' Contar eventos sísmicos por nodo del mesh y período temporal
#'
#' Asigna cada evento al nodo del mesh más cercano usando la matriz de proyección
#' A (inla.spde.make.A), y devuelve una matriz de conteos [nodo × período].
#'
#' @param mesh inla.mesh.2d. La malla SPDE.
#' @param sismos sf o data.frame. Catálogo sísmico con columnas X, Y, year_i.
#' @param k integer. Número total de períodos temporales.
#'
#' @return matrix [mesh$n × k]. conteos[i, t] = número de eventos del período t
#'   asignados al nodo i del mesh.
#'
eventos_por_nodo <- function(mesh, sismos, k) {

  conteos <- matrix(0, nrow = mesh$n, ncol = k)

  for (t in 1:k) {
    idx_t <- which(sismos$year_i == t)
    if (length(idx_t) > 0) {
      A_t <- inla.spde.make.A(mesh = mesh, loc = cbind(sismos$X[idx_t], sismos$Y[idx_t]))
      for (i in seq_along(idx_t)) {
        nodo <- which.max(A_t[i, ])   # Nodo del mesh con mayor peso de proyección
        conteos[nodo, t] <- conteos[nodo, t] + 1
      }
    }
  }

  return(conteos)
}

# Calcular la matriz de conteos observados (m × k)
conteos_obs <- eventos_por_nodo(meshSismos, sismosSp, k)

# Crear dataframe para el análisis bootstrap
# observed: conteos reales por nodo y período
# expect  : intensidad predicha × peso de Voronoi = número esperado de eventos
spdf_sf <- data.frame(
  node     = rep(1:m, k),
  time     = rep(1:k, each = m),
  observed = as.vector(conteos_obs),
  expect   = preds$intensity_mean * w[1:m]   # Esperanza por nodo y período
)

# Nota: expect_M2 y expect_M3 necesitarían ajustar dos modelos separados.
# En este script se asume que la comparación se hace con respecto a un modelo base.
# El campo spdf_sf$expect se usa como "expect_M3" (modelo con prior informativo).


# =============================================================================
# SECCIÓN 22: BOOTSTRAP POR PERÍODO TEMPORAL
# =============================================================================

B         <- 5000
periodos  <- 1:k

resultados_bootstrap <- list()
graficos_bootstrap   <- list()

for (t in periodos) {
  cat("\nProcesando período", t, "...")

  delta_t <- bootstrap_por_periodo(spdf_sf, periodo = t, B = B)

  if (!is.null(delta_t)) {
    resultados_bootstrap[[t]] <- data.frame(
      periodo        = t,
      periodo_label  = period_labels[t],
      mean_delta     = mean(delta_t),
      p_M3_better    = mean(delta_t > 0),
      q025           = quantile(delta_t, 0.025),
      q975           = quantile(delta_t, 0.975)
    )

    graficos_bootstrap[[t]] <- grafico_bootstrap(delta_t, titulo = period_labels[t])
  }
}

# Tabla resumen del bootstrap por período
tabla_resultados <- do.call(rbind, resultados_bootstrap)
print(tabla_resultados)


# =============================================================================
# SECCIÓN 23: PANEL DE GRÁFICOS BOOTSTRAP
# =============================================================================

library(patchwork)
library(cowplot)

# Panel con todos los gráficos bootstrap (un histograma por período)
if (length(graficos_bootstrap) > 0) {
  panel_bootstrap <- wrap_plots(graficos_bootstrap, ncol = 2) +
    plot_annotation(
      title = "Comparación de Modelos por Bootstrap – Por Período",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  print(panel_bootstrap)
}

# Gráfico resumen: diferencia media por período con IC 95%
# Coloreado según P(M3 > M2): verde = M3 claramente mejor, rojo = M2 mejor
p_resumen <- ggplot(tabla_resultados, aes(x = periodo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.2, color = "gray50") +
  geom_point(aes(y = mean_delta, color = p_M3_better), size = 4) +
  geom_line(aes(y = mean_delta), color = "#4682B4", linewidth = 0.8) +
  scale_color_gradient2(
    low      = "#D73027", mid = "#FFFFBF", high = "#1A9850",
    midpoint = 0.5,
    name     = "P(M3 > M2)"
  ) +
  scale_x_continuous(breaks = periodos, labels = period_labels) +
  labs(
    title    = "Diferencia en Log-Score por Período",
    subtitle = "Valores positivos indican que M3 es mejor que M2",
    x        = "Período",
    y        = expression(Delta ~ "Log-Score (M3 - M2)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
print(p_resumen)


# =============================================================================
# SECCIÓN 24: BOOTSTRAP GLOBAL (TODOS LOS PERÍODOS JUNTOS)
# =============================================================================

set.seed(2)
B            <- 5000
delta_global <- numeric(B)

ok <- is.finite(spdf_sf$observed) &
  is.finite(spdf_sf$expect_M2) &
  is.finite(spdf_sf$expect_M3) &
  spdf_sf$expect_M2 > 0 &
  spdf_sf$expect_M3 > 0

O  <- spdf_sf$observed[ok]
E1 <- spdf_sf$expect_M2[ok]
E2 <- spdf_sf$expect_M3[ok]
n  <- length(O)

for (b in 1:B) {
  idx          <- sample.int(n, n, replace = TRUE)
  ls1          <- sum(dpois(O[idx], lambda = E1[idx], log = TRUE))
  ls2          <- sum(dpois(O[idx], lambda = E2[idx], log = TRUE))
  delta_global[b] <- ls2 - ls1
}

cat("\n========================================\n")
cat("BOOTSTRAP GLOBAL (TODOS LOS PERÍODOS)\n")
cat("========================================\n")
print(c(
  mean_delta  = mean(delta_global),
  p_M3_better = mean(delta_global > 0),
  q025        = quantile(delta_global, 0.025),
  q975        = quantile(delta_global, 0.975)
))

# Gráfico global de la distribución bootstrap
p_global <- grafico_bootstrap(delta_global, titulo = "Comparación Global de Modelos")
print(p_global)
