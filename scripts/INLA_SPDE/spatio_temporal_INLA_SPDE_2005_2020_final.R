# =============================================================================
# SCRIPT: spatio_temporal_INLA_SPDE_2005_2020_final.R
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
#   con estructura de covarianza espacio-temporal dada por el producto de Kronecker:
#
#     Q = Q_T ⊗ Q_S
#
#   donde Q_S es la matriz de precisión SPDE y Q_T es la precisión AR(1).
#
#     log(λ(s,t)) = β₀ + β₁·volcanes + β₂·falla_dextral + β₃·falla_normal + u(s,t)
#
# PARÁMETROS A PRIORI: estimados por la CNN (modelo "CNN + descriptores");
#   ver Sección 4.
#
# AUTOR: Jason Mauricio Romero Ríos
# UNIVERSIDAD: Universidad Distrital Francisco José de Caldas
# TESIS: Maestría en Ciencias de la Información y Comunicaciones – Geomática
# =============================================================================


# =============================================================================
# SECCIÓN 1: CARGA DE PAQUETES
# =============================================================================

# Instalar fmesher si no está disponible (repositorio de INLA)
if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable")
}

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
library(patchwork)   # Composición de paneles (wrap_plots, operador |)
library(cowplot)     # Composición de paneles (plot_grid, get_legend)
library(grid)        # Grobs de bajo nivel (rasterGrob)
library(dplyr)       # Manipulación de datos


# =============================================================================
# SECCIÓN 2: CONFIGURACIÓN DE RUTAS Y PARÁMETROS GLOBALES
# =============================================================================

setwd("/home/jasonromeroia/Documents/personal/Tesis_MCIC/lgcp_inla_spde_sismos/")

# Directorio de covariables preprocesadas
files_rds <- "covariables_rds"

# Cargar funciones auxiliares del libro SPDE y utilitarias
source("R/spde-book-functions.R")   # book.mesh.dual()
source("R/discrete_gradient.R")
source("R/utils.R")   # create_ppp(), tema_tesis(), grafico_bootstrap(), etc.

path_file_seismic <- "Data/gdf_espacial_2005_2020.gpkg"

# --- Directorio de salida de figuras ---
path_image_results <- "imagenes_doc"
dir_out_st <- file.path(path_image_results, "modelos_INLA_2005_2020")
dir.create(dir_out_st, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECCIÓN 3: FUNCIONES AUXILIARES DE FIGURAS
# =============================================================================

#' Guardar una figura en dir_out_st con dpi consistente (300)
save_fig <- function(p, name, w, h)
  ggsave(file.path(dir_out_st, name), plot = p, width = w, height = h, dpi = 300)

#' Leer un PNG y normalizarlo a RGB/RGBA (maneja escala de grises y gris+alpha)
read_png_rgba <- function(path) {
  x <- png::readPNG(path)
  d <- dim(x)
  if (length(d) == 2) {
    # gris -> replicar a RGB
    x <- array(rep(x, 3), dim = c(d[1], d[2], 3))
  } else if (d[3] == 2) {
    # gris + alpha -> expandir a RGBA
    rgb   <- x[, , 1]
    alpha <- x[, , 2]
    x2 <- array(0, dim = c(d[1], d[2], 4))
    x2[, , 1] <- rgb
    x2[, , 2] <- rgb
    x2[, , 3] <- rgb
    x2[, , 4] <- alpha
    x <- x2
  }
  x
}

# Logo Universidad Distrital
logo      <- read_png_rgba("imagenes_doc/logo_ud.png")
logo_grob <- rasterGrob(logo, interpolate = TRUE)

#' Extraer la leyenda de un ggplot 
extraer_leyenda <- function(p) {
  g <- tryCatch(cowplot::get_legend(p), error = function(e) NULL)
  if (is.null(g) || inherits(g, "zeroGrob")) {
    g <- cowplot::get_plot_component(p, "guide-box-right", return_all = TRUE)
    if (is.list(g) && !inherits(g, "grob")) g <- g[[1]]
  }
  g
}

# --- Leyenda salidas graficas
header_map_basic <- paste0(
  "Universidad Distrital\nFrancisco José de Caldas\n",
  "Maestría en Ciencias de la Información\n",
  "y Comunicaciones, Énfasis en Geomática\n",
  "Elaborado por: Jason Romero\n",
  "Sistema de Referencia: EPSG:3116"
)

# Eje vertical común: logo, texto y leyenda se centran sobre este valor.
CX <- 0.50

ANCHO_MAPA <- 0.84   # 

#' Componer un mapa facetado con panel lateral (logo + encabezado + leyenda)
#'
#'
#' @param p_mapa ggplot. Mapa
#' @return cowplot::plot_grid con el mapa (sin leyenda) y el panel lateral.
componer_mapa_panel_lateral <- function(p_mapa) {
  leyenda_grob <- extraer_leyenda(p_mapa + theme(legend.position = "right"))
  
  panel_lateral <- cowplot::ggdraw() +
    cowplot::draw_grob(
      logo_grob,
      x = CX, y = 0.99, width = 0.55, height = 0.22,
      hjust = 0.5, vjust = 1
    ) +
    cowplot::draw_text(
      header_map_basic,
      x = CX, y = 0.75,
      hjust = 0.5, vjust = 1,          # centra el bloque Y cada línea
      size = 6.8, lineheight = 1.15
    ) +
    cowplot::draw_grob(
      leyenda_grob,
      x = CX, y = 0.65, width = 0.92, height = 0.48,
      hjust = 0.5, vjust = 1
    ) +
    theme(plot.margin = margin(2, 2, 2, 0))
  
  cowplot::plot_grid(
    p_mapa + theme(legend.position = "none"), panel_lateral,
    ncol = 2, rel_widths = c(ANCHO_MAPA, 1 - ANCHO_MAPA)
  )
}

#' Alto (pulgadas) de una figura facetada según el aspecto del bbox
#'
#'
#' @param sf_obj sf. Objeto cuyo bbox define la relación de aspecto.
#' @param w_fig numeric. Ancho total de la figura (pulgadas).
#' @param ncols,nfilas integer. Disposición de las facetas.
#' @return numeric. Alto en pulgadas (+0.8 ~ franjas de facetas).
alto_figura_facetada <- function(sf_obj, w_fig, ncols, nfilas) {
  bb  <- sf::st_bbox(sf_obj)
  asp <- as.numeric((bb["ymax"] - bb["ymin"]) / (bb["xmax"] - bb["xmin"]))
  w_fig * ANCHO_MAPA * (asp * nfilas / ncols) + 0.8
}


# =============================================================================
# SECCIÓN 4: PARÁMETROS A PRIORI (estimados mediante CNN)
# =============================================================================

# --- Estimaciones puntuales de la CNN (modelo "CNN + descriptores", catálogo 2020) ---
# El campo espacial se ancla con la estimación de la CNN sobre 2020; el modelo
# espacio-temporal comparte esa estructura espacial.
scale_hat <- 72346     # scale Matern de spatstat (metros)
var_hat   <- 6.6284    # varianza del campo latente Y

# --- Conversión a la parametrización de INLA-SPDE (nu = 1) ---
# scale de spatstat != range de INLA: rGRFmatern usa z = (h/scale)*sqrt(2*nu),
# de donde kappa = sqrt(2*nu)/scale y range = sqrt(8*nu)/kappa = 2*scale.
# prior.sigma va sobre la desviación estándar marginal, no la varianza => sigma = sqrt(var).
range_cnn <- 3 * scale_hat     # 144692 m acá se multiplica por 3 para dar holgura
sigma_cnn <- sqrt(var_hat)     # 2.4679

range_simulated <- range_cnn
sigma_simulated <- sigma_cnn * 1.5 # acá se multiplica por 1.5 para dar holgura

# --- Hiperparámetros del PC-prior: ------------------
#   flexibilidad para el campo latente espacio-temporal.
#     P(range < range_prior_u) = alpha_pc
#   Sigma: anclaje en la cola con holgura sobre la estimación de la CNN.
#     P(sigma > sigma_prior_u) = alpha_pc_sigma

range_prior_u   <- range_simulated   # 217038 m
alpha_pc        <- 0.9               # P(range < range_prior_u) = 0.9
sigma_prior_u   <- sigma_simulated
alpha_pc_sigma  <- 0.01              # P(sigma > sigma_prior_u) = 0.01


# =============================================================================
# SECCIÓN 5: CARGA DEL CATÁLOGO SÍSMICO Y CONSTRUCCIÓN DEL ÍNDICE TEMPORAL
# =============================================================================

# Cargar el polígono del área de estudio (Colombia continental, EPSG:3116)
shapeZona_sp <- readRDS(paste0(files_rds, "/shapeZona_sp"))

# Simplificar geometría para acelerar intersecciones (tolerancia 5 km)
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = TRUE)

# Leer el catálogo sísmico completo (2005-2020)
sismosSp <- st_read(path_file_seismic)
sismosSp$year <- sismosSp$YEAR
year_min <- 2005
sismosSp <- subset(sismosSp, year >= year_min)
sismosSp <- subset(sismosSp, year <= 2020)

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
sismosSp$year_i <- ((sismosSp$year - year_min) %/% 2) + 1

# Verificar distribución de eventos por período
table(sismosSp$year_i)

sismosSp$year_i <- as.integer(as.character(sismosSp$year_i))


# =============================================================================
# SECCIÓN 6: CREACIÓN DEL PATRÓN PUNTUAL Y MALLA SPDE
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

max.edge_params <- max(c(diff(range(xy[, 1])), diff(range(xy[, 2])))) / 80

# Malla 2D de Delaunay SIN loc: el mesh se construye solo con la geometría de
# la zona (mesh independiente de los datos; los eventos se proyectan con A).
meshSismos <- inla.mesh.2d(
  boundary = shapeZona,
  max.edge = c(1, 3) * max.edge_params,
  cutoff   = param_cutoff,
  offset   = offset_param,
  crs      = st_crs(shapeZona)
)
cat("Número de vértices del mesh:", meshSismos$n, "\n")  # ~12900 con /80

plot(meshSismos, main = "Malla SPDE – Colombia 2005-2020")


# =============================================================================
# SECCIÓN 7: PESOS DE INTEGRACIÓN (malla dual de Voronoi)
# =============================================================================

nvSismos <- meshSismos$n   # Número de vértices del mesh

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
# SECCIÓN 8: MALLA TEMPORAL 1D Y EXPOSICIÓN ESPACIO-TEMPORAL
# =============================================================================

# Malla 1D sobre los índices temporales (1 a k, con boundary = "free")
# Cada nodo de la malla temporal corresponde a un período bienal.
# boundary = "free" (sin condición de frontera periódica) es más rápido que "cyclic".
tmesh <- inla.mesh.1d(
  loc      = 1:length(unique(sismosSp$year_i)),
  boundary = "free"
)
k <- length(tmesh$loc)   # Número de períodos temporales (k = 8)

# --- Exposición espacio-temporal = VOLUMEN espacio-temporal ------------------
# La verosimilitud LGCP aproxima ∫∫ λ(s,t) ds dt, así que el peso de cada
# pseudo-observación es (área de Voronoi del vértice i) × (ancho del nodo
# temporal t), NO solo el área.
#
# Orden: idx$s = rep(1:m, times = k) (nodo espacial varía primero) y
#        idx$s.group = rep(1:k, each = m), por eso rep(w, k) * rep(w.t, each = nv).
w.t    <- diag(inla.mesh.fem(tmesh)$c0)
st.vol <- rep(w, k) * rep(w.t, each = nvSismos)

cat("Anchos de nodo temporal (w.t):", paste(round(w.t, 3), collapse = " "), "\n")
cat("Volumen espacio-temporal total:", sum(st.vol),
    "| área x rango temporal:", sum(w) * diff(range(tmesh$loc)), "\n")


# =============================================================================
# SECCIÓN 9: CARGA Y EXTRACCIÓN DE COVARIABLES
# =============================================================================

# Cargar covariables preprocesadas como objetos im (imagen raster de spatstat)
topografia_im_scaled      <- readRDS(paste0(files_rds, "/topografia_im_scaled.rds"))
isostasia_im_scaled       <- readRDS(paste0(files_rds, "/isostasia_im_scaled.rds"))
volcanes_im_scaled        <- readRDS(paste0(files_rds, "/volcanes_im_scaled.rds"))
falla_sinestral_im_scaled <- readRDS(paste0(files_rds, "/sinestral_im_scaled.rds"))
falla_dextral_im_scaled   <- readRDS(paste0(files_rds, "/dextral_im_scaled.rds"))
falla_normal_im_scaled    <- readRDS(paste0(files_rds, "/normal_im_scaled.rds"))
falla_inversa_im_scaled   <- readRDS(paste0(files_rds, "/inversa_im_scaled.rds"))

# Catálogo de covariables DISPONIBLES (se extraen todas; cuáles entran al modelo
# Nota: falla_sinestral se excluye por multicolinealidad (VIF > 10).
covar <- list(
  isostasia    = isostasia_im_scaled,      # Anomalía isostática (mGal), escalada
  volcanes     = volcanes_im_scaled,        # Distancia a volcanes, escalada
  falla_inversa = falla_inversa_im_scaled, # Distancia a fallas inversas, escalada
  falla_normal  = falla_normal_im_scaled,  # Distancia a fallas normales, escalada
  msnm          = topografia_im_scaled,     # Elevación (m s.n.m.), escalada
  falla_dextral = falla_dextral_im_scaled
)

covar <- lapply(covar, function(X) blur(X, sigma = 5000, bleed = FALSE, normalise = TRUE))

# --- Extracción de covariables en los vértices del mesh ---
# nearest.pixel() devuelve la celda de la imagen más cercana a cada punto.
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
for (nm in names(covar_pts_all)) sismosSp[[nm]] <- covar_pts_all[[nm]]


# =============================================================================
# SECCIÓN 10: ESPECIFICACIÓN DEL MODELO SPDE CON PC-PRIORS
# =============================================================================

# SPDE de Matérn con PC-priors informativos derivados de la CNN.
spdesismos <- inla.spde2.pcmatern(
  mesh        = meshSismos,
  alpha       = 2,
  prior.range = c(range_prior_u, alpha_pc),
  prior.sigma = c(sigma_prior_u, alpha_pc_sigma)
)

m <- spdesismos$n.spde   # Número de vértices del mesh (= nvSismos)

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
# SECCIÓN 11: BLOQUE DE INTEGRACIÓN Y STACKS DE OBSERVACIÓN POR PERÍODO
# =============================================================================
# La verosimilitud LGCP aproximada sobre la malla dual tiene dos bloques:
#   (a) k*m pseudo-observaciones y = 0 en los nodos del mesh (todos los períodos),
#       con exposición E = área de Voronoi. Aproxima la integral ∫∫ λ(s,t) ds dt.
#   (b) n pseudo-observaciones y = 1 en los eventos reales, con E = 0.

covs_int    <- covar_mesh_rep    # covariables en los nodos, replicadas k veces
covs_int$a0 <- rep(1, k * m)     # intercepto

stk_int <- inla.stack(
  data    = list(y = rep(0, k * m), expect = st.vol),
  A       = list(Diagonal(k * m), 1),
  effects = list(idx, covs_int),
  tag     = "integration"
)

#' Construir el stack de observaciones de un período temporal
#'
#' Contiene únicamente los eventos reales del período (y = 1, E = 0). Los nodos
#' de integración viven en stk_int y no se repiten aquí.
#'
#' @param anio integer. Índice del período temporal (1 a k).
#'
#' @return Un objeto inla.stack con:
#'   - data$y      : 1 en cada evento del período
#'   - data$expect : 0 (los eventos no aportan al término integral)
#'   - A[[1]]      : Ast, proyección n × (k*m) de los eventos al campo
#'   - A[[2]]      : 1 para efectos fijos
#'   - tag         : "year_X" (lo usa el LCPO por período de la Sección 30)
#'
build_obs_year <- function(anio) {
  
  filas <- which(sismosSp$year_i == anio)  # Índices de los sismos del período anio
  n     <- length(filas)                   # Número de eventos en el período
  
  # Matriz de proyección A para los eventos del período anio
  Ast <- inla.spde.make.A(
    mesh       = meshSismos,
    loc        = cbind(sismosSp$X[filas], sismosSp$Y[filas]),
    group      = sismosSp$year_i[filas],   # Asigna cada evento a su grupo temporal
    n.group    = k,
    group.mesh = tmesh
  )
  
  covs    <- purrr::map(covar_pts_all, ~ .x[filas])
  covs$a0 <- rep(1, n)   # Intercepto
  
  inla.stack(
    data    = list(y = rep(1, n), expect = rep(0, n)),
    A       = list(Ast, 1),
    effects = list(idx, covs),
    tag     = paste0("year_", anio)
  )
}


# =============================================================================
# SECCIÓN 12: CONSTRUCCIÓN DEL STACK GLOBAL
# =============================================================================

years <- sort(unique(sismosSp$year_i))   # Períodos: 1, 2, ..., k

# Integración (una sola vez) + observaciones de cada período
stk <- Reduce(inla.stack, c(list(stk_int), purrr::map(years, build_obs_year)))

# Verificación: las filas de integración deben ser k*m exactamente, no k * (k*m)
n_int <- length(inla.stack.index(stk, "integration")$data)


stk$tag <- "all_years"  # Etiqueta global del stack combinado

# Extraer el vector de exposición (áreas de Voronoi en integración, 0 en eventos)
E <- inla.stack.data(stk)$expect


# =============================================================================
# SECCIÓN 13: GRID DE PREDICCIÓN (15 km) Y STACKS DE PREDICCIÓN POR AÑO
# =============================================================================

interval1 <- 5000

shapeZona_sp_interno <- st_buffer(shapeZona_sp, -1000)
grid_pred <- st_make_grid(shapeZona_sp_interno, cellsize = interval1,
                          what = "centers", square = TRUE)
grid_in   <- grid_pred[st_within(grid_pred, shapeZona_sp_interno, sparse = FALSE)]
pts.pred  <- as.data.frame(st_coordinates(grid_in))
colnames(pts.pred) <- c("x", "y")
n_pred    <- nrow(pts.pred)

# Ventana ampliada para evaluar covariables sin NA en el borde
expanded_window <- grow.rectangle(as.rectangle(p$window), 50000)
ppp.pred <- ppp(pts.pred$x, pts.pred$y, window = expanded_window)

# Un stack de predicción por año (mismo índice SPDE que las observaciones)
stacks_pred <- purrr::map(years, function(anio) {
  A.pred <- inla.spde.make.A(
    mesh    = meshSismos, loc = as.matrix(pts.pred),
    group   = rep(anio, n_pred), n.group = k, group.mesh = tmesh)
  # nearest.pixel: MISMO método que en el ajuste (Sección 9). im[ppp.pred]
  # devuelve NA fuera del soporte de la imagen y esas celdas se caerían de las
  # métricas con un patrón distinto por modelo.
  covs_pred <- purrr::map(covar, \(im) {
    px <- spatstat.geom::nearest.pixel(ppp.pred$x, ppp.pred$y, im)
    im[cbind(px$row, px$col)]
  })
  covs_pred$a0 <- rep(1, n_pred)
  # expect = 0 explícito: si se omite, inla.stack rellena con NA al combinar con
  # stk y E_join queda con NAs. Con y = NA la fila no aporta verosimilitud, pero
  # un E numérico evita depender de cómo INLA maneje NA en la exposición.
  inla.stack(data = list(y = rep(NA_real_, n_pred), expect = rep(0, n_pred)),
             A = list(A.pred, 1),
             effects = list(idx, covs_pred), tag = paste0("pred_year_", anio))
})
stk_pred <- Reduce(inla.stack, stacks_pred)

join.stack <- inla.stack(stk, stk_pred)
E_join     <- inla.stack.data(join.stack)$expect



# =============================================================================
# SECCIÓN 14: FÓRMULA Y AJUSTE DEL MODELO M1 (CON PRIOR INFORMATIVO CNN)
# =============================================================================

use_covariables_model <- c("volcanes", "falla_inversa", "falla_normal", "isostasia", "msnm")

cat("Covariables del modelo:", paste(use_covariables_model, collapse = " + "), "\n")

# --- Covariables NO lineales (OPCIONAL) --------
# Las covariables listadas aquí entran como f(inla.group(x), model = "rw2") en
# lugar de efecto lineal, en M0 y M1 por igual para que la comparación de
# priors siga teneindo sentido
use_nonlinear <- character(0)   # p.ej. c("msnm", "volcanes"); 

pc_prec_rw2 <- 'list(prec = list(prior = "pc.prec", param = c(1, 0.01)))'
terminos_nl <- sprintf(
  'f(inla.group(%s, n = 25, method = "quantile"), model = "rw2", scale.model = TRUE, constr = TRUE, hyper = %s)',
  use_nonlinear, pc_prec_rw2
)

covs_lineales <- setdiff(use_covariables_model, use_nonlinear)

# --- PC-prior sobre la correlación temporal AR(1) (pc.cor1, modelo base rho = 1) ---
# param = c(rho_u, rho_alpha) significa  P(rho > rho_u) = rho_alpha.
#   c(0.7, 0.7): 70% de probabilidad a priori de que rho > 0.7 -> favorece
#   persistencia temporal alta sin imponerla. Es parte del componente informativo
#   de M1; M0 mantiene el prior por defecto sobre rho.
rho_prior_params <- c(0.7, 0.7)
rho_u     <- rho_prior_params[1]
rho_alpha <- rho_prior_params[2]

# --- Fórmula del modelo espacio-temporal ---
# y ~ 1: intercepto global
# volcanes, falla_normal, isostasia, msnm: efectos fijos geológicos
# f(s, model = spdesismos, group = s.group, control.group = list(model = 'ar1')):
#   campo Gaussiano latente SPDE con dependencia temporal AR(1) entre grupos
#   - s.group sigue el índice temporal (year_i = 1..k)
#   - hyper$rho: PC-prior informativo sobre la correlación temporal
# El término f(s, ...) se inyecta como texto para que reformulate() lo respete.
h.spec <- list(rho = list(prior = "pc.cor1", param = c(rho_u, rho_alpha)))
termino_spde_m1 <- paste0(
  'f(s, model = spdesismos, ',
  'group = s.group, ',
  'control.group = list(model = "ar1", hyper = ',
  paste(deparse(h.spec, width.cutoff = 500), collapse = ""),
  '))'
)

formula_st <- reformulate(c(covs_lineales, terminos_nl, termino_spde_m1),
                          response = "y", intercept = TRUE)
print(formula_st)

# (para el caché RDS de los modelos) ---------------------
firma_nucleo <- paste(c(
  paste(sort(use_covariables_model), collapse = "-"),
  if (length(use_nonlinear) > 0) paste0("nl", paste(sort(use_nonlinear), collapse = "-")),
  sprintf("k%d_m%d_n%d", k, m, nrow(sismosSp))
), collapse = "_")
firma_stack <- sprintf("%s_i%d", firma_nucleo, interval1)
ruta_M1_rds <- file.path(files_rds, sprintf("pp.resM1_2005_2020_%s.rds", firma_stack))

control_fixed_M0 <- list(mean = 0, prec = 1)

# --- Ajuste del modelo M1: CON prior informativo (valores CNN)
# int.strategy = "eb" y cmin = 0: mismos controles que el modelo espacial 2020.
start.time <- Sys.time()
mu_simulated <- -21.6429
control_fixed_list <- list(
  mean.intercept = mu_simulated,
  prec.intercept = 0.01           # varianza = 10 (prior débilmente informativo)
)
control_fixed_M1 <- c(control_fixed_M0, control_fixed_list)
print(start.time)
if (file.exists(ruta_M1_rds)) {
  obj_M1 <- readRDS(ruta_M1_rds)
  # saveRDS guarda el objeto inla directo; $pp.resM1 sobre un objeto inla
  # devuelve NULL y el modelo quedaba NULL al releer el caché. Se acepta
  # también el formato antiguo list(pp.resM1 = ...).
  t_M1 <- attr(obj_M1, "tiempo_ajuste_s")
  pp.resM1 <- if (!is.null(obj_M1$pp.resM1)) obj_M1$pp.resM1 else obj_M1
  segundos_M1 <- if (is.null(t_M1)) NA_real_ else round(t_M1, 2)
  cat("Modelo M1 leido de", ruta_M1_rds, "\n")
} else {
  time_M1 <- system.time({
    pp.resM1 <- inla(
      formula           = formula_st,
      family            = "poisson",
      data              = inla.stack.data(join.stack),
      control.predictor = list(
        A       = inla.stack.A(join.stack),
        compute = TRUE,
        link    = 1    # log-link
      ),
      control.fixed     = control_fixed_M1,   # Prior beta N(0,1): acota la divergencia EB
      control.inla      = list(int.strategy = "eb", cmin = 0),
      control.compute   = list(config = TRUE, cpo = TRUE),
      E                 = E_join                 # Exposición (obs) + NA en predicción
    )
  })
  segundos_M1 <- time_M1[["elapsed"]]
  attr(pp.resM1, "tiempo_ajuste_s") <- segundos_M1
  saveRDS(pp.resM1, ruta_M1_rds)
  cat("Modelo M1 guardado en", ruta_M1_rds, "\n")
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

cat("Tiempo de ajuste M1 (con prior):", segundos_M1, "segundos\n")


# =============================================================================
# SECCIÓN 15: AJUSTE DEL MODELO M0 (SIN PRIOR INFORMATIVO)
# =============================================================================
# Modelo de referencia para cuantificar el aporte del prior de la CNN: misma
# estructura (covariables + AR(1)) pero con un SPDE de Matérn estándar (prior
# vago). Comparte el mismo stack.

# Mismas covariables que M1 (use_covariables_model), mismo AR(1).
spde_noprior <- inla.spde2.matern(
  mesh        = meshSismos,
  alpha       = 2
)

termino_spde_m0 <- paste0(
  'f(s, model = spde_noprior, ',
  'group = s.group, control.group = list(model = "ar1"))'
)

formula_st_m0 <- reformulate(c(covs_lineales, terminos_nl, termino_spde_m0),
                             response = "y", intercept = TRUE)
print(formula_st_m0)

ruta_M0_rds <- file.path(files_rds, sprintf("pp.resM0_2005_2020_%s.rds", firma_stack))
start.time <- Sys.time()
print(start.time)
if (file.exists(ruta_M0_rds)) {
  obj_M0 <- readRDS(ruta_M0_rds)
  t_M0 <- attr(obj_M0, "tiempo_ajuste_s")
  pp.resM0 <- if (!is.null(obj_M0$pp.resM0)) obj_M0$pp.resM0 else obj_M0
  segundos_M0 <- if (is.null(t_M0)) NA_real_ else round(t_M0, 2)
  cat("Modelo M0 leido de", ruta_M0_rds, "\n")
} else {
  time_M0 <- system.time({
    pp.resM0 <- inla(
      formula           = formula_st_m0,
      family            = "poisson",
      data              = inla.stack.data(join.stack),
      control.predictor = list(
        A       = inla.stack.A(join.stack),
        compute = TRUE,
        link    = 1
      ),
      control.fixed     = control_fixed_M0,   # Prior beta N(0,1), sin intercepto CNN
      control.inla      = list(int.strategy = "eb", cmin = 0),
      control.compute   = list(config = TRUE, cpo = TRUE),
      E                 = E_join
    )
  })
  segundos_M0 <- time_M0[["elapsed"]]
  cat("Tiempo de ajuste M0 (sin prior):", segundos_M0, "segundos\n")
  
  attr(pp.resM0, "tiempo_ajuste_s") <- segundos_M0
  saveRDS(pp.resM0, ruta_M0_rds)
  cat("Modelo M0 (sin prior) guardado en", ruta_M0_rds, "\n")
}

pp.resM1$summary.fixed
pp.resM0$summary.fixed


# =============================================================================
# SECCIÓN 16: TIEMPO DE EJECUCIÓN Y RESÚMENES DE EFECTOS FIJOS
# =============================================================================

df_tiempos <- data.frame(
  Modelo   = factor(c("M0", "M1"), levels = c("M0", "M1")),
  segundos = c(segundos_M0, segundos_M1)
)
df_tiempos$minutos <- df_tiempos$segundos / 60

p_tiempos <- ggplot(df_tiempos, aes(x = Modelo, y = minutos)) +
  geom_col(fill = "steelblue", width = 0.85) +
  labs(x = "Modelo", y = "Tiempo de ejecución (minutos)") +
  theme_bw(base_size = 12)

save_fig(p_tiempos, "time_execution_temporal.png", 4.5, 4.5)

# --- Exportación de resúmenes de efectos fijos a .txt ------------------------
# Ambos modelos incluyen covariables: se exporta la tabla completa
# (Variable, Media, IC 95 %) para M0 y M1.

#' Exportar el resumen de efectos fijos de un modelo INLA a un .txt
#'
#' @param res objeto inla.
#' @param modelo character. Etiqueta del modelo ("M0", "M1").
#' @param archivo character. Ruta del .txt de salida.
#' @param solo_intercepto logical. Si TRUE exporta únicamente el intercepto.
#' @param digitos integer. Decimales de redondeo.
exportar_resumen_fijos <- function(res, modelo, archivo,
                                   solo_intercepto = FALSE, digitos = 3) {
  sf_ <- as.data.frame(res$summary.fixed)
  sf_$Variable <- rownames(sf_)
  if (solo_intercepto) {
    sf_ <- sf_[sf_$Variable %in% c("(Intercept)", "b0", "intercept"), , drop = FALSE]
  }
  fmt <- paste0("%.", digitos, "f")
  tabla <- data.frame(
    Variable  = sf_$Variable,
    Media     = sprintf(fmt, sf_$mean),
    `IC 95 %` = sprintf(paste0("[", fmt, ", ", fmt, "]"),
                        sf_[["0.025quant"]], sf_[["0.975quant"]]),
    check.names = FALSE
  )
  con <- file(archivo, open = "w", encoding = "UTF-8")
  writeLines(c(paste("Modelo", modelo),
               capture.output(print(tabla, row.names = FALSE))), con)
  close(con)
  cat("Resumen de efectos fijos exportado:", archivo, "\n")
}

exportar_resumen_fijos(pp.resM0, "M0",
                       file.path(dir_out_st, "summary_fijos_M0_temporal.txt"))
exportar_resumen_fijos(pp.resM1, "M1",
                       file.path(dir_out_st, "summary_fijos_M1_temporal.txt"))


# =============================================================================
# SECCIÓN 17: RESUMEN POSTERIOR DE LOS HIPERPARÁMETROS
# =============================================================================
# Los hiperparámetros del campo (rango, sigma, rho AR1).

get_hyper_row <- function(res, patron) {
  nm <- grep(patron, rownames(res$summary.hyperpar), value = TRUE)
  stopifnot(length(nm) >= 1)
  res$summary.hyperpar[nm[1], , drop = FALSE]
}
range_summary <- get_hyper_row(pp.resM1, "^Range")
sigma_summary <- get_hyper_row(pp.resM1, "Stdev|Stddev|Sigma")
rho_summary   <- get_hyper_row(pp.resM1, "GroupRho")

tabla_hiperparametros <- data.frame(
  Parametro      = c("Range (m)", "Sigma (σ_ω)", "AR(1) coef (a)"),
  Ancla_prior_M1 = c(range_prior_u, sigma_prior_u, rho_u),
  Post_media_M1  = c(range_summary$mean, sigma_summary$mean, rho_summary$mean),
  Post_q025      = c(range_summary$`0.025quant`, sigma_summary$`0.025quant`,
                     rho_summary$`0.025quant`),
  Post_q975      = c(range_summary$`0.975quant`, sigma_summary$`0.975quant`,
                     rho_summary$`0.975quant`)
)
print(tabla_hiperparametros)

tabla_final_hiperparametros <- tabla_hiperparametros %>%
  mutate(IC_95 = paste0("[", format(round(Post_q025, 4), nsmall = 4), ", ",
                        format(round(Post_q975, 4), nsmall = 4), "]")) %>%
  select(Parametro, Media = Post_media_M1, IC_95)

tabla_final_hiperparametros$Media <- round(tabla_final_hiperparametros$Media, 3)
write.table(tabla_final_hiperparametros,
            file.path(dir_out_st, "tabla_hiperparametros_temporal.txt"),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)


# =============================================================================
# SECCIÓN 18: EXTRACCIÓN DEL CAMPO ESPACIAL LATENTE POR PERÍODO
# =============================================================================

# El campo latente espacio-temporal u(s,t) tiene m × k entradas:
#   spatial_effect$mean[(t-1)*m + i] = media posterior de u(s_i, t)
spatial_effect <- pp.resM1$summary.random$s

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
# Filas del bloque de integración (los k*m vértices del mesh). Se piden por tag
# en vez de asumir 1:(k*m): así sigue siendo correcto aunque cambie el orden con
# que se apilan los stacks.
mesh_idx <- inla.stack.index(join.stack, "integration")$data
stopifnot(length(mesh_idx) == k * m)
fitted_intensity <- pp.resM1$summary.fitted.values[mesh_idx, ]

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

# --- Etiquetas de períodos para los gráficos (extraídas del catálogo real) ---
period_labels <- paste0("Periodo ", 1:k)
if (exists("sismosSp") && "year" %in% names(sismosSp)) {
  # st_drop_geometry() ANTES del group_by: sin él, dplyr::summarise sobre un sf
  # ejecuta un st_union por grupo (decenas de miles de puntos x k grupos) solo
  # para construir etiquetas de texto.
  years_by_period <- sismosSp %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(year_i) %>%
    dplyr::summarise(years = paste(range(year), collapse = "-"), .groups = "drop") %>%
    dplyr::arrange(year_i)
  period_labels <- years_by_period$years
}


# =============================================================================
# SECCIÓN 19: MAPA DEL CAMPO LATENTE POR PERÍODO (plot_spatial_field)
# =============================================================================

#' Graficar el campo latente espacial por período temporal
#'
#' Produce un panel de mapas (uno por período temporal) con el campo Gaussiano
#' latente u(s,t)
#'
plot_spatial_field <- function(data, shape, limits, title_main,
                               breaks = seq(-7, 7, 2),
                               ncol = 2,
                               period_labels = NULL) {
  
  data$time <- as.numeric(data$time)
  periodos  <- sort(unique(data$time))
  
  # Si no se proporcionan period_labels, usar los valores por defecto
  if (is.null(period_labels)) {
    period_labels <- paste("Período", periodos)
  }
  
  # Función para crear un panel individual
  make_panel <- function(t, label) {
    df_t <- data[data$time == t, ]
    
    p <- ggplot() +
      geom_sf(data = df_t, aes(fill = spatial_mean), color = NA) +
      geom_sf(data = shape, fill = NA, color = "black", linewidth = 0.3) +
      scale_fill_viridis_c(
        limits = limits,
        breaks = breaks,
        oob = scales::squish,
        guide = "none"  # Sin leyenda en cada panel
      ) +
      coord_sf(expand = FALSE) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "none"
      ) +
      ggtitle(label)
    
    return(p)
  }
  
  # Crear todos los paneles con sus labels
  plots_list <- mapply(make_panel, periodos, period_labels, SIMPLIFY = FALSE)
  
  # Combinar paneles y agregar leyenda colectada
  spatial_plot <- wrap_plots(plots_list, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") &
    scale_fill_viridis_c(
      name = "Campo\nlatente",
      limits = limits,
      breaks = breaks,
      oob = scales::squish,
      guide = guide_colourbar(
        nbin = 500,
        raster = TRUE,
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 0.5,
        barwidth = 1.2,
        barheight = 15,
        direction = "vertical",
        title.position = "top",
        title.theme = element_text(hjust = 0.5, size = 9)
      )
    ) &
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  # Agregar título general
  spatial_plot <- spatial_plot +
    plot_annotation(
      title = title_main,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    )
  
  return(spatial_plot)
}

# Ejecutar: campo latente u(s,t) para los k períodos
p_spatial <- plot_spatial_field(
  data       = preds,
  shape      = shapeZona_sp,
  limits = c(-7, 7),
  title_main = "",
  ncol = 3,
  breaks = seq(-7, 7, 2),
  period_labels = period_labels
)
print(p_spatial)
# Figura usada en la tesis (fig:sptial_effectM3_temporal)
save_fig(p_spatial, "sptial_effectM3_temporal.png", w = 10, h = 12)


# =============================================================================
# SECCIÓN 20: MAPAS DE INTENSIDAD ESPERADA λ̂(s,t) POR PERÍODO
# =============================================================================
# Panel de mapas facetados con la intensidad predicha, recortada a los
# percentiles 2-98 para el color, y compuesto con el panel lateral
# institucional (logo + encabezado + leyenda).

# --- 1. Recorte de la intensidad para la escala de color ---------------------
CLIP_GLOBAL <- TRUE

if (CLIP_GLOBAL) {
  q <- quantile(preds$intensity_mean, c(0.02, 0.98), na.rm = TRUE)
  preds_clip <- preds %>%
    mutate(intensity_plot = pmin(pmax(intensity_mean, q[1]), q[2]))
} else {
  preds_clip <- preds %>%
    group_by(time) %>%
    mutate(
      intensity_plot = pmin(
        pmax(intensity_mean, quantile(intensity_mean, 0.02, na.rm = TRUE)),
        quantile(intensity_mean, 0.98, na.rm = TRUE)
      )
    ) %>%
    ungroup()
}

lims_pred <- range(preds_clip$intensity_plot, na.rm = TRUE)

# Suponiendo que preds_clip$time toma valores 1..8
period_labels <- c("2005-2006","2007-2008","2009-2010","2011-2012",
                   "2013-2014","2015-2016","2017-2018","2019-2020")

names(period_labels) <- sort(unique(as.character(preds_clip$time)))
# --- 2. Mapa facetado --------------------------------------------------------
p_preds_mapa <- ggplot(preds_clip) +
  geom_sf(aes(fill = intensity_plot), color = NA) +
  geom_sf(data = shapeZona_sp, fill = NA, color = "grey25", linewidth = 0.3) +
  facet_wrap(
    ~ time, ncol = 4,
    labeller = as_labeller(period_labels)
  ) +
  scale_fill_viridis_c(
    name   = "Intensidad\nesperada (λ)",
    limits = lims_pred,
    oob    = scales::squish,
    labels = scales::label_scientific(digits = 2),
    guide  = guide_colourbar(
      barwidth     = 0.8,
      barheight    = 9,
      frame.colour = "black",
      ticks.colour = "black",
      title.position = "top",
      title.theme  = element_text(size = 8.5, face = "bold", hjust = 0),
      label.theme  = element_text(size = 7)
    )
  ) +
  coord_sf(expand = FALSE)  +
  tema_tesis() +
  theme(
    strip.text    = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.15, "lines"),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.margin   = margin(2, 0, 2, 0)
  )

# --- 3. Composición con panel lateral y exportación --------------------------
p_final <- componer_mapa_panel_lateral(p_preds_mapa)

w_fig <- 15   # pulgadas
h_fig <- alto_figura_facetada(preds_clip, w_fig, ncols = 4, nfilas = 2)

# Figura usada en la tesis (fig:intensidad_temporal)
save_fig(p_final, "result_intensidades_tempora.png", w_fig, h_fig)


# =============================================================================
# SECCIÓN 21: EVOLUCIÓN TEMPORAL DE LA MEDIA DEL CAMPO LATENTE
# =============================================================================

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
# Figura usada en la tesis (fig:intensidad_campolatente_temporal)
save_fig(p_combined_temporal, "campo_latente_temporal.png", w = 11, h = 5)


# =============================================================================
# SECCIÓN 22: ANÁLISIS DEL COEFICIENTE AR(1) TEMPORAL
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
# SECCIÓN 23: NÚMERO ESPERADO DE EVENTOS POR PERÍODO
# =============================================================================

# Integrar la intensidad sobre el ESPACIO para obtener el nº esperado de eventos
# del período: N_esperado(t) = Σ_i λ̂(s_i, t) × w[i]
#
# Se usa `w` (áreas de Voronoi), NO st.vol: st.vol ya incluye el ancho temporal
# w.t y aquí solo se integra en el espacio dentro de cada período. Antes se
# escribía st.vol[1:m], que coincidía con w solo porque st.vol no llevaba w.t.
eventos_esperados <- preds %>%
  st_drop_geometry() %>%
  group_by(time) %>%
  summarise(
    N_esperado = sum(intensity_mean * w, na.rm = TRUE),
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
# Validación observado vs esperado por período
save_fig(p_eventos, "eventos_obs_vs_esperado.png", w = 9, h = 5)


# =============================================================================
# SECCIÓN 24: PRIOR (CNN) vs POSTERIOR DE LOS HIPERPARÁMETROS
# =============================================================================
# Muestra si los datos actualizaron el PC-prior anclado con la CNN o si el prior
# domina. Tres paneles: rango espacial, sigma (ambos con la densidad analítica del
# PC-prior superpuesta) y el coeficiente AR(1) temporal (solo posterior: pcmatern
# no fija su prior). Modelo único: pp.resM1.
# Con hiperparámetros FIJOS (Sección 14) no existe posterior de hiperparámetros:
# la sección completa se omite.
if (length(pp.resM1$marginals.hyperpar) > 0) {
  
  # --- Densidades analíticas del PC-prior (las mismas que usa inla.spde2.pcmatern) ---
  #   Rango (d = 2):  pi(r)     = lambda1 * r^-2 * exp(-lambda1 / r),  lambda1 = -log(alpha)*r0
  #   Sigma:          pi(sigma) = lambda2 * exp(-lambda2 * sigma),     lambda2 = -log(alpha)/sigma0
  # IMPORTANTE: se usan range_prior_u / sigma_prior_u (Sección 4), que son EXACTAMENTE
  # los valores pasados a inla.spde2.pcmatern(). Antes esta curva se dibujaba con
  # sigma_cnn en vez de sigma_prior_u = 2*sigma_cnn, así que la figura mostraba un
  # prior el doble de concentrado que el realmente usado en el ajuste.
  lambda1_range <- -log(alpha_pc)       * range_prior_u
  lambda2_sigma <- -log(alpha_pc_sigma) / sigma_prior_u
  pc_prior_range <- function(r) lambda1_range * r^(-2) * exp(-lambda1_range / r)
  pc_prior_sigma <- function(s) lambda2_sigma * exp(-lambda2_sigma * s)
  
  # --- Helper: marginal posterior de un hiperparámetro por patrón de nombre ---
  get_hyper_marg <- function(res, patron) {
    nm <- grep(patron, names(res$marginals.hyperpar), value = TRUE)
    if (length(nm) == 0) return(NULL)
    as.data.frame(res$marginals.hyperpar[[nm[1]]])   # columnas x, y
  }
  
  # Paleta: prior gris discontinuo vs posterior azul de marca; se refuerza con linetype.
  col_hp <- c("PC-prior (CNN)" = "#7D8A96", "Posterior" = "#4682B4")
  lty_hp <- c("PC-prior (CNN)" = "dashed",  "Posterior" = "solid")
  
  # --- RANGE (se grafica en km; densidad reescalada por el cambio de unidad) ---
  mr <- get_hyper_marg(pp.resM1, "^Range")
  r_max  <- max(mr$x, 3 * range_prior_u, na.rm = TRUE)
  r_grid <- seq(1, r_max, length.out = 400)
  df_range <- rbind(
    data.frame(x = r_grid / 1000, y = pc_prior_range(r_grid) * 1000, curva = "PC-prior (CNN)"),
    data.frame(x = mr$x / 1000,   y = mr$y * 1000,                   curva = "Posterior")
  )
  
  p_range <- ggplot(df_range, aes(x, y, color = curva, linetype = curva)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = range_cnn / 1000, color = "#2d2d2d", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x = range_cnn / 1000, y = Inf, label = "",
             vjust = 1.5, hjust = -0.05, size = 3, color = "#2d2d2d") +
    scale_color_manual(values = col_hp) +
    scale_linetype_manual(values = lty_hp) +
    coord_cartesian(xlim = c(0, min(r_max / 1000, 600))) +
    labs(title = "", x = "Rango (km)", y = "Densidad", color = NULL, linetype = NULL) +
    tema_tesis()
  
  # --- SIGMA (desviación estándar marginal del campo) ---
  ms <- get_hyper_marg(pp.resM1, "Stdev|Stddev|Sigma")
  s_max  <- max(ms$x, 1.5 * sigma_prior_u, na.rm = TRUE)
  s_grid <- seq(1e-3, s_max, length.out = 400)
  df_sigma <- rbind(
    data.frame(x = s_grid, y = pc_prior_sigma(s_grid), curva = "PC-prior (CNN)"),
    data.frame(x = ms$x,   y = ms$y,                   curva = "Posterior")
  )
  
  p_sigma <- ggplot(df_sigma, aes(x, y, color = curva, linetype = curva)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = sigma_cnn, color = "#2d2d2d", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x = sigma_cnn, y = Inf, label = "",
             vjust = 1.5, hjust = -0.05, size = 3, color = "#2d2d2d") +
    scale_color_manual(values = col_hp) +
    scale_linetype_manual(values = lty_hp) +
    labs(title = "",
         x = expression(sigma), y = "Densidad", color = NULL, linetype = NULL) +
    tema_tesis()
  
  # --- AR(1) temporal (GroupRho): solo posterior, sin prior analítico ---
  mrho <- get_hyper_marg(pp.resM1, "GroupRho")
  rho_med <- rho_summary$`0.5quant`
  p_rho <- ggplot(mrho, aes(x, y)) +
    geom_line(color = "#4682B4", linewidth = 0.8) +
    geom_vline(xintercept = rho_med, color = "#2d2d2d", linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = rho_med, y = Inf, label = paste0("Mediana ≈ ", round(rho_med, 3)),
             vjust = 1.5, hjust = 1.05, size = 3, color = "#2d2d2d") +
    labs(title = "", x = expression(a), y = "Densidad") +
    tema_tesis()
  
  p_hyper <- (p_range | p_sigma | p_rho) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "",
                    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom")
  
  print(p_hyper)
  save_fig(p_hyper, "prior_vs_posterior_hyper_temporal.png", w =12, h = 5)
  
} else {
  cat("[aviso] Hiperparámetros fijos: se omite la figura prior vs posterior (Sección 24).\n")
}


# =============================================================================
# SECCIÓN 25: FOREST PLOT DE EFECTOS FIJOS (covariables geológicas)
# =============================================================================
# Media posterior e IC 95% de cada covariable. Se excluye el intercepto (vive en
# otra escala y aplastaría los IC de las covariables).

sf_ <- as.data.frame(pp.resM1$summary.fixed)
sf_$efecto <- rownames(sf_)
sf_ <- sf_[!sf_$efecto %in% c("(Intercept)", "b0", "intercept"), , drop = FALSE]
df_fixed <- data.frame(
  efecto = sf_$efecto,
  media  = sf_$mean,
  lo     = sf_[["0.025quant"]],
  hi     = sf_[["0.975quant"]]
)

p_forest <- ggplot(df_fixed, aes(x = media, y = efecto)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#2d2d2d", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  color = "#4682B4", linewidth = 0.7, fatten = 3) +
  labs(title = "Efectos fijos de las covariables geológicas",
       subtitle = "Media posterior e IC 95% (log-escala); intercepto omitido",
       x = "Coeficiente", y = NULL) +
  tema_tesis()

print(p_forest)
save_fig(p_forest, "forest_efectos_fijos.png", w = 9, h = 5)


# =============================================================================
# SECCIÓN 26: EFECTOS NO LINEALES ESTIMADOS (solo si use_nonlinear no es vacío)
# =============================================================================
# Curva f(x) posterior (media + IC 95%) de cada covariable rw2. El ID del
# efecto es el punto medio del bin de inla.group (covariable escalada).
if (length(use_nonlinear) > 0) {
  plots_nl <- lapply(use_nonlinear, function(nm) {
    id <- grep(nm, names(pp.resM1$summary.random), value = TRUE, fixed = TRUE)
    id <- setdiff(id, "s")[1]
    df_nl <- pp.resM1$summary.random[[id]]
    ggplot(df_nl, aes(x = ID, y = mean)) +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  fill = "#4682B4", alpha = 0.2) +
      geom_line(color = "#4682B4", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#2d2d2d") +
      labs(title = nm, x = paste(nm, "(escalada)"),
           y = "Efecto en log-intensidad") +
      tema_tesis()
  })
  p_nl <- wrap_plots(plots_nl, ncol = min(2, length(plots_nl)))
  print(p_nl)
  save_fig(p_nl, "efectos_no_lineales.png",
           w = 9, h = 4.5 * ceiling(length(plots_nl) / 2))
}


# =============================================================================
# SECCIÓN 27: FUNCIÓN DE CORRELACIÓN MATÉRN ESPACIAL POSTERIOR
# =============================================================================
# Convierte el "rango" (abstracto) en una curva de correlación vs distancia.
# nu = 1:  rho(h) = (kappa*h) * K_1(kappa*h),  con kappa = sqrt(8)/range.

matern_corr_nu1 <- function(h, rng) {
  kappa <- sqrt(8) / rng
  z <- kappa * h
  ifelse(h <= 0, 1, z * besselK(z, 1))
}
rng_med <- range_summary$`0.5quant`
rng_lo  <- range_summary$`0.025quant`
rng_hi  <- range_summary$`0.975quant`

h_grid <- seq(0, 1.5 * rng_hi, length.out = 400)
c_lo   <- matern_corr_nu1(h_grid, rng_lo)
c_hi   <- matern_corr_nu1(h_grid, rng_hi)
df_corr <- data.frame(
  h   = h_grid / 1000,
  med = matern_corr_nu1(h_grid, rng_med),
  lo  = pmin(c_lo, c_hi),
  hi  = pmax(c_lo, c_hi)
)

p_corr <- ggplot(df_corr, aes(h)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4682B4", alpha = 0.20) +
  geom_line(aes(y = med), color = "#4682B4", linewidth = 0.9) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "#2d2d2d") +
  geom_vline(xintercept = rng_med / 1000, linetype = "dashed", color = "#2d2d2d") +
  annotate("text", x = rng_med / 1000, y = 0.92,
           label = paste0("Rango ≈ ", round(rng_med / 1000), " km"),
           hjust = -0.05, size = 3.2, color = "#2d2d2d") +
  labs(title = "Función de correlación Matérn espacial posterior",
       subtitle = "ν = 1; banda = IC 95% del rango; línea punteada: ρ = 0.1",
       x = "Distancia (km)", y = expression(rho(h))) +
  tema_tesis()

print(p_corr)
save_fig(p_corr, "correlacion_matern.png", w = 8, h = 5)


# =============================================================================
# SECCIÓN 28: DECAIMIENTO DE LA CORRELACIÓN TEMPORAL AR(1)
# =============================================================================
# El AR(1) del campo latente implica corr(t, t+l) = a^l entre períodos bianuales.
# Traduce el coeficiente 'a' (Sección 22) en una curva de persistencia con IC 95%
# y la vida media (l tal que a^l = 1/e => l = -1/log(a)).

rho_lo <- rho_summary$`0.025quant`
rho_hi <- rho_summary$`0.975quant`
lags   <- 0:(k - 1)
df_ar1 <- data.frame(
  lag = lags,
  med = rho_ar1^lags,
  lo  = pmin(rho_lo^lags, rho_hi^lags),
  hi  = pmax(rho_lo^lags, rho_hi^lags)
)
vida_media <- -1 / log(rho_ar1)   # en períodos (cada período = 2 años)

p_ar1 <- ggplot(df_ar1, aes(x = lag)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4682B4", alpha = 0.20) +
  geom_line(aes(y = med), color = "#4682B4", linewidth = 0.9) +
  geom_point(aes(y = med), color = "#4682B4", size = 2.5) +
  geom_hline(yintercept = exp(-1), linetype = "dotted", color = "#2d2d2d") +
  geom_vline(xintercept = vida_media, linetype = "dashed", color = "#E69F00", linewidth = 0.7) +
  annotate("text", x = vida_media, y = 0.95,
           label = paste0("Vida media ≈ ", round(vida_media, 1), " períodos (",
                          round(2 * vida_media), " años)"),
           hjust = -0.03, size = 3.2, color = "#2d2d2d") +
  scale_x_continuous(breaks = lags) +
  labs(title = "Persistencia temporal del campo latente (AR(1))",
       subtitle = expression(paste("corr(t, t+l) = ", a^l, ";  banda = IC 95%;  línea punteada: ", rho == 1/e)),
       x = "Rezago l (períodos bianuales)", y = "Correlación") +
  tema_tesis()

print(p_ar1)
save_fig(p_ar1, "correlacion_temporal_ar1.png", w = 8, h = 5)


# =============================================================================
# SECCIÓN 29: EXCEDENCIA DE lambda(s,t) POR PERÍODO
# =============================================================================
# Mapa de riesgo P(lambda(s,t) > umbral) por período.
# Reutiliza `preds` (sf, m x k) de la Sección 18 (media y sd de la intensidad).

# Etiquetas de período para las facetas
period_lab_fun <- function(v) period_labels[as.integer(v)]

borde <- st_transform(st_geometry(shapeZona_sp), st_crs(preds))

# --- Probabilidad de excedencia P(lambda > lambda0) por período ---

m_int   <- preds$intensity_mean
sd_int  <- preds$intensity_sd
s2_ln   <- log1p((sd_int / m_int)^2)
mu_ln   <- log(m_int) - s2_ln / 2
lambda0 <- as.numeric(quantile(m_int, 0.90, na.rm = TRUE))   # decil superior global
preds$exceed <- 1 - pnorm((log(lambda0) - mu_ln) / sqrt(s2_ln))

p_exceed <- ggplot(preds) +
  geom_sf(aes(fill = exceed), color = NA) +
  geom_sf(data = borde, fill = NA, color = "grey30", linewidth = 0.2) +
  facet_wrap(~ time, ncol = 4, labeller = as_labeller(period_lab_fun)) +
  scale_fill_viridis_c(option = "rocket", name = "P", direction = -1, limits = c(0, 1)) +
  labs(title = "") +
  tema_tesis() +
  theme(axis.text = element_blank(), panel.grid = element_blank())


p_final_exceed <- componer_mapa_panel_lateral(p_exceed)

w_fig <- 15   # pulgadas
h_fig <- alto_figura_facetada(preds_clip, w_fig, ncols = 4, nfilas = 2)

save_fig(p_final_exceed, "excedencia_periodo.png", w_fig, h_fig)


# =============================================================================
# SECCIÓN 30: LOG-SCORE Y LCPO POR AÑO (M0 sin prior vs M1 con prior)
# =============================================================================
# Métricas calculadas sobre el GRID DE PREDICCIÓN de 15 km (Sección 13), igual que
# el script original. Para cada año:
#   - intensidad λ̂ = exp(predictor lineal) en el tag "pred_year_X",
#   - observado = nº de eventos del año por celda (st_intersects),
#   - esperado  = λ̂ × área de la celda,
#   - Log-Score = sum(dpois(obs, esp, log))  (mayor es mejor),
#   - LCPO      = -sum(log(cpo)) de las filas de observación (tag "year_X"), menor mejor.

poisson_log_score <- function(observed, expected, mask = NULL, eps = 1e-12) {
  ok <- is.finite(observed) & is.finite(expected) & expected > 0
  if (!is.null(mask)) ok <- ok & mask   # máscara común entre modelos (comparabilidad)
  sum(dpois(observed[ok], lambda = pmax(expected[ok], eps), log = TRUE))
}

# Estadísticos por año sobre el grid de predicción, para un modelo ya ajustado
calcular_estadisticos_por_anio <- function(result_inla, grid_in, join_stack, years) {
  pts.pred_ <- as.data.frame(st_coordinates(grid_in))
  names(pts.pred_) <- c("x", "y")
  por_anio <- list()
  for (anio in years) {
    idx_year  <- inla.stack.index(join_stack, paste0("pred_year_", anio))$data
    # exp(E[eta]) es la MEDIANA posterior de lambda, no la media: falta el factor
    # exp(sigma_eta^2/2) de Jensen. Por eso n_expected queda sistemáticamente por
    # debajo de n_observed, y ese sesgo NO debe leerse como falta de ajuste.
    # Se conserva porque el sesgo es el mismo en M0 y M1, así que el Log-Score
    # comparado y el bootstrap de delta siguen siendo válidos.
    pred_mean <- exp(result_inla$summary.linear.predictor[idx_year, "mean"])
    
    pts_year <- pts.pred_
    pts_year$pred_intensity <- pred_mean
    coordinates(pts_year) <- ~x + y
    proj4string(pts_year) <- st_crs(shapeZona_sp)$proj4string
    gridded(pts_year) <- TRUE
    spdf_sf_year <- st_as_sf(as(as(pts_year, "SpatialPixelsDataFrame"),
                                "SpatialPolygonsDataFrame"))
    spdf_sf_year <- st_transform(spdf_sf_year, crs = st_crs(sismosSp))
    
    sismos_year <- sismosSp[sismosSp$year_i == anio, ]
    spdf_sf_year$observed <- lengths(st_intersects(spdf_sf_year, sismos_year))
    spdf_sf_year$area     <- as.numeric(st_area(spdf_sf_year))
    spdf_sf_year$expected <- spdf_sf_year$pred_intensity * spdf_sf_year$area
    
    por_anio[[as.character(anio)]] <- list(
      year       = anio,
      spdf_sf    = spdf_sf_year,
      n_observed = sum(spdf_sf_year$observed),
      n_expected = sum(spdf_sf_year$expected, na.rm = TRUE),
      log_score  = poisson_log_score(spdf_sf_year$observed, spdf_sf_year$expected)
    )
  }
  # LCPO por año: -sum(log(cpo)) sobre las filas de observación (tag "year_X").
  # Se excluyen los CPO poco fiables marcados por INLA (cpo$failure > 0). La máscara
  # COMÚN entre modelos se impone después (ver bucle de máscara común), aquí basta el
  # filtro de fiabilidad por modelo.
  # Con cpo desactivado (corrida de prueba) el modelo no trae $cpo: el LCPO queda
  # NA por año y las secciones aguas abajo lo omiten.
  if (is.null(result_inla$cpo$cpo)) {
    lcpo_por_anio <- rep(NA_real_, length(years))
  } else {
    lcpo_por_anio <- sapply(years, function(anio) {
      idx_obs <- inla.stack.index(join_stack, paste0("year_", anio))$data
      cpo  <- result_inla$cpo$cpo[idx_obs]
      fail <- result_inla$cpo$failure[idx_obs]
      ok   <- is.finite(cpo) & cpo > 0
      if (!is.null(fail)) ok <- ok & is.finite(fail) & fail == 0
      -sum(log(cpo[ok]))
    })
  }
  list(por_anio = por_anio, lcpo_por_anio = lcpo_por_anio)
}

stats_M0 <- calcular_estadisticos_por_anio(pp.resM0, grid_in, join.stack, years)
stats_M1 <- calcular_estadisticos_por_anio(pp.resM1, grid_in, join.stack, years)

# --- Máscara COMÚN M0/M1 por año para Log-Score y LCPO ---------------------------
# Tanto -Σ log(CPO) como Σ log dpois() escalan con el nº de puntos, así que para que
# la comparación M0 vs M1 sea justa se evalúan ambos modelos sobre los MISMOS índices
# válidos por año. Para el LCPO se excluyen además los CPO poco fiables (failure > 0).
for (i in seq_along(years)) {
  a <- as.character(years[i])
  
  # (a) Log-Score: celdas del grid con esperanza > 0 en AMBOS modelos
  sf0 <- stats_M0$por_anio[[a]]$spdf_sf
  sf1 <- stats_M1$por_anio[[a]]$spdf_sf
  ok_ls <- is.finite(sf0$observed) &
    is.finite(sf0$expected) & sf0$expected > 0 &
    is.finite(sf1$expected) & sf1$expected > 0
  stats_M0$por_anio[[a]]$log_score <- poisson_log_score(sf0$observed, sf0$expected, mask = ok_ls)
  stats_M1$por_anio[[a]]$log_score <- poisson_log_score(sf1$observed, sf1$expected, mask = ok_ls)
  
  # (b) LCPO: pseudo-observaciones (tag "year_X") válidas y fiables en AMBOS
  # modelos. Solo aplica si los modelos se ajustaron con cpo (corrida final).
  if (!is.null(pp.resM0$cpo$cpo) && !is.null(pp.resM1$cpo$cpo)) {
    idx  <- inla.stack.index(join.stack, paste0("year_", years[i]))$data
    cpo0 <- pp.resM0$cpo$cpo[idx];  fail0 <- pp.resM0$cpo$failure[idx]
    cpo1 <- pp.resM1$cpo$cpo[idx]; fail1 <- pp.resM1$cpo$failure[idx]
    ok0  <- is.finite(cpo0) & cpo0 > 0
    ok1  <- is.finite(cpo1) & cpo1 > 0
    if (!is.null(fail0)) ok0 <- ok0 & is.finite(fail0) & fail0 == 0
    if (!is.null(fail1)) ok1 <- ok1 & is.finite(fail1) & fail1 == 0
    ok_cpo <- ok0 & ok1
    stats_M0$lcpo_por_anio[i] <- -sum(log(cpo0[ok_cpo]))
    stats_M1$lcpo_por_anio[i] <- -sum(log(cpo1[ok_cpo]))
    
    cat(sprintf("Año %s | LS celdas comunes: %d/%d | CPO comunes: %d/%d (failures M0=%d, M1=%d)\n",
                a, sum(ok_ls), length(ok_ls), sum(ok_cpo), length(ok_cpo),
                if (is.null(fail0)) 0L else sum(fail0 > 0, na.rm = TRUE),
                if (is.null(fail1)) 0L else sum(fail1 > 0, na.rm = TRUE)))
  } else {
    cat(sprintf("Año %s | LS celdas comunes: %d/%d | LCPO omitido (cpo desactivado, corrida de prueba)\n",
                a, sum(ok_ls), length(ok_ls)))
  }
}

# --- Tabla resumen por año y modelo ---
tabla_resumen <- function(stats, model_name) {
  data.frame(
    year       = as.integer(names(stats$por_anio)),
    n_observed = sapply(stats$por_anio, `[[`, "n_observed"),
    n_expected = round(sapply(stats$por_anio, `[[`, "n_expected"), 1),
    log_score  = round(sapply(stats$por_anio, `[[`, "log_score"), 2),
    lcpo       = round(stats$lcpo_por_anio, 2),
    model      = model_name, row.names = NULL)
}
tabla_completa <- rbind(tabla_resumen(stats_M0, "M0"), tabla_resumen(stats_M1, "M1"))
print(tabla_completa)

# --- Barras agrupadas por año: Log-Score y LCPO, M0 gris vs M1 azul ---
col_mod <- c("M0" = "#7D8A96", "M1" = "#4682B4")
tabla_completa$model <- factor(tabla_completa$model, levels = c("M0", "M1"))
lab_mod <- c("M0 (Sin Prior)", "M1 (Con Prior)")

# Log-Score: valores NEGATIVOS -> la etiqueta cuelga bajo la punta de la barra
p_ls <- ggplot(tabla_completa, aes(x = factor(year), y = log_score, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "white", linewidth = 0.4, alpha = 0.9) +
  geom_text(aes(label = round(log_score)), position = position_dodge(0.8),
            angle = 90, hjust = 1.15, size = 2.6, fontface = "bold", color = "#2d2d2d") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.02))) +  # aire abajo para el texto
  scale_fill_manual(values = col_mod, labels = lab_mod, name = NULL) +
  labs(title = "Log Score", subtitle = "", x = "Año", y = "Log Score") +
  tema_tesis()

# LCPO: valores POSITIVOS -> la etiqueta sale hacia arriba desde la punta
p_lcpo <- ggplot(tabla_completa, aes(x = factor(year), y = lcpo, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "white", linewidth = 0.4, alpha = 0.9) +
  geom_text(aes(label = round(lcpo, 1)), position = position_dodge(0.8),
            angle = 90, hjust = -0.1, size = 2.6, fontface = "bold", color = "#2d2d2d") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +  # aire arriba para el texto
  scale_fill_manual(values = col_mod, labels = lab_mod, name = NULL) +
  labs(title = "LCPO", subtitle = "", x = "Año", y = "LCPO") +
  tema_tesis()

paneles_val <- p_ls | p_lcpo

p_val <- paneles_val +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))) &
  theme(legend.position = "bottom")

print(p_val)
# Figura usada en la tesis (fig:log_score_lcpo_temporal)
save_fig(p_val, "log_score_lcpo_temporal.png", w = 12, h = 6)


# =============================================================================
# SECCIÓN 31: BOOTSTRAP DEL LOG-SCORE POR AÑO (M1 vs M0)
# =============================================================================
# δ = LogScore_M1 - LogScore_M0 por año, remuestreando CELDAS del grid de predicción
# (con reemplazo), igual que el script original. δ > 0 favorece al modelo con prior.
# Estilo tesis: grafico_bootstrap() de R/utils.R (histograma + densidad + banda IC +
# caja de estadísticos con μ, IC y P(M1 > M0)).

set.seed(1)
B <- 5000   # Número de réplicas bootstrap

boot_list  <- vector("list", length(years))
tabla_boot <- vector("list", length(years))
for (i in seq_along(years)) {
  anio <- years[i]
  sf0 <- stats_M0$por_anio[[as.character(anio)]]$spdf_sf
  sf1 <- stats_M1$por_anio[[as.character(anio)]]$spdf_sf
  ok  <- is.finite(sf0$observed) & is.finite(sf0$expected) & is.finite(sf1$expected) &
    sf0$expected > 0 & sf1$expected > 0
  O <- sf0$observed[ok]; E0 <- sf0$expected[ok]; E1 <- sf1$expected[ok]; n <- length(O)
  delta <- numeric(B)
  for (b in 1:B) {
    j <- sample.int(n, n, replace = TRUE)
    delta[b] <- sum(dpois(O[j], E1[j], log = TRUE)) -   # Log-Score M1 (con prior)
      sum(dpois(O[j], E0[j], log = TRUE))     # Log-Score M0 (sin prior)
  }
  gb <- grafico_bootstrap(delta,
                          modelo_base      = "M0",
                          modelo_comparado = "M1",
                          titulo           = period_labels[i],
                          stats_size       = 2.6)
  boot_list[[i]]  <- gb$plot
  tabla_boot[[i]] <- cbind(periodo = period_labels[i], gb$stats)
}


panel_boot <- wrap_plots(boot_list, ncol = 4) +
  plot_annotation(
    title = "",
    theme = theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5)))

print(panel_boot)

save_fig(panel_boot, "bootstrap_temporal.png", w = 16, h = 8)
