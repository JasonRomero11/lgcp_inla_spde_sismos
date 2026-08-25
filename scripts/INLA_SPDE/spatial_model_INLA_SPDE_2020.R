# Modelo LGCP espacial con INLA-SPDE sobre el catálogo sísmico de Colombia, 2020.
# Se ajustan cuatro variantes para separar el efecto del prior informativo
# (estimado con la CNN) del efecto de las covariables geológicas:
#
#   M0  prior vago,  sin covariables      M1  prior CNN,  sin covariables
#   M2  prior vago,  con covariables      M3  prior CNN,  con covariables
#
# Jason Mauricio Romero Ríos
# Universidad Distrital Francisco José de Caldas
# Maestría en Ciencias de la Información y las Comunicaciones – Geomática


library(INLA)
library(sp)
library(ggplot2)
library(sf)
library(spatial)
library(spData)
library(spdep)
library(maps)
library(gridExtra)
library(spatstat)
library(deldir)
library(raster)
library(viridis)
library(terra)
library(lubridate)
library(patchwork)
library(fmesher)
library(tidyr)
library(cowplot)

if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable")
}


# ---- Rutas y parámetros globales --------------------------------------------

setwd("/home/jasonromeroia/Documents/personal/Tesis_MCIC/lgcp_inla_spde_sismos/")

files_rds          <- "covariables_rds"
path_image_results <- "imagenes_doc"

source("R/spde-book-functions.R")
source("R/discrete_gradient.R")
source("R/utils.R")

path_file_seismic <- "Data/gdf_espacial_2020.gpkg"


# ---- Priors a partir de la CNN ----------------------------------------------
# scale y varianza estimados por la CNN sobre simulaciones LGCP del catálogo 2020.

scale_hat <- 72346     # scale Matérn de spatstat, en metros
var_hat   <- 6.0905    # varianza del campo latente

# spatstat parametriza con scale y INLA con range. Para nu = 1,
# kappa = sqrt(2)/scale y range = sqrt(8)/kappa = 2*scale. prior.sigma va sobre
# la desviación estándar, no sobre la varianza.
range_cnn <- 3 * scale_hat      # 2*scale (144692 m) con 1.5 de holgura
sigma_cnn <- sqrt(var_hat)      # 2.4679

range_simulated <- range_cnn
sigma_simulated <- sigma_cnn
mu_simulated    <- -21.6429

control_fixed_list <- list(
  mean.intercept = mu_simulated,
  prec.intercept = 0.1
)

# Hiperparámetros del PC-prior que recibe inla.spde2.pcmatern():
#   P(range < range_prior_u) = alpha_pc_range
#   P(sigma > sigma_prior_u) = alpha_pc_sigma
range_prior_u  <- range_simulated
alpha_pc_range <- 0.6
sigma_prior_u  <- sigma_simulated * 1.5
alpha_pc_sigma <- 0.1

save_graphics <- TRUE

dir_out_sp <- file.path(path_image_results, "modelos_INLA_2020")
dir.create(dir_out_sp, recursive = TRUE, showWarnings = FALSE)


# ---- Catálogo sísmico -------------------------------------------------------

shapeZona_sp <- readRDS(paste0(files_rds, "/shapeZona_sp"))
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = TRUE)

sismosSp <- st_read(path_file_seismic)
sismosSp$year <- sismosSp$YEAR

sismosSp <- subset(sismosSp, year >= 2020)
sismosSp <- subset(sismosSp, year <= 2020)

sismosSp$X <- st_coordinates(sismosSp)[, 1]
sismosSp$Y <- st_coordinates(sismosSp)[, 2]

if (any(is.na(sismosSp$X) | is.na(sismosSp$Y))) {
  warning("Hay valores NA en las coordenadas. Eliminando...")
  sismosSp <- sismosSp[!is.na(sismosSp$X) & !is.na(sismosSp$Y), ]
}


# ---- Patrón puntual y malla SPDE --------------------------------------------

creating_ppp <- create_ppp(shapeZona_sp, sismosSp)

p  <- creating_ppp$point_pattern
xy <- creating_ppp$coordinates_eventos

param_cutoff <- 5000
offset_param <- c(100, 20000)

shapeZona <- st_make_valid(shapeZona_sp)
shapeZona <- as_Spatial(shapeZona)
shapeZona <- st_as_sf(shapeZona)

# El mesh se construye sin loc, así que la cuadratura de la integral del LGCP
# son solo los vértices. Dividir por 40 daba celdas de Voronoi de ~20 km, los
# nodos no resolvían las covariables y M2/M3 degeneraban (betas ~ +-200). Con
# /80 las celdas quedan en ~10 km y los cuatro modelos ajustan bien.
max.edge_params <- max(c(diff(range(xy[, 1])), diff(range(xy[, 2])))) / 80

meshSismos <- inla.mesh.2d(
  boundary = shapeZona,
  max.edge = c(1, 3) * max.edge_params,
  cutoff   = param_cutoff,
  offset   = offset_param,
  crs      = st_crs(shapeZona)
)

cat("Número de vértices del mesh:", meshSismos$n)
plot(meshSismos, main = "Malla de Delaunay – Colombia 2020")

p_mesh <- ggplot() +
  geom_fm(data = meshSismos) +
  coord_sf() +
  theme_minimal()

ggsave(paste0(path_image_results, "/Triangulacion.png"),
       plot = p_mesh, width = 14, height = 8, dpi = 300)


# ---- Stack INLA del proceso puntual -----------------------------------------

ptsSismos      <- as.matrix(xy)
mesh.ptsSismos <- as.matrix(meshSismos$loc[, 1:2])
allptsSismos   <- rbind(mesh.ptsSismos, ptsSismos)

nvSismos <- meshSismos$n
nSismos  <- nrow(ptsSismos)

# Malla dual de Voronoi: sus áreas son los pesos de integración.
dmesh <- st_as_sf(book.mesh.dual(meshSismos))

st_crs(dmesh) <- st_crs(shapeZona_sp)
p_dmesh <- ggplot() +
  geom_sf(data = dmesh, fill = NA, color = "black", linewidth = 0.15) +
  geom_sf(data = shapeZona_sp, fill = NA, color = "blue", linewidth = 0.5) +
  theme_minimal() +
  theme(axis.text = element_blank(), panel.grid = element_blank())

ggsave(paste0(path_image_results, "/TriangulacionVoronoi.png"),
       plot = p_dmesh, width = 14, height = 8, dpi = 300)

intersect_idx <- st_intersects(dmesh, shapeZona_sp, sparse = FALSE)[, 1]
dmesh_in      <- dmesh[intersect_idx, ]
intersections <- st_intersection(dmesh_in, shapeZona_sp)

areas <- st_area(intersections)

# w[i] = área del polígono de Voronoi del vértice i dentro de la zona; 0 fuera.
w <- numeric(nrow(dmesh))
w[intersect_idx] <- as.numeric(areas)

total_area  <- as.numeric(st_area(shapeZona_sp))
sum_wSismos <- sum(w)
cat("Área total de Colombia continental:", total_area, "m²\n")
cat("Suma de áreas de Voronoi:", sum_wSismos, "m²\n")

# Formulación Poisson del LGCP: y = 0 con E = área en los vértices,
# y = 1 con E = 0 en los eventos.
y.ppSismos      <- rep(0:1, c(nvSismos, nSismos))
wSismos_numeric <- as.numeric(unlist(w))
e.ppSismos      <- as.numeric(c(w, rep(0, nSismos)))

lmatSismos <- inla.spde.make.A(meshSismos, ptsSismos)
imaSismos  <- Diagonal(nvSismos, rep(1, nvSismos))
A.ppSismos <- rbind(imaSismos, lmatSismos)


# ---- Covariables en vértices y eventos --------------------------------------

expanded_window  <- grow.rectangle(as.rectangle(p$window), 50000)
allpts.pppSismos <- ppp(allptsSismos[, 1], allptsSismos[, 2], expanded_window)

topografia_im_scaled      <- readRDS(paste0(files_rds, "/topografia_im_scaled.rds"))
isostasia_im_scaled       <- readRDS(paste0(files_rds, "/isostasia_im_scaled.rds"))
volcanes_im_scaled        <- readRDS(paste0(files_rds, "/volcanes_im_scaled.rds"))
falla_sinestral_im_scaled <- readRDS(paste0(files_rds, "/sinestral_im_scaled.rds"))
falla_dextral_im_scaled   <- readRDS(paste0(files_rds, "/dextral_im_scaled.rds"))
falla_normal_im_scaled    <- readRDS(paste0(files_rds, "/normal_im_scaled.rds"))
falla_inversa_im_scaled   <- readRDS(paste0(files_rds, "/inversa_im_scaled.rds"))

covar <- list(
  msnm            = topografia_im_scaled,
  isostasia       = isostasia_im_scaled,
  volcanes        = volcanes_im_scaled,
  falla_sinestral = falla_sinestral_im_scaled,
  falla_dextral   = falla_dextral_im_scaled,
  falla_inversa   = falla_inversa_im_scaled,
  falla_normal    = falla_normal_im_scaled
)

# Los nodos de integración muestrean la covariable a la resolución del píxel,
# mientras los eventos caen sobre valores extremos (p.ej. distancia mínima a la
# falla). Esa asimetría genera cuasi-separación y los beta se desbocan. Suavizar
# a la escala de la cuadratura hace que nodos, eventos y celdas de predicción
# vean la misma superficie; el efecto se interpreta a esa resolución.
covar <- lapply(covar, function(X) blur(X, sigma = 1000, bleed = FALSE, normalise = TRUE))

# Indexado matricial en lugar de un sapply punto a punto: evita ~(nv + n)
# llamadas a `[.im` por covariable.
covs100 <- lapply(covar, function(X) {
  px <- nearest.pixel(allpts.pppSismos$x, allpts.pppSismos$y, X)
  X[cbind(px$row, px$col)]
})

covs100$b0 <- rep(1, nvSismos + nSismos)

# Un NA en una covariable vuelve NA el predictor lineal de esa fila.
stopifnot(all(lengths(covs100) == nvSismos + nSismos))
na_ajuste <- sapply(covs100, function(v) sum(is.na(v)))
if (any(na_ajuste > 0)) {
  print(na_ajuste[na_ajuste > 0])
  warning("Covariables con NA en los puntos de ajuste (mesh + eventos).")
} else {
  cat("Covariables en puntos de ajuste: 0 NAs\n")
}


# ---- Ajuste del modelo ------------------------------------------------------
# interval1 es el tamaño de celda (m) de la cuadrícula de predicción,
# use_covariables_model un subconjunto de names(covs100) (c() para el modelo
# nulo) y use_aprior_information activa el PC-prior de la CNN.
# Devuelve el objeto inla, la cuadrícula de predicción y el stack combinado.

modelo_lgcp_col <- function(interval1, use_covariables_model, use_aprior_information) {

  if (use_aprior_information) {
    cat("-------------------------------------\n")
    cat("Modelo con PC-priors informativos\n")
    spdesismos <- inla.spde2.pcmatern(
      mesh        = meshSismos,
      alpha       = 2,                     # nu = 1, d = 2
      prior.range = c(range_prior_u, alpha_pc_range),
      prior.sigma = c(sigma_prior_u, alpha_pc_sigma)
    )
  } else {
    spdesismos <- inla.spde2.matern(meshSismos, alpha = 2)
  }
  n_spde <- spdesismos$n.spde

  spde.indexSismos <- inla.spde.make.index(
    name  = "spatial.field",
    n.spde = n_spde
  )

  spde.stack <- inla.stack(
    data    = list(y = y.ppSismos, e = e.ppSismos),
    A       = list(A.ppSismos, 1),
    effects = list(spde.indexSismos, covs100),
    tag     = "pp"
  )

  # Buffer negativo para no predecir sobre el borde.
  shapeZona_sp_interno <- st_buffer(shapeZona_sp, -1000)
  grid     <- st_make_grid(shapeZona_sp_interno, cellsize = interval1, what = "centers", square = TRUE)
  grid_in  <- grid[st_within(grid, shapeZona_sp_interno, sparse = FALSE)]
  pts.pred <- as.data.frame(st_coordinates(grid_in))

  # expanded_window y no p$window: ppp() descarta en silencio los puntos que
  # caen fuera de la ventana y eso desalinearía el stack.
  ppp.pred <- ppp(pts.pred[, 1], pts.pred[, 2], window = expanded_window)
  stopifnot(npoints(ppp.pred) == nrow(pts.pred))

  A.pred <- inla.spde.make.A(mesh = meshSismos, loc = as.matrix(pts.pred))

  # Mismo método que en el ajuste: X[ppp.pred] devolvería NA fuera del soporte
  # de la imagen y esas celdas se caerían del Log-Score.
  covs100.pred <- lapply(covar, function(X) {
    px <- nearest.pixel(ppp.pred$x, ppp.pred$y, X)
    X[cbind(px$row, px$col)]
  })
  covs100.pred$b0 <- rep(1, nrow(pts.pred))
  stopifnot(all(lengths(covs100.pred) == nrow(pts.pred)))

  # y = NA marca las filas a predecir. e = 0 explícito: si se omite, inla.stack
  # rellena con NA al combinar.
  spde.stack.pred <- inla.stack(
    data    = list(y = rep(NA_real_, nrow(pts.pred)), e = rep(0, nrow(pts.pred))),
    A       = list(A.pred, 1),
    effects = list(spde.indexSismos, covs100.pred),
    tag     = "pred"
  )

  join.stack <- inla.stack(spde.stack, spde.stack.pred)

  formula <- as.formula(
    paste("y ~ 1 +",
          paste(use_covariables_model, collapse = " + "),
          "+ f(spatial.field, model = spdesismos)")
  )

  # Con el prior por defecto (prec = 0.001) la optimización EB diverge en las
  # direcciones de las covariables. Con covariables en [-1,1], prec = 1 es
  # débilmente informativo. Se aplica a los cuatro modelos para no sesgar la
  # comparación M2 vs M3.
  control_fixed_modelo <- list(mean = 0, prec = 1)
  if (use_aprior_information & length(use_covariables_model)>0 ) {
    # Se añade a la lista: asignar control_fixed_list directo pisaría el prior
    # de los beta.
    control_fixed_modelo$mean.intercept <- control_fixed_list$mean.intercept
    control_fixed_modelo$prec.intercept <- control_fixed_list$prec.intercept
  }


  args_inla <- list(
    formula  = formula,
    family   = "poisson",
    data     = inla.stack.data(join.stack),
    # cmin = 0 según los autores de INLA ante los warnings
    # "GMRFLib_2order_approx: NAN/INF in logl" con verosimilitud Poisson.
    control.inla = list(int.strategy = "eb", cmin = 0),
    control.fixed = control_fixed_modelo,
    control.predictor = list(
      A       = inla.stack.A(join.stack),
      compute = TRUE,
      link    = 1
    ),
    control.compute = list(
      config = TRUE,
      cpo    = TRUE
    ),
    E = inla.stack.data(join.stack)$e
  )

  stopifnot(!any(is.na(args_inla$E)), all(args_inla$E >= 0))

  pp.res <- do.call(inla, args_inla)

  return(list(
    result_inla = pp.res,
    grid        = grid_in,
    full_stack  = join.stack
  ))
}


# ---- Ajuste de M0 a M3 ------------------------------------------------------

cat("Lado equivalente de Voronoi medio:", sqrt(mean(w)), "m\n")
interval1 <- 5000

# M2 y M3 comparten covariables por construcción; listarlas dos veces dejaría de
# aislar el efecto del prior.
covariables_modelo <- c("volcanes", "falla_inversa", "falla_normal", "isostasia","msnm")
stopifnot(all(covariables_modelo %in% names(covar)))

# La firma codifica todo lo que define el stack. Si algo cambia, cambia el
# nombre del archivo y el modelo se reajusta en vez de cargarse obsoleto.
firma_cache <- sprintf("i%d_m%d_n%d", interval1, nvSismos, nSismos)
firma_cov   <- paste(sort(covariables_modelo), collapse = "-")

ruta_M0_rds <- file.path(files_rds, sprintf("pp.resM0_2020_%s.rds", firma_cache))
ruta_M2_rds <- file.path(files_rds, sprintf("pp.resM2_2020_%s_%s.rds", firma_cache, firma_cov))
cat("Caché M0:", ruta_M0_rds, "\nCaché M2:", ruta_M2_rds, "\n")

# El tiempo de ajuste se guarda como atributo del .rds para recuperar el de la
# corrida original y no el de la lectura. Cachés antiguos lo dejan en NA.
if (file.exists(ruta_M0_rds)) {
  cat("Cargando modelo M0 (sin prior, sin covariables) desde", ruta_M0_rds, "\n")
  pp.resM0 <- readRDS(ruta_M0_rds)
  segundos_M0 <- attr(pp.resM0, "tiempo_ajuste_s")
  if (is.null(segundos_M0)) segundos_M0 <- NA_real_
} else {
  time_M0 <- system.time({
    pp.resM0 <- modelo_lgcp_col(interval1, c(), FALSE)
  })
  segundos_M0 <- time_M0[["elapsed"]]
  attr(pp.resM0, "tiempo_ajuste_s") <- segundos_M0
  saveRDS(pp.resM0, ruta_M0_rds)
}

time_M1 <- system.time({
  pp.resM1 <- modelo_lgcp_col(interval1, c(), TRUE)
})
segundos_M1 <- time_M1[["elapsed"]]

if (file.exists(ruta_M2_rds)) {
  cat("Cargando modelo M2 (sin prior, con covariables) desde", ruta_M2_rds, "\n")
  pp.resM2 <- readRDS(ruta_M2_rds)
  segundos_M2 <- attr(pp.resM2, "tiempo_ajuste_s")
  if (is.null(segundos_M2)) segundos_M2 <- NA_real_
} else {
  time_M2 <- system.time({
    pp.resM2 <- modelo_lgcp_col(interval1, covariables_modelo, FALSE)
  })
  segundos_M2 <- time_M2[["elapsed"]]
  attr(pp.resM2, "tiempo_ajuste_s") <- segundos_M2
  saveRDS(pp.resM2, ruta_M2_rds)
}

time_M3 <- system.time({
  pp.resM3 <- modelo_lgcp_col(interval1, covariables_modelo, TRUE)
})
segundos_M3 <- time_M3[["elapsed"]]

# Con covariables en [0,1], |beta| > ~10 o sd desmesuradas indican que la
# optimización EB divergió.
cat("\n[diag] Efectos fijos M2:\n"); print(pp.resM2$result_inla$summary.fixed)
cat("\n[diag] Efectos fijos M3:\n"); print(pp.resM3$result_inla$summary.fixed)
idx_diag   <- inla.stack.index(pp.resM3$full_stack, "pred")$data
sd_pred_M3 <- pp.resM3$result_inla$summary.linear.predictor[idx_diag, "sd"]
cat("\n[diag] sd del predictor lineal en las celdas de predicción (M3):\n")
print(summary(sd_pred_M3))
cat("\n[diag] mode.status (0 = Hessiano OK) M2:",
    pp.resM2$result_inla$mode$mode.status,
    "| M3:", pp.resM3$result_inla$mode$mode.status, "\n\n")


# ---- Tiempo de ejecución por modelo -----------------------------------------

df_tiempos <- data.frame(
  Modelo   = factor(c("M0", "M1", "M2", "M3"), levels = c("M0", "M1", "M2", "M3")),
  segundos = c(segundos_M0, segundos_M1, segundos_M2, segundos_M3)
)
if (any(is.na(df_tiempos$segundos))) {
  cat("[aviso] Tiempos NA en:",
      paste(df_tiempos$Modelo[is.na(df_tiempos$segundos)], collapse = ", "),
      "- el caché .rds es anterior al registro de tiempos; borrarlo y reajustar.\n")
}

p_tiempos <- ggplot(df_tiempos, aes(x = Modelo, y = segundos)) +
  geom_col(fill = "steelblue", width = 0.85) +
  labs(x = "Modelo", y = "Tiempo de ejecución (segundos)") +
  theme_bw(base_size = 12)

if (save_graphics) {
  ggsave(file.path(dir_out_sp, "time_execution.png"),
         plot = p_tiempos, width = 4.5, height = 4.5, dpi = 300)
}
print(df_tiempos)


# ---- Resúmenes de efectos fijos ---------------------------------------------

exportar_resumen_fijos <- function(res, modelo, archivo,
                                   solo_intercepto = FALSE, digitos = 4) {
  sf_ <- as.data.frame(res$summary.fixed)
  sf_$Variable <- rownames(sf_)
  if (solo_intercepto) {
    sf_ <- sf_[sf_$Variable %in% c("(Intercept)", "b0", "intercept"), , drop = FALSE]
  }
  fmt <- paste0("%.", digitos, "f")
  tabla <- data.frame(
    Variable     = sf_$Variable,
    Media        = sprintf(fmt, sf_$mean),
    `IC 2.5 %`   = sprintf(fmt, sf_[["0.025quant"]]),
    `IC 97.5 %`  = sprintf(fmt, sf_[["0.975quant"]]),
    check.names  = FALSE
  )
  con <- file(archivo, open = "w", encoding = "UTF-8")
  writeLines(c(paste("Modelo", modelo),
               capture.output(print(tabla, row.names = FALSE))), con)
  close(con)
  cat("Resumen de efectos fijos exportado:", archivo, "\n")
}

exportar_resumen_fijos(pp.resM0$result_inla, "M0",
                       file.path(dir_out_sp, "summary_fijos_M0.txt"),
                       solo_intercepto = TRUE)
exportar_resumen_fijos(pp.resM1$result_inla, "M1",
                       file.path(dir_out_sp, "summary_fijos_M1.txt"),
                       solo_intercepto = TRUE)
exportar_resumen_fijos(pp.resM2$result_inla, "M2",
                       file.path(dir_out_sp, "summary_fijos_M2.txt"))
exportar_resumen_fijos(pp.resM3$result_inla, "M3",
                       file.path(dir_out_sp, "summary_fijos_M3.txt"))

results_model_to_tableLatex(pp.resM0$result_inla)
results_model_to_tableLatex(pp.resM1$result_inla)
results_model_to_tableLatex(pp.resM2$result_inla)
results_model_to_tableLatex(pp.resM3$result_inla)


# ---- Campo latente espacial -------------------------------------------------

if (save_graphics) {
  out_M0_M1 <- plot_spatial_effects_1x2(
    resA          = pp.resM0$result_inla,
    resB          = pp.resM1$result_inla,
    titles        = c("", ""),
    output_path   = file.path(dir_out_sp, "spatial_effects_M0_M1.png"),
    show_stats    = TRUE,
    stats_position    = "topright",
    use_quantile_limits = TRUE,
    q             = c(0.02, 0.98)
  )

  out_M1_M2 <- plot_spatial_effects_1x2(
    resA          = pp.resM2$result_inla,
    resB          = pp.resM3$result_inla,
    titles        = c("", ""),
    output_path   = file.path(dir_out_sp, "spatial_effects_M2_M3.png"),
    show_stats    = TRUE,
    stats_position    = "topright",
    use_quantile_limits = TRUE,
    q             = c(0.02, 0.98)
  )
}

hist(pp.resM3$result_inla$summary.random$spatial.field$mean)


# ---- Intensidades predichas -------------------------------------------------
# Los cuatro modelos comparten cuadrícula de predicción, así que basta con la de M1.

grid_in    <- pp.resM1$grid
join.stack <- pp.resM1$full_stack

pts.pred_ <- as.data.frame(st_coordinates(grid_in))
names(pts.pred_) <- c("x", "y")
pts.pred_$dummy <- 1
coordinates(pts.pred_) <- ~x + y
proj4string(pts.pred_) <- st_crs(shapeZona_sp)$proj4string
gridded(pts.pred_)     <- TRUE

spdf       <- as(pts.pred_, "SpatialPixelsDataFrame")
spdf_poly  <- as(spdf, "SpatialPolygonsDataFrame")
spdf_sf    <- st_as_sf(spdf_poly)
spdf_sf    <- st_transform(spdf_sf, crs = st_crs(sismosSp))

idx <- inla.stack.index(join.stack, "pred")$data

spdf$M0 <- exp(pp.resM0$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M1 <- exp(pp.resM1$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M2 <- exp(pp.resM2$result_inla$summary.linear.predictor[idx, "mean"])
spdf$M3 <- exp(pp.resM3$result_inla$summary.linear.predictor[idx, "mean"])

spdf_sf$M0 <- spdf$M0
spdf_sf$M1 <- spdf$M1
spdf_sf$M2 <- spdf$M2
spdf_sf$M3 <- spdf$M3

stopifnot(st_crs(spdf_sf) == st_crs(sismosSp))

spdf_sf$ID       <- seq_len(nrow(spdf_sf))
spdf_sf$observed <- lengths(st_intersects(spdf_sf, sismosSp))
spdf_sf$area     <- as.numeric(st_area(spdf_sf))

spdf_sf$expect_M0 <- spdf_sf$M0 * spdf_sf$area
spdf_sf$expect_M1 <- spdf_sf$M1 * spdf_sf$area
spdf_sf$expect_M2 <- spdf_sf$M2 * spdf_sf$area
spdf_sf$expect_M3 <- spdf_sf$M3 * spdf_sf$area

# La suma de esperados debe aproximar el total observado; una desviación grande
# delata un ajuste degenerado aunque el mapa se vea razonable.
cat(sprintf("[diag] Eventos observados en grid: %d | Esperados M0: %.0f  M1: %.0f  M2: %.0f  M3: %.0f\n",
            sum(spdf_sf$observed),
            sum(spdf_sf$expect_M0, na.rm = TRUE), sum(spdf_sf$expect_M1, na.rm = TRUE),
            sum(spdf_sf$expect_M2, na.rm = TRUE), sum(spdf_sf$expect_M3, na.rm = TRUE)))


# ---- Residuos de Pearson ----------------------------------------------------

pearson_resid <- function(O, E, eps = 1e-12) {
  ok  <- is.finite(O) & is.finite(E) & E > 0
  r   <- rep(NA_real_, length(O))
  E2  <- pmax(E, eps)
  r[ok] <- (O[ok] - E2[ok]) / sqrt(E2[ok])
  r
}

spdf_sf$pearson_M0 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M0)
spdf_sf$pearson_M1 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M1)
spdf_sf$pearson_M2 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M2)
spdf_sf$pearson_M3 <- pearson_resid(spdf_sf$observed, spdf_sf$expect_M3)

summary(spdf_sf$pearson_M0)
summary(spdf_sf$pearson_M1)
summary(spdf_sf$pearson_M2)
summary(spdf_sf$pearson_M3)

exportar_descriptivos_pearson <- function(residuos, modelo, archivo, digitos = 3) {
  v <- residuos[is.finite(residuos)]
  fmt <- paste0("%.", digitos, "f")
  tabla <- data.frame(
    Estadístico = c("Mínimo", "Q1", "Mediana", "Media", "Q3", "Máximo"),
    Valor = sprintf(fmt, c(min(v), quantile(v, 0.25), median(v), mean(v),
                           quantile(v, 0.75), max(v))),
    check.names = FALSE
  )
  con <- file(archivo, open = "w", encoding = "UTF-8")
  writeLines(c(paste("Modelo", modelo),
               capture.output(print(tabla, row.names = FALSE))), con)
  close(con)
  cat("Descriptivos de residuos exportados:", archivo, "\n")
}

exportar_descriptivos_pearson(spdf_sf$pearson_M0, "M0",
                              file.path(dir_out_sp, "residuos_pearson_M0.txt"))
exportar_descriptivos_pearson(spdf_sf$pearson_M1, "M1",
                              file.path(dir_out_sp, "residuos_pearson_M1.txt"))
exportar_descriptivos_pearson(spdf_sf$pearson_M2, "M2",
                              file.path(dir_out_sp, "residuos_pearson_M2.txt"))
exportar_descriptivos_pearson(spdf_sf$pearson_M3, "M3",
                              file.path(dir_out_sp, "residuos_pearson_M3.txt"))

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


# ---- Log-Score de Poisson ---------------------------------------------------
# LS = sum(log P(O_i | lambda = E_i)); mayor es mejor.

poisson_log_score <- function(observed, expected, eps = 1e-12) {
  ok <- is.finite(observed) & is.finite(expected) & expected > 0
  O  <- observed[ok]
  E  <- pmax(expected[ok], eps)
  sum(dpois(O, lambda = E, log = TRUE))
}

LS_M0 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M0)
LS_M1 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M1)
LS_M2 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M2)
LS_M3 <- poisson_log_score(spdf_sf$observed, spdf_sf$expect_M3)

cat("Log-Scores (mayor es mejor):\n")
print(c(logscore_M0 = LS_M0, logscore_M1 = LS_M1, logscore_M2 = LS_M2, logscore_M3 = LS_M3))


# ---- LCPO -------------------------------------------------------------------
# El stack "pp" tiene primero los nodos de integración (y = 0) y después los
# eventos (y = 1). El LCPO se evalúa solo sobre los eventos: incluir los ~12900
# ceros de cuadratura hace que dominen la métrica.

idx_pp      <- inla.stack.index(pp.resM1$full_stack, "pp")$data
idx_eventos <- idx_pp[(nvSismos + 1):(nvSismos + nSismos)]
stopifnot(length(idx_eventos) == nSismos)

# LCPO = -sum(log(CPO_i)); menor es mejor. Al ser una suma, los modelos deben
# evaluarse sobre los mismos índices, descartando los CPO que INLA marca como
# poco fiables en cualquiera de ellos.
lcpo_scores_comunes <- function(modelos, idx) {
  cpos <- lapply(modelos, function(r) r$cpo$cpo[idx])
  ok   <- rep(TRUE, length(idx))
  for (nm in names(modelos)) {
    fail <- modelos[[nm]]$cpo$failure
    o <- is.finite(cpos[[nm]]) & cpos[[nm]] > 0
    if (!is.null(fail)) o <- o & is.finite(fail[idx]) & fail[idx] == 0
    ok <- ok & o
  }
  cat(sprintf("LCPO sobre eventos comunes y fiables: %d de %d\n", sum(ok), length(ok)))
  vapply(cpos, function(cp) -sum(log(cp[ok])), numeric(1))
}

modelos_inla <- list(
  M0 = pp.resM0$result_inla, M1 = pp.resM1$result_inla,
  M2 = pp.resM2$result_inla, M3 = pp.resM3$result_inla
)
LCPO_todos <- lcpo_scores_comunes(modelos_inla, idx_eventos)
LCPO_M0 <- LCPO_todos[["M0"]]; LCPO_M1 <- LCPO_todos[["M1"]]
LCPO_M2 <- LCPO_todos[["M2"]]; LCPO_M3 <- LCPO_todos[["M3"]]

cat("LCPO (menor es mejor):\n")
print(c(resLCPO_M0 = LCPO_M0, resLCPO_M1 = LCPO_M1, resLCPO_M2 = LCPO_M2, resLCPO_M3 = LCPO_M3))

p1 <- grafico_comparacion_metricas(LS_M0, LS_M1, LS_M2, LS_M3,
                                   LCPO_M0, LCPO_M1, LCPO_M2, LCPO_M3)

ggsave(file.path(dir_out_sp, "log_score_lcpo.png"), plot = p1, width = 10, height = 7)


# ---- Bootstrap M2 vs M3 -----------------------------------------------------
# En cada una de las B remuestras se calcula el Log-Score de ambos modelos sobre
# las mismas celdas; delta > 0 significa que M3 es mejor en esa remuestra.

B <- 5000
delta <- numeric(B)

ok <- is.finite(spdf_sf$observed) & is.finite(spdf_sf$expect_M2) &
  is.finite(spdf_sf$expect_M3) &
  spdf_sf$expect_M2 > 0 & spdf_sf$expect_M3 > 0

O  <- spdf_sf$observed[ok]
E1 <- spdf_sf$expect_M2[ok]
E2 <- spdf_sf$expect_M3[ok]
n  <- length(O)

# La variable del bucle es b_idx y no idx: el idx global (índices de predicción)
# se vuelve a usar más abajo.
for (b in 1:B) {
  b_idx    <- sample.int(n, n, replace = TRUE)
  ls1      <- sum(dpois(O[b_idx], lambda = E1[b_idx], log = TRUE))
  ls2      <- sum(dpois(O[b_idx], lambda = E2[b_idx], log = TRUE))
  delta[b] <- ls2 - ls1
}

cat("Bootstrap (M3 vs M2):\n")
print(c(
  mean_delta   = mean(delta),
  p_M3_better  = mean(delta > 0),
  q025         = quantile(delta, 0.025),
  q975         = quantile(delta, 0.975)
))

delta_bootstrap <- delta

p3 <- grafico_bootstrap(delta_bootstrap)
ggsave(file.path(dir_out_sp, "boostrap.png"), plot = p3$plot, width = 10, height = 7)


# ---- Panel de residuos ------------------------------------------------------

spplot(
  spdf,
  c("pearson_M0", "pearson_M1", "pearson_M2", "pearson_M3"),
  names.attr  = c("M0", "M1", "M2", "M3"),
  col.regions = viridis::plasma(32),
  main        = "Residuos de Pearson por modelo"
)


# ---- Mapa de intensidad estimada (M3) ---------------------------------------

library(grid)
ruta_logo_ud <- file.path(
  "/home/jasonromeroia/Documents/personal/Tesis_MCIC",
  "MCIC_TESIS_2026_JasonRomero/images/logo_ud.png"
)
read_png_rgba <- function(path) {
  x <- png::readPNG(path)
  d <- dim(x)
  if (length(d) == 2) {
    # gris -> RGB
    x <- array(rep(x, 3), dim = c(d[1], d[2], 3))
  } else if (d[3] == 2) {
    # gris + alpha -> RGBA
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

logo <- read_png_rgba("imagenes_doc/logo_ud.png")
logo_grob <- rasterGrob(logo, interpolate = TRUE)

lims_M3 <- as.numeric(quantile(spdf_sf$M3[is.finite(spdf_sf$M3)],
                               probs = c(0.02, 0.98), na.rm = TRUE))

borde_zona <- st_transform(st_geometry(shapeZona_sp), st_crs(spdf_sf))

p_M3_mapa <- ggplot(spdf_sf) +
  geom_sf(aes(fill = M3), color = NA) +
  geom_sf(data = borde_zona, fill = NA, color = "grey25", linewidth = 0.3) +
  scale_fill_viridis_c(
    name   = "Intensidad\nesperada (λ)",
    limits = lims_M3,
    oob    = scales::squish,
    guide  = guide_colourbar(
      barwidth = 1, barheight = 10,
      frame.colour = "black", ticks.colour = "black",
      title.position = "top",
      title.theme = element_text(size = 10, face = "bold", hjust = 0)
    )
  ) +
  coord_sf(expand = FALSE) +
  ggspatial::annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    style       = ggspatial::north_arrow_fancy_orienteering(),
    height      = unit(1.2, "cm"),
    width       = unit(1.2, "cm")
  ) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.30) +
  labs(caption = "") +
  tema_tesis() +
  theme(
    plot.caption = element_text(size = 8, color = "#64748b", hjust = 0)
  )

leyenda_grob <- cowplot::get_legend(
  p_M3_mapa + theme(legend.position = "right")
)
p_mapa_sin_leyenda <- p_M3_mapa + theme(legend.position = "none")

header_map_basic <- paste0(
  "Universidad Distrital\nFrancisco José de Caldas\n",
  "Maestría en Ciencias de la Información\n",
  "y Comunicaciones, Énfasis en Geomática\n",
  "Elaborado por: Jason Romero\n",
  "Sistema de Referencia: EPSG:3116"
)

panel_lateral <- cowplot::ggdraw() +
  cowplot::draw_grob(logo_grob,
                     x = 0.35, y = 0.99, width = 0.55, height = 0.22,
                     hjust = 0.5, vjust = 1) +
  cowplot::draw_text(header_map_basic,
                     x = 0.35, y = 0.74,
                     hjust = 0.5, vjust = 1,
                     size = 8.5, lineheight = 1.25) +
  cowplot::draw_grob(leyenda_grob,
                     x = 0.35, y = 0.65, width = 1, height = 0.45,
                     hjust = 0.5, vjust = 1)

p_M3_final <- cowplot::plot_grid(
  p_mapa_sin_leyenda, panel_lateral,
  ncol = 2, rel_widths = c(0.78, 0.22)
)

ggsave(file.path(dir_out_sp, "intensidades_modelos.png"),
       p_M3_final, width = 10.5, height = 8, dpi = 300, bg = "white")


# ---- Prior vs posterior de los hiperparámetros ------------------------------
# Solo M1 y M3, que son los que usan el PC-prior: la curva gris es su prior y la
# comparación es directa. M0 y M2 usan prior vago y no entran en la figura.
#
# Densidades del PC-prior, iguales a las de inla.spde2.pcmatern():
#   rango (d = 2):  pi(r) = l1 * r^-2 * exp(-l1/r),  l1 = -log(alpha)*r0
#   sigma:          pi(s) = l2 * exp(-l2*s),         l2 = -log(alpha)/s0

lambda1_range <- -log(alpha_pc_range) * range_prior_u
lambda2_sigma <- -log(alpha_pc_sigma) / sigma_prior_u
pc_prior_range <- function(r) lambda1_range * r^(-2) * exp(-lambda1_range / r)
pc_prior_sigma <- function(s) lambda2_sigma * exp(-lambda2_sigma * s)

get_hyper_marg <- function(res, patron) {
  nm <- grep(patron, names(res$marginals.hyperpar), value = TRUE)
  if (length(nm) == 0) return(NULL)
  as.data.frame(res$marginals.hyperpar[[nm[1]]])
}

# Rango en km: la densidad se reescala por el cambio de unidad.
mr_M1 <- get_hyper_marg(pp.resM1$result_inla, "^Range")
mr_M3 <- get_hyper_marg(pp.resM3$result_inla, "^Range")
r_max <- max(mr_M1$x, mr_M3$x, 3 * range_prior_u, na.rm = TRUE)
r_grid <- seq(1, r_max, length.out = 400)
df_range <- rbind(
  data.frame(x = r_grid / 1000,   y = pc_prior_range(r_grid) * 1000, curva = "PC-prior (CNN)"),
  data.frame(x = mr_M1$x / 1000,  y = mr_M1$y * 1000,                curva = "Posterior M1"),
  data.frame(x = mr_M3$x / 1000,  y = mr_M3$y * 1000,                curva = "Posterior M3")
)

# Azul y naranja: par seguro para daltonismo, reforzado con linetype por si se
# imprime en blanco y negro.
col_hp <- c("PC-prior (CNN)" = "#7D8A96", "Posterior M1" = "#E69F00", "Posterior M3" = "#4682B4")
lty_hp <- c("PC-prior (CNN)" = "dashed",  "Posterior M1" = "solid",   "Posterior M3" = "solid")

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

ms_M1 <- get_hyper_marg(pp.resM1$result_inla, "Stdev|Stddev|Sigma")
ms_M3 <- get_hyper_marg(pp.resM3$result_inla, "Stdev|Stddev|Sigma")
s_max <- max(ms_M1$x, ms_M3$x, 1.5 * sigma_prior_u, na.rm = TRUE)
s_grid <- seq(1e-3, s_max, length.out = 400)
df_sigma <- rbind(
  data.frame(x = s_grid,    y = pc_prior_sigma(s_grid), curva = "PC-prior (CNN)"),
  data.frame(x = ms_M1$x,   y = ms_M1$y,                curva = "Posterior M1"),
  data.frame(x = ms_M3$x,   y = ms_M3$y,                curva = "Posterior M3")
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

p_hyper <- (p_range | p_sigma) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))) &
  theme(legend.position = "bottom")

ggsave(file.path(dir_out_sp, "prior_vs_posterior_hyper.png"),
       plot = p_hyper, width = 11, height = 5, dpi = 300)


# ---- Forest plot de efectos fijos -------------------------------------------
# M2 vs M3, sin el intercepto: vive en otra escala y aplastaría los IC.

forest_df <- function(res, modelo) {
  sf_ <- as.data.frame(res$summary.fixed)
  sf_$efecto <- rownames(sf_)
  sf_ <- sf_[!sf_$efecto %in% c("(Intercept)", "b0", "intercept"), , drop = FALSE]
  data.frame(
    efecto = sf_$efecto,
    media  = sf_$mean,
    lo     = sf_[["0.025quant"]],
    hi     = sf_[["0.975quant"]],
    modelo = modelo
  )
}
df_fixed <- rbind(
  forest_df(pp.resM2$result_inla, "M2"),
  forest_df(pp.resM3$result_inla, "M3")
)

p_forest <- ggplot(df_fixed, aes(x = media, y = efecto, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#2d2d2d", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), linewidth = 0.7, fatten = 3) +
  scale_color_manual(values = c("M2" = "#7D8A96", "M3" = "#4682B4")) +
  labs(title = "",
       subtitle = "",
       x = "Coeficiente", y = NULL, color = "Modelo") +
  tema_tesis()

ggsave(file.path(dir_out_sp, "forest_efectos_fijos.png"),
       plot = p_forest, width = 9, height = 5, dpi = 300)


# ---- Correlación Matérn posterior -------------------------------------------
# Traduce el rango a una curva de correlación vs distancia.
# Para nu = 1: rho(h) = (kappa*h) * K_1(kappa*h), con kappa = sqrt(8)/range.

matern_corr_nu1 <- function(h, rng) {
  kappa <- sqrt(8) / rng
  z <- kappa * h
  ifelse(h <= 0, 1, z * besselK(z, 1))
}
hyp_M3  <- pp.resM3$result_inla$summary.hyperpar
rng_i   <- grep("^Range", rownames(hyp_M3))[1]
rng_med <- hyp_M3[rng_i, "0.5quant"]
rng_lo  <- hyp_M3[rng_i, "0.025quant"]
rng_hi  <- hyp_M3[rng_i, "0.975quant"]

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
  labs(title = "Función de correlación Matérn posterior (M3)",
       subtitle = "ν = 1; banda = IC 95% del rango; línea punteada: ρ = 0.1",
       x = "Distancia (km)", y = expression(rho(h))) +
  tema_tesis()

ggsave(file.path(dir_out_sp, "correlacion_matern.png"),
       plot = p_corr, width = 8, height = 5, dpi = 300)


# ---- Probabilidad de excedencia (M3) ----------------------------------------
# El idx global quedó sobrescrito por el bucle bootstrap, así que se recalcula
# desde el mismo stack que se usó para las intensidades.

idx_pred <- inla.stack.index(pp.resM1$full_stack, "pred")$data

lp_mean <- pp.resM3$result_inla$summary.linear.predictor[idx_pred, "mean"]
lp_sd   <- pp.resM3$result_inla$summary.linear.predictor[idx_pred, "sd"]

# Umbral: decil superior de la intensidad de M3.
lambda0 <- as.numeric(quantile(spdf$M3, 0.90, na.rm = TRUE))
spdf$exceed_M3    <- 1 - pnorm((log(lambda0) - lp_mean) / lp_sd)
spdf_sf$exceed_M3 <- spdf$exceed_M3

borde <- st_transform(st_geometry(shapeZona_sp), st_crs(spdf_sf))

p_exceed_mapa <- ggplot(spdf_sf) +
  geom_sf(aes(fill = exceed_M3), color = NA) +
  geom_sf(data = borde, fill = NA, color = "grey25", linewidth = 0.3) +
  scale_fill_viridis_c(
    option    = "rocket",
    direction = -1,
    name      = "Probabilidad de\nexcedencia (P)",
    limits    = c(0, 1),
    guide     = guide_colourbar(
      barwidth = 1, barheight = 10,
      frame.colour = "black", ticks.colour = "black",
      title.position = "top",
      title.theme = element_text(size = 10, face = "bold", hjust = 0)
    )
  ) +
  coord_sf(expand = FALSE) +
  ggspatial::annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    style       = ggspatial::north_arrow_fancy_orienteering(),
    height      = unit(1.2, "cm"),
    width       = unit(1.2, "cm")
  ) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.30) +
  labs(caption = "") +
  tema_tesis() +
  theme(
    plot.caption = element_text(size = 8, color = "#64748b", hjust = 0)
  )

leyenda_exceed_grob <- cowplot::get_legend(
  p_exceed_mapa + theme(legend.position = "right")
)
p_exceed_sin_leyenda <- p_exceed_mapa + theme(legend.position = "none")

panel_lateral_exceed <- cowplot::ggdraw() +
  cowplot::draw_grob(logo_grob,
                     x = 0.35, y = 0.99, width = 0.55, height = 0.22,
                     hjust = 0.5, vjust = 1) +
  cowplot::draw_text(header_map_basic,
                     x = 0.35, y = 0.74,
                     hjust = 0.5, vjust = 1,
                     size = 8.5, lineheight = 1.25) +
  cowplot::draw_grob(leyenda_exceed_grob,
                     x = 0.35, y = 0.65, width = 1, height = 0.45,
                     hjust = 0.5, vjust = 1)

p_exceed_final <- cowplot::plot_grid(
  p_exceed_sin_leyenda, panel_lateral_exceed,
  ncol = 2, rel_widths = c(0.78, 0.22)
)

ggsave(file.path(dir_out_sp, "excedencia_M3.png"),
       p_exceed_final, width = 10.5, height = 8, dpi = 300, bg = "white")
