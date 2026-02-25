################################################################################
# SCRIPT: script_rds.R
# PROPÓSITO: Preprocesamiento y serialización de covariables del modelo LGCP
#
# DESCRIPCIÓN:
#   Carga las fuentes de datos originales (rasters y shapefiles), los reproyecta
#   al CRS objetivo (EPSG:3116 – Magna Sirgas Colombia Bogotá), calcula mapas de
#   distancia para las estructuras geológicas lineales (fallas) y puntuales
#   (volcanes), escala cada covariable al rango [-1, 1] y guarda los objetos
#   resultantes como archivos .rds para ser consumidos por los modelos INLA-SPDE.
#
# ENTRADAS REQUERIDAS:
#   - data_new/clip_zona_continental.geojson         : polígono de zona de estudio
#   - data/raster/Topography_Colombia_2000m.tif      : DEM (elevación en metros)
#   - data/raster/mosaico_isostasia_cor.tif          : anomalía isostática (mGal)
#   - data/shapefile/volcanes/volcanes_61.shp         : puntos de volcanes activos
#   - data/shapefile/fallas/falla_rumbo_sinestral.shp : fallas sinestrales
#   - data/shapefile/fallas/falla_rumbo_dextral.shp   : fallas dextrales
#   - data/shapefile/fallas/falla_normal.shp          : fallas normales
#   - data/shapefile/fallas/falla_inversa.shp         : fallas inversas
#
# SALIDAS (en covariables_rds/):
#   - shapeZona_sp               : zona de estudio como objeto sf en EPSG:3116
#   - topografia_im.rds          : imagen spatstat de elevación (sin escalar)
#   - topografia_im_scaled.rds   : elevación escalada a [-1, 1]
#   - isostasia_im.rds           : imagen spatstat de isostasia (sin escalar)
#   - isostasia_im_scaled.rds    : isostasia escalada a [-1, 1]
#   - volcanes_im.rds            : distancia a volcanes (sin escalar)
#   - volcanes_im_scaled.rds     : distancia a volcanes escalada a [-1, 1]
#   - sinestral_im.rds / _scaled : distancia a fallas sinestrales
#   - dextral_im.rds / _scaled   : distancia a fallas dextrales
#   - normal_im.rds / _scaled    : distancia a fallas normales
#   - inversa_im.rds / _scaled   : distancia a fallas inversas
#
# NOTA: Ajustar setwd() a la ruta local del proyecto antes de ejecutar.
################################################################################

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
library("deldir")
library(terra)
library(lubridate)
library(fmesher)

if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable")
}

# CONFIGURACIÓN: cambiar esta ruta al directorio raíz del proyecto local
setwd("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/")
outpath_rds = "covariables_rds"
source("R/spde-book-functions.R")

target_crs = 3116

################################################################################
# FUNCIONES AUXILIARES
################################################################################

# reproj_sf_verbose: Reproyecta un objeto sf a EPSG:3116 imprimiendo mensajes
# Parámetros:
#   obj        : objeto sf a reproyectar
#   target_crs : código EPSG destino (por defecto 3116)
# Retorna: objeto sf reproyectado (invisiblemente)
reproj_sf_verbose <- function(obj, target_crs = 3116) {
  obj_name <- deparse(substitute(obj))
  crs_now <- st_crs(obj)
  msg_now <- if (is.na(crs_now)) "sin CRS" else crs_now$input
  message(obj_name, " → ", msg_now)
  if (is.na(crs_now) || crs_now$epsg != target_crs) {
    obj <- st_transform(obj, target_crs)
    message("   reproyectado a EPSG:", target_crs)
  }
  invisible(obj)
}

# reproj_rast: Reproyecta un raster (Raster* o SpatRaster) a EPSG objetivo.
# Acepta tanto objetos {raster} como {terra}. Devuelve en la misma clase.
# Parámetros:
#   r            : objeto raster a reproyectar
#   target_epsg  : código EPSG destino (por defecto 3116)
#   method       : método de interpolación ("bilinear" para continuas)
reproj_rast <- function(r, target_epsg = 3116, method = "bilinear") {
  is_raster <- inherits(r, "Raster")
  if (is_raster) r <- terra::rast(r)
  tgt_crs <- paste0("EPSG:", target_epsg)
  needs_proj <- !terra::same.crs(r, tgt_crs)
  if (needs_proj) {
    message("   reproyectando a ", tgt_crs, " …")
    r <- terra::project(r, tgt_crs, method = method)
  } else {
    message("   CRS ya es ", tgt_crs)
  }
  if (is_raster) r <- raster::brick(r)
  return(r)
}

# convertRasterToImg: Convierte un objeto Raster* al formato im de {spatstat}.
# Realiza la transposición y el volteo necesarios para alinear filas/columnas
# entre la representación raster (origen arriba-izquierda) y la de spatstat
# (origen abajo-izquierda).
# Parámetros:
#   input_raster : objeto de clase Raster*
# Retorna: objeto im de spatstat
convertRasterToImg = function(input_raster){
  matriz <- as.matrix(input_raster)
  ext <- extent(input_raster)
  xrange <- c(ext@xmin, ext@xmax)
  yrange <- c(ext@ymin, ext@ymax)
  ventana <- owin(xrange = xrange, yrange = yrange)
  # Transposición + inversión de filas para corregir orientación
  matriz_ajustada_ <- apply(t(matriz), 1, rev)
  raster_to_im <- as.im(matriz_ajustada_, W = ventana)
  return(raster_to_im)
}

# scale_im_to_minus1_1: Escala los valores de un objeto im al rango [-1, 1]
# usando el mínimo y máximo globales (sin recorte de percentiles).
# Parámetros:
#   im_obj : objeto im de spatstat
# Retorna: lista con $im (objeto im escalado)
scale_im_to_minus1_1 <- function(im_obj) {
  values <- im_obj$v
  min_val <- min(values, na.rm = TRUE)
  max_val <- max(values, na.rm = TRUE)
  scaled_values <- 2 * ((values - min_val) / (max_val - min_val)) - 1
  im_obj$v <- scaled_values
  return(list(im = im_obj))
}

# scale_im_to_unit_range: Escala los valores de un im a [-1, 1] aplicando
# primero un recorte por percentiles (lower_q, upper_q) para reducir la
# influencia de valores extremos (outliers). FUNCIÓN USADA EN EL MODELO.
# Parámetros:
#   im_obj  : objeto im de spatstat
#   lower_q : cuantil inferior para recorte (por defecto 2%)
#   upper_q : cuantil superior para recorte (por defecto 98%)
# Retorna: lista con $im, $min, $max, $q_low, $q_high
scale_im_to_unit_range <- function(im_obj, lower_q = 0.02, upper_q = 0.98) {
  if (!("v" %in% names(im_obj))) {
    stop("El objeto 'im_obj' debe contener un campo 'v'")
  }
  values <- im_obj$v
  if (all(is.na(values))) stop("Los valores en 'im_obj$v' son todos NA")
  vals_ok <- values[!is.na(values)]
  if (length(unique(vals_ok)) <= 1) {
    stop("Los valores son constantes — no se puede escalar a [-1,1]")
  }
  q_low  <- quantile(vals_ok, probs = lower_q, na.rm = TRUE)
  q_high <- quantile(vals_ok, probs = upper_q, na.rm = TRUE)
  values_clipped <- pmin(pmax(values, q_low), q_high)
  min_val <- min(values_clipped, na.rm = TRUE)
  max_val <- max(values_clipped, na.rm = TRUE)
  values_01  <- (values_clipped - min_val) / (max_val - min_val)
  values_m11 <- values_01 * 2 - 1
  im_obj$v <- values_m11
  return(list(im = im_obj, min = min_val, max = max_val,
              q_low = q_low, q_high = q_high))
}

# scale_im_to_standard: Estandariza los valores de un im a media 0 y sd 1.
# (No se usa en el modelo final; se deja como alternativa de escalado.)
scale_im_to_standard <- function(im_obj) {
  if (!exists("v", im_obj)) stop("El objeto 'im_obj' debe contener un campo 'v'")
  values <- im_obj$v
  if (all(is.na(values))) stop("Los valores en 'im_obj$v' son todos NA")
  if (length(unique(values[!is.na(values)])) <= 1) {
    stop("Los valores en 'im_obj$v' son constantes, no se puede estandarizar")
  }
  mean_val <- mean(values, na.rm = TRUE)
  sd_val   <- sd(values, na.rm = TRUE)
  im_obj$v <- (values - mean_val) / sd_val
  return(list(im = im_obj, mean = mean_val, sd = sd_val))
}

# sf_to_psp: Convierte un objeto sf de líneas (LINESTRING o MULTILINESTRING)
# al formato psp (planar segment pattern) de {spatstat}.
# Parámetros:
#   sf_lines    : objeto sf con geometría de líneas
#   win         : ventana owin donde se define el patrón
#   epsg_target : CRS destino para reproyección (por defecto 3116)
# Retorna: objeto psp de spatstat
sf_to_psp <- function(sf_lines, win, epsg_target = 3116) {
  sf_proj <- st_transform(sf_lines, crs = epsg_target)
  lines_list <- st_geometry(sf_proj)
  x0 <- y0 <- x1 <- y1 <- c()
  for (line in lines_list) {
    coords <- st_coordinates(line)
    if ("L2" %in% colnames(coords)) {
      line_ids <- unique(coords[, "L2"])
      for (id in line_ids) {
        sub_coords <- coords[coords[, "L2"] == id, ]
        if (nrow(sub_coords) < 2) next
        for (i in 1:(nrow(sub_coords) - 1)) {
          x0 <- c(x0, sub_coords[i, 1]);     y0 <- c(y0, sub_coords[i, 2])
          x1 <- c(x1, sub_coords[i+1, 1]);   y1 <- c(y1, sub_coords[i+1, 2])
        }
      }
    } else {
      if (nrow(coords) < 2) next
      for (i in 1:(nrow(coords) - 1)) {
        x0 <- c(x0, coords[i, 1]);   y0 <- c(y0, coords[i, 2])
        x1 <- c(x1, coords[i+1, 1]); y1 <- c(y1, coords[i+1, 2])
      }
    }
  }
  return(psp(x0, y0, x1, y1, window = win))
}

################################################################################
# CARGA Y GUARDADO DE LA ZONA DE ESTUDIO
################################################################################

# shapeZona_sp: polígono de la zona de estudio continental de Colombia
# (EPSG:3116 – Magna Sirgas Colombia Bogotá), usado como dominio D del modelo
shapeZona <- st_read("data_new/clip_zona_continental.geojson")
shapeZona <- st_make_valid(shapeZona)
shapeZona_sp <- as_Spatial(shapeZona)
shapeZona_sp <- st_as_sf(shapeZona_sp)

saveRDS(shapeZona_sp, paste0(outpath_rds, "/shapeZona_sp"))

################################################################################
# CARGA DE DATOS FUENTE
################################################################################

# Rasters continuos
topography <- raster("data/raster/Topography_Colombia_2000m.tif")  # DEM (m)
isostatic  <- raster("data/raster/mosaico_isostasia_cor.tif")       # isostasia (mGal)

# Estructuras geológicas (vectoriales)
volcanes              <- st_read("data/shapefile/volcanes/volcanes_61.shp")
falla_rumbo_sinestral <- st_read("data/shapefile/fallas/falla_rumbo_sinestral.shp")
falla_rumbo_dextral   <- st_read("data/shapefile/fallas/falla_rumbo_dextral.shp")
falla_normal          <- st_read("data/shapefile/fallas/falla_normal.shp")
falla_inversa         <- st_read("data/shapefile/fallas/falla_inversa.shp")

################################################################################
# REPROYECCIÓN A EPSG:3116
################################################################################

for (nm in c("volcanes", "falla_rumbo_sinestral", "falla_rumbo_dextral",
             "falla_normal", "falla_inversa")) {
  assign(nm, reproj_sf_verbose(get(nm)))
}
topography <- reproj_rast(topography)
isostatic  <- reproj_rast(isostatic)

################################################################################
# RECORTE Y ALINEACIÓN DE RASTERS A LA ZONA DE ESTUDIO
################################################################################

# Buffer de 1 km para evitar efectos de borde en el recorte
shapeZona_sp_buffer <- st_buffer(shapeZona_sp, 1000)

# Elev2: raster de topografía recortado al buffer
Elev2 <- crop(topography, shapeZona_sp_buffer)

# Isost: isostasia recortada y remuestreada a la misma resolución que Elev2
Isost                <- crop(isostatic, shapeZona_sp_buffer)
isostatic_resampled  <- resample(isostatic, Elev2, method = "bilinear")
isostatic_aligned    <- mask(isostatic_resampled, Elev2)

################################################################################
# CÁLCULO DE MAPAS DE DISTANCIA PARA FALLAS Y VOLCANES (via {terra})
################################################################################

# r: raster template vacío sobre el buffer, resolución 2000 m
r <- rast(shapeZona_sp_buffer, resolution = 2000)

# fault_vect: lista de fallas convertidas a SpatVector de terra
fault_vect <- list(
  sinestral = vect(falla_rumbo_sinestral),
  dextral   = vect(falla_rumbo_dextral),
  normal    = vect(falla_normal),
  inversa   = vect(falla_inversa)
)

# fault_dist: lista de rasters de distancia mínima a cada tipo de falla (metros)
fault_dist <- purrr::map(fault_vect, ~ terra::distance(r, .x))

# volcanes_distmap: raster de distancia mínima a volcanes activos (metros)
volcanes_vect    <- vect(volcanes)
volcanes_distmap <- distance(r, volcanes_vect)

################################################################################
# FUNCIÓN prep_and_save: ESCALA UNA COVARIABLE Y GUARDA AMBAS VERSIONES (.rds)
################################################################################

# prep_and_save: Convierte un raster a im de spatstat, aplica escalado [-1,1]
# por percentiles y guarda dos archivos: la imagen cruda (_im.rds) y la escalada
# (_im_scaled.rds). Ambas versiones se guardan para permitir interpretación
# posterior de los coeficientes del modelo.
# Parámetros:
#   rast     : objeto SpatRaster o Raster* de entrada
#   name     : nombre base para los archivos de salida
#   path_out : directorio de destino
prep_and_save <- function(rast, name, path_out) {
  if (inherits(rast, "SpatRaster"))
    rast <- raster::raster(rast)
  im_obj        <- convertRasterToImg(rast)
  result_scaled <- scale_im_to_unit_range(im_obj)   # escalado [-1, 1] con recorte
  result_scaled <- result_scaled$im
  if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)
  saveRDS(result_scaled, file = file.path(path_out, paste0(name, "_im_scaled.rds")))
  saveRDS(im_obj,        file = file.path(path_out, paste0(name, "_im.rds")))
}

################################################################################
# EJECUCIÓN: ESCALAR Y GUARDAR TODAS LAS COVARIABLES
################################################################################

# all_covs: lista unificada de todas las covariables del modelo
all_covs <- c(
  fault_dist,                        # sinestral, dextral, normal, inversa
  topografia = Elev2,                # elevación (MSNM)
  isostasia  = isostatic_resampled,  # anomalía isostática
  volcanes   = volcanes_distmap      # distancia a volcanes
)

# Escala y serializa todas las covariables en un solo paso
cov_imgs_scaled <- purrr::imap(all_covs, prep_and_save, path_out = outpath_rds)
