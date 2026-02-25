################################################################################
# SCRIPT: simulations_rglcp.R
# PROPÓSITO: Simulación de procesos Cox Log-Gaussianos (LGCP) para entrenamiento
#            de una red neuronal CNN que estime los parámetros a priori del modelo
#
# DESCRIPCIÓN:
#   Implementa la metodología de Verönneau-Iphigenie et al. (2022), extendida con
#   extracción de características funcionales adicionales de la función L̂(r) de
#   Besag. Se generan 10,000 realizaciones LGCP con parámetros (μ, σ², scale)
#   muestreados aleatoriamente, se calcula L̂(r) para cada patrón puntual, y se
#   extraen 12 características que sirven como entradas a la CNN.
#
# METODOLOGÍA:
#   1. Muestrear parámetros (μ, σ², scale) de distribuciones uniformes
#   2. Simular un proceso LGCP con kernel Matérn (ν=1) sobre la ventana de Colombia
#   3. Calcular la función L centrada: L̂(r) - r (corrección de borde "border")
#   4. Extraer 12 características de la curva L̂(r)
#   5. Guardar los resultados por bloques (chunks) de 100 simulaciones
#
# ENTRADAS REQUERIDAS:
#   - data_new/clip_zona_continental_simplificado.geojson : zona de estudio
#     (versión simplificada para acelerar el cálculo de owin)
#
# SALIDAS:
#   - results_lgcp_features_2/Data_LGCP_batch_XXXX.rds : archivos por chunk
#   - Cada chunk contiene un tibble con columnas:
#       mu, var, scale, N          : parámetros del proceso y número de puntos
#       r, L                       : vectores numéricos (listas) de la curva L̂(r)
#       L_max, L_mean, L_var, ...  : características funcionales (ver abajo)
#
# NOTA: Ejecutar con Rscript o en sesión interactiva; usa 10 cores por defecto.
#       Ajustar setwd() y ncores según el entorno local.
################################################################################

library(pbmcapply)   # paralelización con barra de progreso
library(tidyverse)
library(spatstat)    # análisis de procesos puntuales
library(sf)

set.seed(123)
ncores <- 10

# CONFIGURACIÓN: ajustar esta ruta al directorio raíz del proyecto
setwd("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/")

################################################################################
# CARGA DE LA VENTANA DE SIMULACIÓN
################################################################################

# shapeZona_sp: polígono simplificado de Colombia continental (EPSG:3116)
# Se simplifica agresivamente (dTolerance=40km) para agilizar las 10,000 sims
shapeZona_sp <- st_read("data_new/clip_zona_continental_simplificado.geojson")
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 40000, preserveTopology = TRUE)

# r_iso_owin: ventana de observación en formato owin de spatstat
# Define el dominio D sobre el que se simulan los procesos puntuales
r_iso_owin <- as.owin(shapeZona_sp)

# area_win: área total de la ventana en m² (referencia para calcular μ)
area_win <- area.owin(r_iso_owin)
cat("Área de la ventana:", area_win, "m²\n")

################################################################################
# DEFINICIÓN DE RANGOS DE PARÁMETROS
################################################################################

ntrain <- 10000   # número total de simulaciones objetivo

# N_min, N_max: rango de número esperado de puntos por patrón
# Se usa escala log-uniforme para cubrir órdenes de magnitud variados
N_min <- 500
N_max <- 80000

log_N_expected <- runif(ntrain, log(N_min), log(N_max))
expected_N     <- exp(log_N_expected)

# var: varianza del campo Gaussiano latente σ² — controla el grado de clustering
# Rango [0.3, 3.0]: valores realistas para patrones sísmicos agregados
var <- runif(ntrain, 0.3, 3.0)

# scale: escala de correlación espacial del kernel Matérn (en metros)
# Rango [80km, 200km]: basado en evidencia empírica de la sismicidad colombiana
scale <- runif(ntrain, 80000, 200000)

# mu: intensidad media del campo Gaussiano latente en escala logarítmica
# Calculado analíticamente para obtener el número esperado de puntos target:
#   E[N] = exp(μ + σ²/2) · área  →  μ = log(E[N]/área) - σ²/2
mu <- log(expected_N / area_win) - var / 2

cat("\n===== RANGOS DE PARÁMETROS =====\n")
cat("mu:", round(range(mu), 2), "\n")
cat("var:", round(range(var), 2), "\n")
cat("scale (km):", round(range(scale)/1000, 1), "\n")
cat("E[N]:", round(range(expected_N)), "\n")

################################################################################
# FUNCIÓN extract_L_features
# Extrae 12 características de la curva L̂(r) centrada
################################################################################

# extract_L_features: Calcula estadísticos funcionales de L̂(r) que permiten
# romper parcialmente la no-identificabilidad entre var y scale del LGCP.
# Parámetros:
#   r : vector de distancias (radios de evaluación de L)
#   L : vector de L̂(r) centrada (L_border - r)
#   N : número de puntos del patrón
# Retorna: lista nombrada con 12 características
extract_L_features <- function(r, L, N) {

  idx_pos <- r > 0
  r_pos   <- r[idx_pos]
  L_pos   <- L[idx_pos]

  # 1. L_max: valor máximo de L̂(r) — intensidad del clustering
  L_max <- max(L, na.rm = TRUE)

  # 2. L_mean: media de L̂(r) — nivel promedio de agregación
  L_mean <- mean(L, na.rm = TRUE)

  # 3. L_var: varianza de L̂(r) — heterogeneidad del clustering a distintas escalas
  L_var <- var(L, na.rm = TRUE)

  # 4. L_min: valor mínimo de L̂(r)
  L_min <- min(L, na.rm = TRUE)

  # 5. r_at_Lmax: radio donde L̂(r) es máximo — directamente relacionado con scale
  r_at_Lmax <- r[which.max(L)]

  # 6. AUC: área bajo la curva L̂(r) — clustering total integrado
  AUC <- sum(diff(r) * (head(L, -1) + tail(L, -1)) / 2, na.rm = TRUE)

  # 7. slope_init: pendiente de la regresión lineal en los primeros 10 radios
  #    Captura el clustering a corta distancia
  n_init <- min(10, length(r_pos))
  if (n_init >= 3) {
    slope_init <- tryCatch(
      coef(lm(L_pos[1:n_init] ~ r_pos[1:n_init]))[2],
      error = function(e) 0
    )
  } else {
    slope_init <- 0
  }

  # 8. slope_mid: pendiente en la zona media [30%, 70%] — decaimiento del clustering
  n_mid      <- length(r_pos)
  mid_start  <- floor(n_mid * 0.3)
  mid_end    <- floor(n_mid * 0.7)
  if (mid_end > mid_start + 2) {
    slope_mid <- tryCatch(
      coef(lm(L_pos[mid_start:mid_end] ~ r_pos[mid_start:mid_end]))[2],
      error = function(e) 0
    )
  } else {
    slope_mid <- 0
  }

  # 9-10. curvature_mean, curvature_max: segunda derivada aproximada de L̂(r)
  #       suavizada con ventana de 5 puntos — forma de la curva
  if (length(L) > 4) {
    L_smooth <- stats::filter(L, rep(1/5, 5), sides = 2)
    L_smooth <- L_smooth[!is.na(L_smooth)]
    if (length(L_smooth) > 2) {
      d2L            <- diff(diff(L_smooth))
      curvature_mean <- mean(d2L, na.rm = TRUE)
      curvature_max  <- max(abs(d2L), na.rm = TRUE)
    } else {
      curvature_mean <- 0; curvature_max <- 0
    }
  } else {
    curvature_mean <- 0; curvature_max <- 0
  }

  # 11. peak_ratio: L_max / r_at_Lmax — relación entre amplitud y alcance del pico
  peak_ratio <- ifelse(r_at_Lmax > 0, L_max / r_at_Lmax, 0)

  # 12. L_skew: asimetría de la distribución de L̂(r)
  L_sd <- sd(L, na.rm = TRUE)
  if (L_sd > 0) {
    L_skew <- mean((L - mean(L, na.rm = TRUE))^3, na.rm = TRUE) / (L_sd^3)
  } else {
    L_skew <- 0
  }

  list(L_max = L_max, L_mean = L_mean, L_var = L_var, L_min = L_min,
       r_at_Lmax = r_at_Lmax, AUC = AUC,
       slope_init = as.numeric(slope_init), slope_mid = as.numeric(slope_mid),
       curvature_mean = curvature_mean, curvature_max = curvature_max,
       peak_ratio = peak_ratio, L_skew = L_skew)
}

################################################################################
# FUNCIÓN run_one_sim
# Simula un proceso LGCP y extrae sus características funcionales
################################################################################

# run_one_sim: Genera una realización LGCP con parámetros dados, calcula la
# función L centrada y extrae las 12 características. Retorna NULL si la
# simulación falla o el número de puntos está fuera del rango permitido.
# Parámetros:
#   mu, var, scale  : parámetros del proceso LGCP
#   win             : ventana owin de simulación
#   min_points      : mínimo de puntos aceptables (default: 30)
#   max_points      : máximo de puntos aceptables (default: 150,000)
# Retorna: lista con parámetros, curva L, características; o NULL si inválido
run_one_sim <- function(mu, var, scale, win, min_points = 30, max_points = 150000) {

  # Simular proceso LGCP con kernel Matérn (ν=1), resolución 128×128
  pp <- tryCatch({
    spatstat.random::rLGCP(
      model = "matern", nu = 1,
      mu = mu, var = var, scale = scale,
      win = win, dimyx = c(128, 128), saveLambda = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(pp)) return(NULL)

  N <- spatstat.geom::npoints(pp)
  if (N < min_points || N > max_points) return(NULL)

  # Calcular L̂(r) con corrección de borde, hasta 200 km, en 128 radios
  rmax  <- 200000
  L_obj <- tryCatch({
    spatstat.explore::Lest(pp, correction = "border", rmax = rmax, nrval = 128)
  }, error = function(e) NULL)

  if (is.null(L_obj)) return(NULL)

  # L centrada: L̂_border(r) - r
  # Bajo CSR la L centrada es ≈ 0; valores positivos indican clustering
  L_centered <- L_obj$border - L_obj$r
  if (any(is.na(L_centered))) return(NULL)

  # Extraer las 12 características funcionales
  features <- extract_L_features(L_obj$r, L_centered, N)

  list(mu = mu, var = var, scale = scale, N = N,
       r = L_obj$r, L = L_centered,
       L_max = features$L_max, L_mean = features$L_mean,
       L_var = features$L_var, L_min = features$L_min,
       r_at_Lmax = features$r_at_Lmax, AUC = features$AUC,
       slope_init = features$slope_init, slope_mid = features$slope_mid,
       curvature_mean = features$curvature_mean, curvature_max = features$curvature_max,
       peak_ratio = features$peak_ratio, L_skew = features$L_skew)
}

################################################################################
# SIMULACIÓN POR BLOQUES (CHUNKS)
# Se procesa en chunks de 100 para gestionar memoria y permitir recuperación
################################################################################

chunk_size <- 100
n_chunks   <- ceiling(ntrain / chunk_size)
out_dir    <- "results_lgcp_features_2"

if (!dir.exists(out_dir)) dir.create(out_dir)

total_valid <- 0

for (k in seq_len(n_chunks)) {
  cat("Chunk", k, "of", n_chunks, "\n")

  idx_start <- (k - 1) * chunk_size + 1
  idx_end   <- min(k * chunk_size, ntrain)
  idx       <- idx_start:idx_end

  mu_k    <- mu[idx]
  var_k   <- var[idx]
  scale_k <- scale[idx]

  # Paralelizar con pbmcmapply (muestra progreso en consola)
  sims_k <- pbmcmapply(
    FUN = function(mu_i, var_i, scale_i) {
      run_one_sim(mu_i, var_i, scale_i, r_iso_owin)
    },
    mu_k, var_k, scale_k,
    SIMPLIFY = FALSE,
    mc.cores = ncores
  )

  # Filtrar simulaciones inválidas (NULL)
  sims_k <- sims_k[!sapply(sims_k, is.null)]

  if (length(sims_k) == 0) {
    cat("  Warning: No valid simulations in chunk", k, "\n")
    next
  }

  # Construir tibble del chunk con parámetros + curva L + características
  Data_LGCP_k <- purrr::map_dfr(sims_k, ~ tibble(
    mu = .x$mu, var = .x$var, scale = .x$scale, N = .x$N,
    r  = list(.x$r), L = list(.x$L),
    # Características funcionales de L̂(r)
    L_max = .x$L_max, L_mean = .x$L_mean, L_var = .x$L_var, L_min = .x$L_min,
    r_at_Lmax = .x$r_at_Lmax, AUC = .x$AUC,
    slope_init = .x$slope_init, slope_mid = .x$slope_mid,
    curvature_mean = .x$curvature_mean, curvature_max = .x$curvature_max,
    peak_ratio = .x$peak_ratio, L_skew = .x$L_skew
  ))

  total_valid <- total_valid + nrow(Data_LGCP_k)
  cat("  Valid simulations:", nrow(Data_LGCP_k), "/", length(idx),
      "| Total acumulado:", total_valid, "\n")

  file_k <- file.path(out_dir, sprintf("Data_LGCP_batch_%04d.rds", k))
  saveRDS(Data_LGCP_k, file_k)
  gc()   # liberar memoria entre chunks
}

cat("\n===== SIMULACIÓN COMPLETADA =====\n")
cat("Total de simulaciones válidas:", total_valid, "\n")

################################################################################
# VERIFICACIÓN POST-SIMULACIÓN
# Carga todos los chunks, verifica distribuciones y calcula correlaciones
################################################################################

files     <- list.files(out_dir, pattern = "Data_LGCP_batch_.*\\.rds$", full.names = TRUE)
Data_LGCP <- purrr::map_dfr(files, readRDS)

cat("\nDistribución de N en las simulaciones:\n")
print(summary(Data_LGCP$N))
cat("\nColumnas disponibles:\n")
print(names(Data_LGCP))

# Histograma de N: verifica cobertura del rango objetivo
library(ggplot2)
p_N <- ggplot(Data_LGCP, aes(x = N)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 63000, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 63000, y = Inf, label = "Tu N real (~63k)",
           vjust = 2, hjust = -0.1, color = "red") +
  scale_x_log10() +
  labs(x = "Número de puntos (N)", y = "Frecuencia") +
  theme_minimal()
print(p_N)

# Cobertura del espacio de parámetros (var vs scale)
p_params <- Data_LGCP %>%
  ggplot() +
  geom_point(aes(x = var, y = scale/1000), alpha = 0.3) +
  labs(x = "var", y = "scale (km)", title = "Cobertura de parámetros") +
  theme_minimal()
print(p_params)

################################################################################
# CORRELACIONES CARACTERÍSTICAS vs PARÁMETROS
# Identifica qué características son más informativas para estimar var y scale
################################################################################

cat("\n===== CORRELACIONES FEATURES vs PARÁMETROS =====\n")
feature_cols <- c("L_max", "L_mean", "L_var", "L_min", "r_at_Lmax", "AUC",
                  "slope_init", "slope_mid", "curvature_mean", "curvature_max",
                  "peak_ratio", "L_skew")

cor_matrix <- cor(Data_LGCP[, c("var", "scale", feature_cols)],
                  use = "pairwise.complete.obs")

cat("\nCorrelación con var:\n")
print(round(sort(cor_matrix["var", feature_cols], decreasing = TRUE), 3))

cat("\nCorrelación con scale:\n")
print(round(sort(cor_matrix["scale", feature_cols], decreasing = TRUE), 3))

# Gráficos de dispersión para las correlaciones más importantes
library(patchwork)

# r_at_Lmax vs scale: debería correlacionar fuerte (scale espacial del proceso)
p1 <- ggplot(Data_LGCP, aes(x = scale/1000, y = r_at_Lmax/1000)) +
  geom_point(alpha = 0.2) + geom_smooth(method = "lm", color = "red") +
  labs(x = "scale (km)", y = "r_at_Lmax (km)",
       title = paste0("r² = ", round(cor(Data_LGCP$scale, Data_LGCP$r_at_Lmax)^2, 3))) +
  theme_minimal()

# L_max vs var: amplitud del clustering relacionada con varianza del campo latente
p2 <- ggplot(Data_LGCP, aes(x = var, y = L_max/1000)) +
  geom_point(alpha = 0.2) + geom_smooth(method = "lm", color = "red") +
  labs(x = "var", y = "L_max (km)",
       title = paste0("r² = ", round(cor(Data_LGCP$var, Data_LGCP$L_max)^2, 3))) +
  theme_minimal()

# AUC vs var: clustering total acumulado relacionado con varianza
p3 <- ggplot(Data_LGCP, aes(x = var, y = AUC)) +
  geom_point(alpha = 0.2) + geom_smooth(method = "lm", color = "red") +
  labs(x = "var", y = "AUC",
       title = paste0("r² = ", round(cor(Data_LGCP$var, Data_LGCP$AUC)^2, 3))) +
  theme_minimal()

# slope_init vs var: pendiente inicial de L̂(r) captura el clustering local
p4 <- ggplot(Data_LGCP, aes(x = var, y = slope_init)) +
  geom_point(alpha = 0.2) + geom_smooth(method = "lm", color = "red") +
  labs(x = "var", y = "slope_init",
       title = paste0("r² = ", round(cor(Data_LGCP$var, Data_LGCP$slope_init,
                                         use = "complete.obs")^2, 3))) +
  theme_minimal()

print((p1 | p2) / (p3 | p4))

cat("\n¡Simulaciones con features completadas!\n")
cat("Directorio de salida:", out_dir, "\n")
