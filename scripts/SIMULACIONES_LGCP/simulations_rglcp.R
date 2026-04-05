################################################################################
# SIMULACIÓN LGCP — Ventana real Colombia (EPSG:3116)
# Versión FINAL: solo features de primer orden (intensidad local)
#
# Features (8):
#   - Quadrat-based: variance, dispersion index (VMR), max/min ratio
#   - Kernel density: variance, skewness, kurtosis, entropy, coef. variation
#
# Justificación: L(r) captura correlación de pares (segundo orden) pero
# NO la heterogeneidad espacial de intensidad. Para un LGCP con campo
# latente Z(s) ~ GP(mu, C), la distribución local de intensidad
# lambda(s) = exp(Z(s)) depende directamente de (mu, sigma², scale).
# Los momentos del campo suavizado capturan esta información desde un
# ángulo complementario a L(r).
################################################################################

library(pbmcapply)
library(tidyverse)
library(spatstat)
library(sf)

set.seed(123)
ncores <- 10

# =============================================================================
# CARGAR VENTANA
# =============================================================================

shapeZona_sp <- st_read("~/Documents/Personal/TesisUDFJCMCIC/solucion2026/data_new/clip_zona_continental_simplificado.geojson")
shapeZona_sp <- st_simplify(shapeZona_sp, dTolerance = 40000,
                            preserveTopology = TRUE)
r_iso_owin <- as.owin(shapeZona_sp)

area_win <- area.owin(r_iso_owin)
cat("Área de la ventana:", area_win, "m²\n")

# =============================================================================
# PARÁMETROS — calibrados a sismicidad Colombia 2020
# =============================================================================

ntrain <- 15000

set.seed(42)
mu    <- runif(ntrain, -19.5, -18.5)
var   <- runif(ntrain,   1.0,   3.0)
scale <- runif(ntrain, 80000, 200000)

N_esp <- area_win * exp(mu + var / 2)
cat("\n===== RANGOS DE PARÁMETROS =====\n")
cat("mu         :", round(range(mu),         3), "\n")
cat("var        :", round(range(var),         3), "\n")
cat("scale (km) :", round(range(scale)/1000,  1), "\n")
cat("E[N]       :", round(range(N_esp)),          "\n\n")

RMAX  <- 200000   # 200 km
NRVAL <- 128

# Bandwidth fijo para kernel density — 50 km, del orden del scale mínimo.
KDE_BANDWIDTH <- 50000
# Grilla para quadratcount — 5x5 sub-regiones dentro del bounding box
NQUAD <- 5

# =============================================================================
# FUNCIÓN: Extraer 8 features de primer orden (intensidad local)
# =============================================================================

extract_intensity_features <- function(pp, nquad = NQUAD, bw = KDE_BANDWIDTH) {

  # --- Quadrat-based features ---
  qc <- tryCatch({
    counts <- as.vector(quadratcount(pp, nx = nquad, ny = nquad))
    counts <- counts[counts >= 0]
    counts
  }, error = function(e) NULL)

  if (is.null(qc) || length(qc) < 4) {
    quad_var         <- 0
    quad_VMR         <- 1
    quad_range_ratio <- 1
  } else {
    quad_var  <- var(qc)
    quad_mean <- mean(qc)
    quad_VMR  <- if (quad_mean > 0) quad_var / quad_mean else 1
    quad_range_ratio <- max(qc) / (min(qc) + 1)
  }

  # --- Kernel density features ---
  kde <- tryCatch(
    density.ppp(pp, sigma = bw, dimyx = c(64, 64)),
    error = function(e) NULL
  )

  if (is.null(kde)) {
    return(list(
      quad_var = quad_var, quad_VMR = quad_VMR,
      quad_range_ratio = quad_range_ratio,
      kde_var = 0, kde_skew = 0, kde_kurt = 0,
      kde_entropy = 0, kde_cv = 0
    ))
  }

  vals <- as.vector(kde$v)
  vals <- vals[!is.na(vals)]
  vals <- pmax(vals, 0)

  kde_mean <- mean(vals)
  kde_sd   <- sd(vals)
  kde_var  <- kde_sd^2

  if (kde_sd > 0) {
    z_vals   <- (vals - kde_mean) / kde_sd
    kde_skew <- mean(z_vals^3)
    kde_kurt <- mean(z_vals^4) - 3
  } else {
    kde_skew <- 0
    kde_kurt <- 0
  }

  vals_sum <- sum(vals)
  if (vals_sum > 0) {
    p <- vals / vals_sum
    p <- p[p > 0]
    kde_entropy <- -sum(p * log(p)) / log(length(p))
  } else {
    kde_entropy <- 0
  }

  kde_cv <- if (kde_mean > 0) kde_sd / kde_mean else 0

  list(
    quad_var         = quad_var,
    quad_VMR         = quad_VMR,
    quad_range_ratio = quad_range_ratio,
    kde_var          = kde_var,
    kde_skew         = kde_skew,
    kde_kurt         = kde_kurt,
    kde_entropy      = kde_entropy,
    kde_cv           = kde_cv
  )
}

# =============================================================================
# FUNCIÓN DE SIMULACIÓN
# =============================================================================

run_one_sim <- function(mu, var, scale, win,
                        min_points = 30, max_points = 150000) {

  pp <- tryCatch(
    spatstat.random::rLGCP(
      model      = "matern",
      nu         = 1,
      mu         = mu,
      var        = var,
      scale      = scale,
      win        = win,
      dimyx      = c(128, 128),
      saveLambda = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(pp)) return(NULL)

  N <- spatstat.geom::npoints(pp)
  if (N < min_points || N > max_points) return(NULL)

  # --- L(r) ---
  L_obj <- tryCatch(
    spatstat.explore::Lest(
      pp,
      correction = "border",
      rmax       = RMAX,
      nrval      = NRVAL
    ),
    error = function(e) NULL
  )
  if (is.null(L_obj)) return(NULL)

  Lc <- L_obj$border - L_obj$r
  if (any(is.na(Lc))) return(NULL)

  # --- Features de primer orden ---
  feats_I <- extract_intensity_features(pp)

  list(
    mu = mu, var = var, scale = scale, N = N,
    r = L_obj$r, L = Lc,
    quad_var         = feats_I$quad_var,
    quad_VMR         = feats_I$quad_VMR,
    quad_range_ratio = feats_I$quad_range_ratio,
    kde_var          = feats_I$kde_var,
    kde_skew         = feats_I$kde_skew,
    kde_kurt         = feats_I$kde_kurt,
    kde_entropy      = feats_I$kde_entropy,
    kde_cv           = feats_I$kde_cv
  )
}

# Helper para construir tibble
sims_to_tibble <- function(sims) {
  purrr::map_dfr(sims, ~ tibble(
    mu    = .x$mu,    var   = .x$var,   scale = .x$scale,
    N     = .x$N,
    r     = list(.x$r),
    L     = list(.x$L),
    quad_var         = .x$quad_var,
    quad_VMR         = .x$quad_VMR,
    quad_range_ratio = .x$quad_range_ratio,
    kde_var          = .x$kde_var,
    kde_skew         = .x$kde_skew,
    kde_kurt         = .x$kde_kurt,
    kde_entropy      = .x$kde_entropy,
    kde_cv           = .x$kde_cv
  ))
}

# =============================================================================
# SIMULACIÓN POR CHUNKS — train
# =============================================================================

chunk_size  <- 50
n_chunks    <- ceiling(ntrain / chunk_size)
out_dir      <- "~/Documents/Personal/TesisUDFJCMCIC/PROPUESTA_EVENTOS_SISMICOS/Resultados/results_lgcp_features/TRAIN_FINAL"
out_dir_test <- "~/Documents/Personal/TesisUDFJCMCIC/PROPUESTA_EVENTOS_SISMICOS/Resultados/results_lgcp_features/TEST_FINAL"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_test, showWarnings = FALSE, recursive = TRUE)

total_valid <- 0

for (k in seq_len(n_chunks)) {
  cat("Train — chunk", k, "of", n_chunks, "\n")

  idx <- ((k - 1) * chunk_size + 1):min(k * chunk_size, ntrain)

  sims_k <- pbmcmapply(
    FUN      = function(mu_i, var_i, scale_i)
      run_one_sim(mu_i, var_i, scale_i, win = r_iso_owin),
    mu[idx], var[idx], scale[idx],
    SIMPLIFY = FALSE, mc.cores = ncores
  )

  sims_k <- sims_k[!sapply(sims_k, is.null)]
  if (length(sims_k) == 0) {
    cat("  Warning: sin simulaciones válidas en chunk", k, "\n"); next
  }

  Data_k      <- sims_to_tibble(sims_k)
  total_valid <- total_valid + nrow(Data_k)
  cat("  Válidas:", nrow(Data_k), "/", length(idx),
      "| Total:", total_valid, "\n")

  saveRDS(Data_k, file.path(out_dir, sprintf("Data_LGCP_train_%04d.rds", k)))
  gc()
}

cat("\n===== TRAIN COMPLETADO — total válidas:", total_valid, "=====\n\n")

# =============================================================================
# SIMULACIÓN POR CHUNKS — test
# =============================================================================
set.seed(1234)
ntest    <- 1000
mu_te    <- runif(ntest, -19.5, -18.5)
var_te   <- runif(ntest,   1.0,   3.0)
scale_te <- runif(ntest, 80000, 200000)

n_chunks_te    <- ceiling(ntest / chunk_size)
total_valid_te <- 0

for (k in seq_len(n_chunks_te)) {
  cat("Test — chunk", k, "of", n_chunks_te, "\n")

  idx <- ((k - 1) * chunk_size + 1):min(k * chunk_size, ntest)

  sims_k <- pbmcmapply(
    FUN      = function(mu_i, var_i, scale_i)
      run_one_sim(mu_i, var_i, scale_i, win = r_iso_owin),
    mu_te[idx], var_te[idx], scale_te[idx],
    SIMPLIFY = FALSE, mc.cores = ncores
  )

  sims_k <- sims_k[!sapply(sims_k, is.null)]
  if (length(sims_k) == 0) {
    cat("  Warning: sin simulaciones válidas en chunk", k, "\n"); next
  }

  Data_k <- sims_to_tibble(sims_k)
  total_valid_te <- total_valid_te + nrow(Data_k)
  cat("  Válidas:", nrow(Data_k), "/", length(idx),
      "| Total:", total_valid_te, "\n")

  saveRDS(Data_k, file.path(out_dir_test, sprintf("Data_LGCP_test_%04d.rds", k)))
  gc()
}

cat("\n===== TEST COMPLETADO — total válidas:", total_valid_te, "=====\n\n")

# =============================================================================
# CARGAR DATASETS COMPLETOS
# =============================================================================

files_tr  <- list.files(out_dir, pattern = "Data_LGCP_train_.*\\.rds$",
                        full.names = TRUE)
Data_LGCP <- purrr::map_dfr(files_tr, readRDS)
cat("Train — filas:", nrow(Data_LGCP), "\n")

files_te       <- list.files(out_dir_test, pattern = "Data_LGCP_test_.*\\.rds$",
                             full.names = TRUE)
Data_LGCP_test <- purrr::map_dfr(files_te, readRDS)
cat("Test  — filas:", nrow(Data_LGCP_test), "\n")
