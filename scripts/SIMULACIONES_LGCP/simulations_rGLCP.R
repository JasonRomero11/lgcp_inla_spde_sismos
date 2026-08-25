#' LGCP Simulation with Intensity Features
#' 
#' Simulate LGCP on Colombia boundary (EPSG:3116)
#' Extract 8 features for ML: quadrat-based + kernel density metrics

# Setup
library(pbmcapply)
library(tidyverse)
library(spatstat)
library(sf)

setwd("/home/jasonromeroia/Documents/personal/paper_lgcp_features_computers-and-geosciences/lgcp-cnn-features")

#source("/home/jasonromeroia/Documents/personal/Tesis_MCIC/lgcp_inla_spde_sismos/scripts/ENTRENANDO_CNN/replot_figures.R")

set.seed(123)
N_CORES <- 11
# Load spatial window
shapeZona_sp <- readRDS("data/shapeZona_sp.rds")
r_iso_owin <- as.owin(shapeZona_sp)
area_win <- area.owin(r_iso_owin)
message(sprintf("Window area: %.2e m²", area_win))


N_TRAIN <- 15000

# Ranges
N_min     <- 500       # smallest expected count
N_max     <- 60000     # largest  expected count
var_min   <- 0.5
var_max   <- 6.0
scale_min <- 20000     # 20 km
scale_max <- 300000    # 300 km

set.seed(42)
log_N_expected <- runif(N_TRAIN, log(N_min), log(N_max))
N_expected     <- exp(log_N_expected)
var            <- runif(N_TRAIN, var_min,   var_max)
scale          <- runif(N_TRAIN, scale_min, scale_max)
mu             <- log(N_expected / area_win) - var / 2

# -------- Sanity checks (training) ------------------------------------------
message("\n===== TRAINING PARAMETER RANGES =====")
message(sprintf("  mu          : [%7.3f , %7.3f]", min(mu),    max(mu)))
message(sprintf("  var         : [%7.3f , %7.3f]", min(var),   max(var)))
message(sprintf("  scale (km)  : [%7.1f , %7.1f]", min(scale)/1000, max(scale)/1000))
message(sprintf("  E[N]        : [%7.0f , %7.0f]", min(N_expected), max(N_expected)))



# Global parameters
RMAX <- 200000

# L(r) grid length: `nrval` is ignored by Lest(), spatstat derives 513 from rmax
NRVAL <- 513

KDE_BANDWIDTH <- 50000
NQUAD <- 5
MIN_POINTS <- 30
MAX_POINTS <- 150000
CHUNK_SIZE <- 50


# Extract intensity features from point pattern
extract_intensity_features <- function(pp, nquad = NQUAD, bw = KDE_BANDWIDTH) {
  # Quadrat counts
  qc <- tryCatch({
    counts <- as.vector(quadratcount(pp, nx = nquad, ny = nquad))
    counts <- counts[counts >= 0]
    counts
  }, error = function(e) NULL)
  
  if (is.null(qc) || length(qc) < 4) {
    quad_var <- 0; quad_VMR <- 1; quad_range_ratio <- 1
  } else {
    quad_var <- var(qc)
    quad_mean <- mean(qc)
    quad_VMR <- if (quad_mean > 0) quad_var / quad_mean else 1
    quad_range_ratio <- max(qc) / (min(qc) + 1)
  }
  
  # Kernel density
  kde <- tryCatch(
    density.ppp(pp, sigma = bw, dimyx = c(64, 64)),
    error = function(e) NULL
  )
  
  if (is.null(kde)) {
    return(list(quad_var = quad_var, quad_VMR = quad_VMR, quad_range_ratio = quad_range_ratio,
                kde_var = 0, kde_skew = 0, kde_kurt = 0, kde_entropy = 0, kde_cv = 0))
  }
  
  vals <- as.vector(kde$v)
  vals <- vals[!is.na(vals)]
  vals <- pmax(vals, 0)
  
  kde_mean <- mean(vals); kde_sd <- sd(vals); kde_var <- kde_sd^2
  
  if (kde_sd > 0) {
    z_vals <- (vals - kde_mean) / kde_sd
    kde_skew <- mean(z_vals^3)
    kde_kurt <- mean(z_vals^4) - 3
  } else {
    kde_skew <- 0; kde_kurt <- 0
  }
  
  vals_sum <- sum(vals)
  if (vals_sum > 0) {
    p <- vals / vals_sum; p <- p[p > 0]
    kde_entropy <- -sum(p * log(p)) / log(length(p))
  } else {
    kde_entropy <- 0
  }
  
  kde_cv <- if (kde_mean > 0) kde_sd / kde_mean else 0
  
  list(quad_var = quad_var, quad_VMR = quad_VMR, quad_range_ratio = quad_range_ratio,
       kde_var = kde_var, kde_skew = kde_skew, kde_kurt = kde_kurt,
       kde_entropy = kde_entropy, kde_cv = kde_cv)
}

# Run single LGCP simulation
run_one_sim <- function(mu, var, scale, win, min_points = MIN_POINTS, max_points = MAX_POINTS) {
  pp <- tryCatch({
    spatstat.random::rLGCP(model = "matern", nu = 1, mu = mu, var = var,
                          scale = scale, win = win, dimyx = c(128, 128), saveLambda = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(pp)) return(NULL)
  
  N <- spatstat.geom::npoints(pp)
  if (N < min_points || N > max_points) return(NULL)
  
  L_obj <- tryCatch(
    spatstat.explore::Lest(pp, correction = "border", rmax = RMAX),
    error = function(e) NULL
  )
  if (is.null(L_obj)) return(NULL)

  stopifnot(length(L_obj$r) == NRVAL)  # all curves must share the same grid

  Lc <- L_obj$border - L_obj$r
  if (any(is.na(Lc))) return(NULL)
  
  feats_I <- extract_intensity_features(pp)
  
  list(mu = mu, var = var, scale = scale, N = N,
       r = L_obj$r, L = Lc,
       quad_var = feats_I$quad_var, quad_VMR = feats_I$quad_VMR,
       quad_range_ratio = feats_I$quad_range_ratio,
       kde_var = feats_I$kde_var, kde_skew = feats_I$kde_skew,
       kde_kurt = feats_I$kde_kurt, kde_entropy = feats_I$kde_entropy,
       kde_cv = feats_I$kde_cv)
}

# Convert simulations to tibble
sims_to_tibble <- function(sims) {
  purrr::map_dfr(sims, ~ tibble(
    mu = .x$mu, var = .x$var, scale = .x$scale, N = .x$N,
    r = list(.x$r), L = list(.x$L),
    quad_var = .x$quad_var, quad_VMR = .x$quad_VMR, quad_range_ratio = .x$quad_range_ratio,
    kde_var = .x$kde_var, kde_skew = .x$kde_skew, kde_kurt = .x$kde_kurt,
    kde_entropy = .x$kde_entropy, kde_cv = .x$kde_cv
  ))
}

# Training simulations
OUT_DIR_TRAIN <- "Results_simulation/TRAIN"
OUT_DIR_TEST  <- "Results_simulation/TEST"

dir.create(OUT_DIR_TRAIN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_TEST, showWarnings = FALSE, recursive = TRUE)

n_chunks_train <- ceiling(N_TRAIN / CHUNK_SIZE)
total_valid_train <- 0

message(sprintf("\nTraining: %d sims in %d chunks\n", N_TRAIN, n_chunks_train))

for (k in seq_len(n_chunks_train)) {
  message(sprintf("Chunk %d / %d", k, n_chunks_train))
  
  idx <- ((k - 1) * CHUNK_SIZE + 1):min(k * CHUNK_SIZE, N_TRAIN)
  
  sims_k <- pbmcmapply(
    FUN = function(mu_i, var_i, scale_i) run_one_sim(mu_i, var_i, scale_i, win = r_iso_owin),
    mu[idx], var[idx], scale[idx],
    SIMPLIFY = FALSE, mc.cores = N_CORES
  )
  
  sims_k <- sims_k[!sapply(sims_k, is.null)]
  if (length(sims_k) == 0) {
    message(sprintf("  No valid sims in chunk %d", k))
    next
  }
  
  Data_k <- sims_to_tibble(sims_k)
  total_valid_train <- total_valid_train + nrow(Data_k)
  message(sprintf("  Valid: %d / %d (total: %d)", nrow(Data_k), length(idx), total_valid_train))
  
  saveRDS(Data_k, file.path(OUT_DIR_TRAIN, sprintf("Data_LGCP_train_%04d.rds", k)))
  gc()
}

message(sprintf("\nTraining completed: %d simulations saved\n", total_valid_train))


# =============================================================================
# Test set — SAME sampling scheme as training (independent seed)
# =============================================================================
N_TEST <- 1500

set.seed(1234)
log_N_expected_test <- runif(N_TEST, log(N_min), log(N_max))
N_expected_test     <- exp(log_N_expected_test)
var_test            <- runif(N_TEST, var_min,   var_max)
scale_test          <- runif(N_TEST, scale_min, scale_max)
mu_test             <- log(N_expected_test / area_win) - var_test / 2

message("\n===== TEST PARAMETER RANGES =====")
message(sprintf("  mu          : [%7.3f , %7.3f]", min(mu_test),    max(mu_test)))
message(sprintf("  var         : [%7.3f , %7.3f]", min(var_test),   max(var_test)))
message(sprintf("  scale (km)  : [%7.1f , %7.1f]", min(scale_test)/1000, max(scale_test)/1000))
message(sprintf("  E[N]        : [%7.0f , %7.0f]", min(N_expected_test), max(N_expected_test)))

n_chunks_test <- ceiling(N_TEST / CHUNK_SIZE)
total_valid_test <- 0

message(sprintf("Test: %d sims in %d chunks\n", N_TEST, n_chunks_test))

for (k in seq_len(n_chunks_test)) {
  message(sprintf("Chunk %d / %d", k, n_chunks_test))
  
  idx <- ((k - 1) * CHUNK_SIZE + 1):min(k * CHUNK_SIZE, N_TEST)
  
  sims_k <- pbmcmapply(
    FUN = function(mu_i, var_i, scale_i) run_one_sim(mu_i, var_i, scale_i, win = r_iso_owin),
    mu_test[idx], var_test[idx], scale_test[idx],
    SIMPLIFY = FALSE, mc.cores = N_CORES
  )
  
  sims_k <- sims_k[!sapply(sims_k, is.null)]
  if (length(sims_k) == 0) {
    message(sprintf("  No valid sims in chunk %d", k))
    next
  }
  
  Data_k <- sims_to_tibble(sims_k)
  total_valid_test <- total_valid_test + nrow(Data_k)
  message(sprintf("  Valid: %d / %d (total: %d)", nrow(Data_k), length(idx), total_valid_test))
  
  saveRDS(Data_k, file.path(OUT_DIR_TEST, sprintf("Data_LGCP_test_%04d.rds", k)))
  gc()
}

message(sprintf("\nTest completed: %d simulations saved\n", total_valid_test))

# Load all results
message("Loading results...")

files_train <- list.files(OUT_DIR_TRAIN, pattern = "Data_LGCP_train_.*\\.rds$", full.names = TRUE)
data_train <- purrr::map_dfr(files_train, readRDS)

files_test <- list.files(OUT_DIR_TEST, pattern = "Data_LGCP_test_.*\\.rds$", full.names = TRUE)
data_test <- purrr::map_dfr(files_test, readRDS)

message(sprintf("Train: %d rows | Test: %d rows", nrow(data_train), nrow(data_test)))


