#' ==============================================================================
#' CNN Training and Parameter Recovery
#' ==============================================================================
#'
#' Two models for estimating LGCP parameters (mu, sigma2, scale):
#'   M1: Baseline CNN — inputs: L(r), N
#'   M2: CNN + 8 intensity features — inputs: L(r), N, intensity features
#'
#' ==============================================================================

#' Needs a Keras 2 backend and a fresh R session (see README).

library(reticulate)
use_virtualenv("~/.virtualenvs/r-tensorflow", required = TRUE)  # must precede library(keras)

library(keras)
Sys.setenv(TF_USE_LEGACY_KERAS = "0")

library(tidyverse)
library(patchwork)
library(spatstat)
library(spatstat.geom)
library(sf)
local({
  kv <- as.character(reticulate::py_get_attr(reticulate::import("keras"), "__version__"))
  if (!startsWith(kv, "2.")) stop("Keras ", kv, " detected; this script needs Keras 2.x. See README.")
  message("Keras backend OK: ", kv)
})

set.seed(123)
tensorflow::set_random_seed(123)

setwd("/home/jasonromeroia/Documents/personal/paper_lgcp_features_computers-and-geosciences/lgcp-cnn-features")

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

out_dir      <- "Results_simulation/TRAIN"
out_dir_test <- "Results_simulation/TEST"

files_train <- list.files(out_dir, pattern = "Data_LGCP_train_.*\\.rds$", full.names = TRUE)
data_train_full <- purrr::map_dfr(files_train, readRDS)
message(sprintf("Training data loaded: %d rows", nrow(data_train_full)))

files_test <- list.files(out_dir_test, pattern = "Data_LGCP_test_.*\\.rds$", full.names = TRUE)
data_test_full <- purrr::map_dfr(files_test, readRDS)
message(sprintf("Test data loaded: %d rows", nrow(data_test_full)))

# Remove records with NA in L curve
data_train_full <- data_train_full[!sapply(data_train_full$L, function(x) any(is.na(x))), ]
data_test_full  <- data_test_full[!sapply(data_test_full$L, function(x) any(is.na(x))), ]

# ============================================================================
# 2. TRAIN / VALIDATION / TEST SPLIT
# ============================================================================


set.seed(123)
idx_train <- sample(seq_len(nrow(data_train_full)), 
                    size = floor(0.8 * nrow(data_train_full)))

data_train <- data_train_full[idx_train, ]
data_val   <- data_train_full[-idx_train, ]

message(sprintf("Internal train: %d | Validation: %d | Test: %d",
                nrow(data_train), nrow(data_val), nrow(data_test_full)))

# ============================================================================
# 3. NORMALIZATION (using train set statistics)
# ============================================================================

# Normalize L(r) curves
m_L  <- mean(unlist(data_train$L))
sd_L <- sd(unlist(data_train$L))
if (sd_L == 0) sd_L <- 1

make_L_array <- function(L_list, m_L, sd_L) {
  L_scaled <- lapply(L_list, function(L) (L - m_L) / sd_L)
  array_reshape(L_scaled, c(length(L_list), length(L_list[[1]]), 1))
}

train_L <- make_L_array(data_train$L, m_L, sd_L)
val_L   <- make_L_array(data_val$L, m_L, sd_L)
test_L  <- make_L_array(data_test_full$L, m_L, sd_L)

# Normalize N (point counts)
train_N_mat <- as.matrix(select(data_train, N))
val_N_mat   <- as.matrix(select(data_val, N))
test_N_mat  <- as.matrix(select(data_test_full, N))

m_N  <- apply(train_N_mat, 2, mean)
sd_N <- apply(train_N_mat, 2, sd)
if (sd_N == 0) sd_N <- 1

train_N <- scale(train_N_mat, center = m_N, scale = sd_N)
val_N   <- scale(val_N_mat, center = m_N, scale = sd_N)
test_N  <- scale(test_N_mat, center = m_N, scale = sd_N)

# Normalize 8 intensity features
feat_cols <- c("quad_var", "quad_VMR", "quad_range_ratio",
               "kde_var", "kde_skew", "kde_kurt", "kde_entropy", "kde_cv")

train_feat_mat <- as.matrix(select(data_train, all_of(feat_cols)))
val_feat_mat   <- as.matrix(select(data_val, all_of(feat_cols)))
test_feat_mat  <- as.matrix(select(data_test_full, all_of(feat_cols)))

m_feat  <- apply(train_feat_mat, 2, mean)
sd_feat <- apply(train_feat_mat, 2, sd)
sd_feat[sd_feat == 0] <- 1

train_feat <- scale(train_feat_mat, center = m_feat, scale = sd_feat)
val_feat   <- scale(val_feat_mat, center = m_feat, scale = sd_feat)
test_feat  <- scale(test_feat_mat, center = m_feat, scale = sd_feat)

# Normalize targets
train_par_raw <- as.matrix(select(data_train, mu:scale))
val_par_raw   <- as.matrix(select(data_val, mu:scale))
test_par_raw  <- as.matrix(select(data_test_full, mu:scale))

m_par  <- apply(train_par_raw, 2, mean)
sd_par <- apply(train_par_raw, 2, sd)
if (any(sd_par == 0)) sd_par[sd_par == 0] <- 1

train_par <- scale(train_par_raw, center = m_par, scale = sd_par)
val_par   <- scale(val_par_raw, center = m_par, scale = sd_par)
test_par  <- scale(test_par_raw, center = m_par, scale = sd_par)

true_par <- test_par_raw

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Convolutional block: 3 conv layers with batch norm and pooling
build_conv_branch <- function(input_layer, input_shape) {
  input_layer %>%
    layer_conv_1d(filters = 64, kernel_size = 7, activation = "relu",
                  input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    layer_max_pooling_1d(pool_size = 5) %>%
    layer_conv_1d(filters = 64, kernel_size = 7, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_max_pooling_1d(pool_size = 5) %>%
    layer_conv_1d(filters = 64, kernel_size = 7, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_flatten()
}

# Training callbacks: early stopping and learning rate reduction
make_callbacks <- function(es_patience = 15, lr_patience = 7) {
  list(
    callback_early_stopping(monitor = "val_loss", patience = es_patience,
                           restore_best_weights = TRUE, verbose = 1),
    callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5,
                                 patience = lr_patience, min_lr = 1e-6, verbose = 1)
  )
}

EPOCHS <- 200
BATCH_SIZE <- 64

# ============================================================================
# 4. MODEL 1: BASELINE CNN
# ============================================================================

message("\n--- Training Model 1: Baseline CNN ---")

main_input <- layer_input(shape = dim(train_L)[-1], name = "L_input")
aux_input  <- layer_input(shape = c(1), name = "N_input")

conv_branch <- build_conv_branch(main_input, dim(train_L)[-1])

output <- layer_concatenate(c(conv_branch, aux_input)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = ncol(train_par), activation = "linear")

model1 <- keras_model(inputs = c(main_input, aux_input), outputs = output)
model1 %>% compile(loss = "mse", optimizer = optimizer_adam(learning_rate = 1e-3),
                   metrics = list("mae"))

history1 <- model1 %>% fit(
  x = list(train_L, train_N), y = train_par,
  epochs = EPOCHS, batch_size = BATCH_SIZE,
  validation_data = list(list(val_L, val_N), val_par),
  callbacks = make_callbacks(), verbose = 1
)

pred1_std <- predict(model1, list(test_L, test_N), verbose = 0)
colnames(pred1_std) <- c("mu", "var", "scale")
pred1 <- sweep(sweep(pred1_std, 2, sd_par, `*`), 2, m_par, `+`)

# ============================================================================
# 5. MODEL 2: CNN + INTENSITY FEATURES (PROPOSED)
# ============================================================================

message("\n--- Training Model 2: CNN + Intensity Features ---")

main_input2 <- layer_input(shape = dim(train_L)[-1], name = "L_input2")
aux_input2  <- layer_input(shape = c(1), name = "N_input2")
feat_input  <- layer_input(shape = c(length(feat_cols)), name = "feat_input")

conv_branch2 <- build_conv_branch(main_input2, dim(train_L)[-1])

# Auxiliary branch for intensity features
feat_branch <- feat_input %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_dense(units = 16, activation = "relu")

output2 <- layer_concatenate(c(conv_branch2, aux_input2, feat_branch)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = ncol(train_par), activation = "linear")

model2 <- keras_model(inputs = c(main_input2, aux_input2, feat_input), outputs = output2)
model2 %>% compile(loss = "mse", optimizer = optimizer_adam(learning_rate = 1e-3),
                   metrics = list("mae"))

history2 <- model2 %>% fit(
  x = list(train_L, train_N, train_feat), y = train_par,
  epochs = EPOCHS, batch_size = BATCH_SIZE,
  validation_data = list(list(val_L, val_N, val_feat), val_par),
  callbacks = make_callbacks(), verbose = 1
)

pred2_std <- predict(model2, list(test_L, test_N, test_feat), verbose = 0)
colnames(pred2_std) <- c("mu", "var", "scale")
pred2 <- sweep(sweep(pred2_std, 2, sd_par, `*`), 2, m_par, `+`)

# ============================================================================
# 6. COMPUTE METRICS
# ============================================================================

# R2 = 1 - SSE/SST (not cor^2, not SSR/SST: they differ outside OLS)
r2_score <- function(true_vec, pred_vec) {
  ss_res <- sum((pred_vec - true_vec)^2)
  ss_tot <- sum((true_vec - mean(true_vec))^2)
  1 - ss_res / ss_tot
}

# Alias used by the scatter plots
r2_ssr_sst <- r2_score

compute_metrics <- function(true_mat, pred_mat, model_name) {
  purrr::map_dfr(c("mu", "var", "scale"), function(param) {
    errors <- pred_mat[, param] - true_mat[, param]
    tibble(
      model = model_name, param = param,
      R2 = round(r2_score(true_mat[, param], pred_mat[, param]), 4),
      RMSE = round(sqrt(mean(errors^2)), 4),
      MAE = round(mean(abs(errors)), 4)
    )
  })
}

# Join key for the LaTeX table, the fill scale and the legend: a mismatch
# silently turns the model into NA
MODEL_BASE <- "CNN de referencia"
MODEL_FEAT <- "CNN + descriptores"
models_order <- c(MODEL_BASE, MODEL_FEAT)

metrics1 <- compute_metrics(test_par, pred1_std, MODEL_BASE)
metrics2 <- compute_metrics(test_par, pred2_std, MODEL_FEAT)
all_metrics <- bind_rows(metrics1, metrics2)

stopifnot(all(all_metrics$model %in% models_order))

message("\nTest metrics (standardized scale):")
print(all_metrics, n = Inf)

# =============================================================================
# 7. LaTeX TABLE
# =============================================================================
param_latex <- c(mu = "$\\mu$", var = "$\\sigma^2$", scale = "$s$")

tex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Métricas de recuperación de parámetros en la escala estandarizada",
  "para el conjunto de prueba --- ventana de Colombia.}",
  "\\label{tab:metrics_final}",
  "\\begin{tabular}{llrrr}",
  "\\hline",
  "Modelo & Parámetro & $R^2$ & RMSE & MAE \\\\",
  "\\hline"
)

for (mod in models_order) {
  rows <- all_metrics %>% filter(model == mod)
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    pl  <- param_latex[row$param]
    tex_lines <- c(tex_lines, sprintf(
      "%s & %s & %.4f & %.4f & %.4f \\\\",
      row$model, pl, row$R2, row$RMSE, row$MAE
    ))
  }
  tex_lines <- c(tex_lines, "\\hline")
}

tex_lines <- c(tex_lines, "\\end{tabular}", "\\end{table}")

tex_out <- file.path(fig_dir, "metrics_table_final.tex")
writeLines(tex_lines, tex_out)
cat("\nLaTeX table saved to:", tex_out, "\n")

# =============================================================================
# 8. SCATTER PLOTS
# =============================================================================

make_scatter_grid <- function(true_mat, pred_mat, model_name) {
  param_labels <- c(mu    = expression(mu),
                    var   = expression(sigma^2),
                    scale = expression(italic(s)))
  plots <- list()
  
  for (p in c("mu", "var", "scale")) {
    df  <- tibble(true = true_mat[, p], pred = pred_mat[, p])
    r2  <- round(r2_ssr_sst(df$true, df$pred), 3)
    rng <- range(c(df$true, df$pred))
    
    plots[[p]] <- ggplot(df, aes(true, pred)) +
      geom_point(alpha = 0.2, size = 0.5, stroke = 0, colour = "steelblue4") +
      geom_abline(slope = 1, intercept = 0,
                  colour = "grey40", linetype = "dashed") +
      annotate("text",
               x     = rng[1], y = rng[2],
               label = paste0("R\u00b2 = ", r2),
               hjust = 0, vjust = 1, size = 3.5) +
      labs(x = "Valor verdadero", y = "Valor estimado", title = param_labels[[p]]) +
      theme_bw(base_size = 10) +
      theme(panel.grid  = element_blank(),
            plot.title  = element_text(hjust = 0.5, size = 11))
  }
  plots[["mu"]] | plots[["var"]] | plots[["scale"]]
}

p_s1 <- make_scatter_grid(true_par, pred1, MODEL_BASE) +
  plot_annotation(title = "CNN de referencia (Vihrs)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

p_s2 <- make_scatter_grid(true_par, pred2, MODEL_FEAT) +
  plot_annotation(title = "CNN + descriptores de intensidad (propuesta)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(fig_dir, "scatter_final_cnn_base.pdf"), p_s1,
       width = 9, height = 3.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "scatter_final_cnn_Ifeat.pdf"), p_s2,
       width = 9, height = 3.5, device = cairo_pdf)

p_combined <- (p_s1 / p_s2)
ggsave(file.path(fig_dir, "scatter_final_combined.pdf"), p_combined,
       width = 9, height = 7, device = cairo_pdf)

cat("Scatter plots saved\n")

# =============================================================================
# 9. LEARNING CURVES
# =============================================================================

loss_df <- bind_rows(
  tibble(epoch = seq_along(history1$metrics$loss),
         Entrenamiento = history1$metrics$loss,
         `Validación` = history1$metrics$val_loss,
         Modelo = MODEL_BASE),
  tibble(epoch = seq_along(history2$metrics$loss),
         Entrenamiento = history2$metrics$loss,
         `Validación` = history2$metrics$val_loss,
         Modelo = MODEL_FEAT)
) %>%
  pivot_longer(c(Entrenamiento, `Validación`),
               names_to = "Conjunto", values_to = "value")

p_loss <- ggplot(loss_df, aes(epoch, value, colour = Conjunto, linetype = Modelo)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Entrenamiento" = "steelblue",
                                "Validación" = "tomato")) +
  labs(x = "Época", y = "MSE", colour = "Conjunto", linetype = "Modelo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "loss_final_combined.pdf"), p_loss,
       width = 7, height = 4.5, device = cairo_pdf)

model_cnn <- subset(loss_df, Modelo == MODEL_BASE)
p_loss_cnn <- ggplot(model_cnn, aes(epoch, value, colour = Conjunto, linetype = Modelo)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Entrenamiento" = "steelblue",
                                "Validación" = "tomato")) +
  labs(x = "Época", y = "MSE", colour = "Conjunto", linetype = "Modelo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "p_loss_cnn.pdf"), p_loss_cnn,
       width = 7, height = 4.5, device = cairo_pdf)

model_cnn_improved <- subset(loss_df, Modelo == MODEL_FEAT)
p_loss_cnn_improved <- ggplot(model_cnn_improved, aes(epoch, value, colour = Conjunto, linetype = Modelo)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Entrenamiento" = "steelblue",
                                "Validación" = "tomato")) +
  labs(x = "Época", y = "MSE", colour = "Conjunto", linetype = "Modelo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "p_loss_cnn_improved.pdf"), p_loss_cnn_improved,
       width = 7, height = 4.5, device = cairo_pdf)

mae_df <- bind_rows(
  tibble(epoch = seq_along(history1$metrics$mae),
         Entrenamiento = history1$metrics$mae,
         `Validación` = history1$metrics$val_mae,
         Modelo = MODEL_BASE),
  tibble(epoch = seq_along(history2$metrics$mae),
         Entrenamiento = history2$metrics$mae,
         `Validación` = history2$metrics$val_mae,
         Modelo = MODEL_FEAT)
) %>%
  pivot_longer(c(Entrenamiento, `Validación`),
               names_to = "Conjunto", values_to = "value")

p_mae <- ggplot(mae_df, aes(epoch, value, colour = Conjunto, linetype = Modelo)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Entrenamiento" = "steelblue",
                                "Validación" = "tomato")) +
  labs(x = "Época", y = "MAE", colour = "Conjunto", linetype = "Modelo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "mae_final_combined.pdf"), p_mae,
       width = 7, height = 4.5, device = cairo_pdf)

cat("Loss and MAE plots saved\n")

# =============================================================================
# 10. R^2 COMPARISON BAR CHART
# =============================================================================

p_r2 <- ggplot(
  all_metrics %>%
    mutate(
      param = factor(param, levels = c("mu", "var", "scale"),
                     labels = c("\u03bc", "\u03c3\u00b2", "s")),
      model = factor(model, levels = models_order)
    ),
  aes(x = param, y = R2, fill = model)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = setNames(c("grey65", "darkorange"), models_order)) +
  labs(x = "Parámetro", y = expression(R^2), fill = "Modelo") +
  ylim(0, 1.08) +
  theme_bw(base_size = 11) +
  theme(panel.grid      = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "r2_comparison_final.pdf"), p_r2,
       width = 5, height = 4, device = cairo_pdf)

cat("R-squared comparison saved\n")
cat("\nDone! All figures and table in:", fig_dir, "\n")


cat("===== Loading seismic catalogue =====\n")

path_file_seismic <- "/home/jasonromeroia/Documents/personal/Tesis_MCIC/lgcp_inla_spde_sismos/Data/gdf_espacial_2020.gpkg" #"Data/EventosColPointsPlanas31162005_2020_continental.gpkg"

#sismos_sf <- readRDS("data/sismos_sf_2020.rds")
sismos_sf <- st_read(path_file_seismic)
#sismos_sf$year <- sismosSp$YEAR

# Filtrar solo el año 2020 (modelo espacial puro, sin estructura temporal)
#sismos_sf <- subset(sismos_sf, year >= 2020)
#sismos_sf <- subset(sismos_sf, year <= 2020)
cat("Total events:", nrow(sismos_sf), "\n")
cat("Columns:", paste(names(sismos_sf), collapse = ", "), "\n")
cat("Available years:", paste(sort(unique(sismos_sf$YEAR)), collapse = ", "), "\n")

# Subset year 2020
sismos_2020 <- subset(sismos_sf, YEAR == 2020)
cat("2020 events:", nrow(sismos_2020), "\n")

# Load window (same as in simulations)
shapeZona_sp <- readRDS("data/shapeZona_sp.rds")
win_colombia <- as.owin(shapeZona_sp)

# Build ppp
pSisAux  <- as.ppp(sismos_2020)
pSismos  <- ppp(pSisAux$x, pSisAux$y, window = win_colombia)
unitname(pSismos) <- c("meter", "meters")

N_obs <- npoints(pSismos)
cat("Points inside window:", N_obs, "\n")

# =============================================================================
# 2. COMPUTE L(r)-r AND 8 INTENSITY FEATURES
# =============================================================================

cat("\n===== Computing summary statistics of the observed pattern =====\n")

RMAX  <- 200000

# L(r) grid length: `nrval` is ignored by Lest(), spatstat derives 513 from rmax
NRVAL <- 513

KDE_BANDWIDTH <- 50000
NQUAD <- 5

# --- L(r) - r ---
L_obj <- Lest(pSismos, correction = "border", rmax = RMAX)
stopifnot(length(L_obj$r) == NRVAL)
Lc_obs <- L_obj$border - L_obj$r
r_obs  <- L_obj$r

cat("L(r) computed: ", length(Lc_obs), "values\n")
cat("D(r) range: [", round(min(Lc_obs), 1), ",", round(max(Lc_obs), 1), "]\n")

# --- 8 intensity features ---
extract_intensity_features <- function(pp, nquad = NQUAD, bw = KDE_BANDWIDTH) {
  qc <- tryCatch({
    counts <- as.vector(quadratcount(pp, nx = nquad, ny = nquad))
    counts <- counts[counts >= 0]
    counts
  }, error = function(e) NULL)
  
  if (is.null(qc) || length(qc) < 4) {
    quad_var <- 0; quad_VMR <- 1; quad_range_ratio <- 1
  } else {
    quad_var  <- var(qc)
    quad_mean <- mean(qc)
    quad_VMR  <- if (quad_mean > 0) quad_var / quad_mean else 1
    quad_range_ratio <- max(qc) / (min(qc) + 1)
  }
  
  kde <- tryCatch(
    density.ppp(pp, sigma = bw, dimyx = c(64, 64)),
    error = function(e) NULL
  )
  
  if (is.null(kde)) {
    return(list(quad_var = quad_var, quad_VMR = quad_VMR,
                quad_range_ratio = quad_range_ratio,
                kde_var = 0, kde_skew = 0, kde_kurt = 0,
                kde_entropy = 0, kde_cv = 0))
  }
  
  vals <- as.vector(kde$v)
  vals <- vals[!is.na(vals)]
  vals <- pmax(vals, 0)
  
  kde_mean <- mean(vals); kde_sd <- sd(vals); kde_var <- kde_sd^2
  
  if (kde_sd > 0) {
    z_vals <- (vals - kde_mean) / kde_sd
    kde_skew <- mean(z_vals^3)
    kde_kurt <- mean(z_vals^4) - 3
  } else { kde_skew <- 0; kde_kurt <- 0 }
  
  vals_sum <- sum(vals)
  if (vals_sum > 0) {
    p <- vals / vals_sum; p <- p[p > 0]
    kde_entropy <- -sum(p * log(p)) / log(length(p))
  } else { kde_entropy <- 0 }
  
  kde_cv <- if (kde_mean > 0) kde_sd / kde_mean else 0
  
  list(quad_var = quad_var, quad_VMR = quad_VMR,
       quad_range_ratio = quad_range_ratio,
       kde_var = kde_var, kde_skew = kde_skew, kde_kurt = kde_kurt,
       kde_entropy = kde_entropy, kde_cv = kde_cv)
}

feats_obs <- extract_intensity_features(pSismos)
cat("\nObserved features:\n")
print(as.data.frame(feats_obs))

# =============================================================================
# 3. NORMALIZE AND PREDICT
# =============================================================================
# Needs the normalization stats and both models from the sections above

cat("\n===== Prediction with trained models =====\n")

# --- Normalize L(r) ---
obs_L <- array(
  (Lc_obs - m_L) / sd_L,
  dim = c(1, length(Lc_obs), 1)
)

# --- Normalize N ---
obs_N <- matrix((N_obs - m_N) / sd_N, nrow = 1)

# --- Normalize features ---
# Index by feat_cols: R subtracts named vectors by position, not by name
feat_vec <- unlist(feats_obs)[feat_cols]
stopifnot(length(feat_vec) == length(feat_cols), !anyNA(feat_vec))
obs_featI <- matrix(
  (feat_vec - m_feat) / sd_feat,
  nrow = 1
)

# --- Baseline CNN prediction ---
pred1_obs_std <- predict(model1, list(obs_L, obs_N))
pred1_obs <- pred1_obs_std * sd_par + m_par
colnames(pred1_obs) <- c("mu", "var", "scale")

# --- CNN + I-feat prediction ---
pred2_obs_std <- predict(model2, list(obs_L, obs_N, obs_featI))
pred2_obs <- pred2_obs_std * sd_par + m_par
colnames(pred2_obs) <- c("mu", "var", "scale")

cat("\n--- Estimated parameters ---\n")
cat("Baseline CNN: mu =", round(pred1_obs[1, "mu"], 4),
    " sigma2 =", round(pred1_obs[1, "var"], 4),
    " s =", round(pred1_obs[1, "scale"], 0), "m\n")
cat("CNN + I-feat: mu =", round(pred2_obs[1, "mu"], 4),
    " sigma2 =", round(pred2_obs[1, "var"], 4),
    " s =", round(pred2_obs[1, "scale"], 0), "m\n")

# Expected E[N]
EN_base  <- area.owin(win_colombia) * exp(pred1_obs[1, "mu"] + pred1_obs[1, "var"] / 2)
EN_feat  <- area.owin(win_colombia) * exp(pred2_obs[1, "mu"] + pred2_obs[1, "var"] / 2)
cat("\nE[N] Baseline CNN: ", round(EN_base), " (observed:", N_obs, ")\n")
cat("E[N] CNN + I-feat: ", round(EN_feat), " (observed:", N_obs, ")\n")

# Table of estimated parameters
params_df <- tibble(
  Modelo = c(MODEL_BASE, MODEL_FEAT),
  mu     = c(pred1_obs[1, "mu"],    pred2_obs[1, "mu"]),
  sigma2 = c(pred1_obs[1, "var"],   pred2_obs[1, "var"]),
  s_m    = c(pred1_obs[1, "scale"], pred2_obs[1, "scale"]),
  EN     = c(EN_base, EN_feat)
)
cat("\n"); print(params_df)

# =============================================================================
# 4. VALIDATION: SIMULATE LGCP WITH ESTIMATED PARAMETERS (proposed model)
# =============================================================================

cat("\n===== Validation: simulations with estimated parameters =====\n")

mu_hat     <- pred2_obs[1, "mu"]
sigma2_hat <- pred2_obs[1, "var"]
s_hat      <- pred2_obs[1, "scale"]

nsim_val <- 199  # for 95% envelope

cat("Simulating", nsim_val, "LGCP realisations with:\n")
cat("  mu =", round(mu_hat, 4), "\n")
cat("  sigma2 =", round(sigma2_hat, 4), "\n")
cat("  s =", round(s_hat, 0), "m\n\n")

# Semilla local: sin esto el estado del RNG depende de TODO lo anterior
# (entrenamiento incluido), asi que la envolvente cambiaba en cada corrida.
set.seed(2020)

sims_val <- lapply(seq_len(nsim_val), function(i) {
  if (i %% 50 == 0) cat("  Simulation", i, "of", nsim_val, "\n")
  tryCatch(
    rLGCP(
      model      = "matern",
      nu         = 1,
      mu         = mu_hat,
      var        = sigma2_hat,
      scale      = s_hat,
      win        = win_colombia,
      dimyx      = c(128, 128),
      saveLambda = FALSE
    ),
    error = function(e) NULL
  )
})

sims_val <- sims_val[!sapply(sims_val, is.null)]
cat("Valid simulations:", length(sims_val), "of", nsim_val, "\n")

# Simulation N counts
N_sims <- sapply(sims_val, npoints)
cat("Simulated N: median =", median(N_sims), ", range = [",
    min(N_sims), ",", max(N_sims), "]\n")
cat("Observed N:", N_obs, "\n")

# =============================================================================
# 5. GLOBAL ENVELOPE UNDER THE FITTED MODEL
# =============================================================================

cat("\n===== Computing L(r) envelope under the fitted model =====\n")

# L(r) for each simulation
L_sims <- lapply(sims_val, function(pp) {
  tryCatch({
    L_obj_s <- Lest(pp, correction = "border", rmax = RMAX)
    L_obj_s$border - L_obj_s$r
  }, error = function(e) NULL)
})
L_sims <- L_sims[!sapply(L_sims, is.null)]

# Pointwise envelope band
L_mat <- do.call(rbind, L_sims)
L_lo  <- apply(L_mat, 2, quantile, probs = 0.025)
L_hi  <- apply(L_mat, 2, quantile, probs = 0.975)
L_med <- apply(L_mat, 2, median)

# =============================================================================
# 6. PAPER FIGURES
# =============================================================================
cat("\n===== Generating figures =====\n")
# --- 6a. Observed L(r) envelope vs fitted model ---
df_env <- tibble(
  r     = r_obs / 1000,  # km
  D_obs = Lc_obs,
  D_med = L_med,
  D_lo  = L_lo,
  D_hi  = L_hi
)

p_envelope <- ggplot(df_env, aes(x = r)) +
  # 95% envelope ribbon
  geom_ribbon(aes(ymin = D_lo, ymax = D_hi, fill = "Envolvente 95%"),
              alpha = 0.3) +
  # Fitted model median (dashed)
  geom_line(aes(y = D_med, colour = "Modelo ajustado (mediana)",
                linetype = "Modelo ajustado (mediana)"),
            linewidth = 0.6) +
  # Observed (solid black)
  geom_line(aes(y = D_obs, colour = "Observado",
                linetype = "Observado"),
            linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted") +
  scale_colour_manual(
    name   = NULL,
    values = c("Observado" = "black",
               "Modelo ajustado (mediana)" = "steelblue")
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Observado" = "solid",
               "Modelo ajustado (mediana)" = "dashed")
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Envolvente 95%" = "steelblue")
  ) +
  labs(
    x = "Distancia r (km)",
    y = expression(hat(L)(r) - r),
    title = expression(paste(
      "Envolvente del 95% bajo el modelo ajustado (",
      hat(mu), ", ", hat(sigma)^2, ", ", hat(s), ")"
    ))
  ) +
  annotate("text",
           x = max(df_env$r) * 0.6,
           y = max(df_env$D_obs) * 0.9,
           label = paste0("N observado = ", N_obs,
                          "\nN simulado (mediana) = ", round(mean(N_sims))),
           hjust = 0, size = 3.5) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 10),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.key.width = unit(1.2, "cm")
  ) +
  guides(
    colour   = guide_legend(order = 1, override.aes = list(linewidth = 0.8)),
    linetype = guide_legend(order = 1),
    fill     = guide_legend(order = 2, override.aes = list(alpha = 0.3))
  )

ggsave(file.path(fig_dir, "envelope_fitted_model.pdf"), p_envelope,
       width = 7, height = 4.5, device = cairo_pdf)

# --- 6b. Histogram of simulated N vs observed ---
p_N <- ggplot(tibble(N = N_sims), aes(x = N)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.6, colour = "white") +
  geom_vline(xintercept = N_obs, colour = "red", linewidth = 1, linetype = "dashed") +
  annotate("text", x = N_obs, y = Inf, label = paste("N observado =", N_obs),
           vjust = 2, hjust = -0.1, colour = "red", size = 3.5) +
  labs(x = "N (número de puntos)", y = "Frecuencia",
       title = "Distribución de N bajo el modelo ajustado") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10))

ggsave(file.path(fig_dir, "hist_N_fitted_model.pdf"), p_N,
       width = 5, height = 4, device = cairo_pdf)

# --- 6c. LaTeX table of estimated parameters ---
tex_params <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Parámetros LGCP estimados para el catálogo sísmico de Colombia 2020",
  "($N=14{,}346$ eventos).}",
  "\\label{tab:params_obs}",
  "\\begin{tabular}{lrrrr}",
  "\\hline",
  "Modelo & $\\hat{\\mu}$ & $\\hat{\\sigma}^2$ & $\\hat{s}$ (m) & $\\mathrm{E}[N]$ \\\\",
  "\\hline",
  sprintf("CNN de referencia & %.4f & %.4f & %s & %s \\\\",
          pred1_obs[1, "mu"], pred1_obs[1, "var"],
          format(round(pred1_obs[1, "scale"]), big.mark = "{,}"),
          format(round(EN_base), big.mark = "{,}")),
  sprintf("CNN + descriptores & %.4f & %.4f & %s & %s \\\\",
          pred2_obs[1, "mu"], pred2_obs[1, "var"],
          format(round(pred2_obs[1, "scale"]), big.mark = "{,}"),
          format(round(EN_feat), big.mark = "{,}")),
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(tex_params, file.path(fig_dir, "params_table_obs.tex"))
cat("LaTeX table saved\n")

