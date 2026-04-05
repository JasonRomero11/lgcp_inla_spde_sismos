###############################################################################
# CNN_resultados_final.R
# Versión FINAL — solo 2 modelos:
#   M1: CNN base (réplica Vihrs) — entrada: L(r), N
#   M2: CNN + 8 features intensidad — entrada: L(r), N, f_I ∈ R^8
#
# Lee datos de TRAIN_FINAL / TEST_FINAL generados por
# simulaciones_final_v2.R
###############################################################################

library(keras)
library(tidyverse)
library(patchwork)

set.seed(123)
tensorflow::set_random_seed(123)

setwd("~/Documents/Personal/TesisUDFJCMCIC/PROPUESTA_EVENTOS_SISMICOS/paper_lgcp_features/")
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

out_dir      <- "~/Documents/Personal/TesisUDFJCMCIC/PROPUESTA_EVENTOS_SISMICOS/Resultados/results_lgcp_features/TRAIN_FINAL"
out_dir_test <- "~/Documents/Personal/TesisUDFJCMCIC/PROPUESTA_EVENTOS_SISMICOS/Resultados/results_lgcp_features/TEST_FINAL"

files_tr <- list.files(out_dir, pattern = "Data_LGCP_train_.*\\.rds$", full.names = TRUE)
Data_LGCP <- purrr::map_dfr(files_tr, readRDS)
cat("Train total rows:", nrow(Data_LGCP), "\n")

files_te <- list.files(out_dir_test, pattern = "Data_LGCP_test_.*\\.rds$", full.names = TRUE)
Data_LGCP_test <- purrr::map_dfr(files_te, readRDS)
cat("Test total rows:", nrow(Data_LGCP_test), "\n")

# Safety filters
Data_LGCP      <- Data_LGCP[!sapply(Data_LGCP$L, function(x) any(is.na(x))), ]
Data_LGCP_test <- Data_LGCP_test[!sapply(Data_LGCP_test$L, function(x) any(is.na(x))), ]

# =============================================================================
# 2. TRAIN / VALIDATION / TEST SPLIT
# =============================================================================

set.seed(123)
idx_train <- sample(seq_len(nrow(Data_LGCP)), size = floor(0.8 * nrow(Data_LGCP)))

Data_train <- Data_LGCP[idx_train, ]
Data_val   <- Data_LGCP[-idx_train, ]

cat("Internal train rows:", nrow(Data_train), "\n")
cat("Validation rows:    ", nrow(Data_val),   "\n")
cat("External test rows: ", nrow(Data_LGCP_test), "\n")

# =============================================================================
# 3. NORMALIZATION — TRAIN STATISTICS ONLY
# =============================================================================

# --- L(r) curve ---
m_L  <- mean(unlist(Data_train$L))
sd_L <- sd(unlist(Data_train$L))
if (sd_L == 0) sd_L <- 1

make_L_array <- function(L_list, m_L, sd_L) {
  L_scaled <- lapply(L_list, function(L) (L - m_L) / sd_L)
  array_reshape(L_scaled, c(length(L_list), length(L_list[[1]]), 1))
}

train_L <- make_L_array(Data_train$L, m_L, sd_L)
val_L   <- make_L_array(Data_val$L,   m_L, sd_L)
test_L  <- make_L_array(Data_LGCP_test$L, m_L, sd_L)

# --- N scalar ---
train_N_mat <- as.matrix(select(Data_train, N))
val_N_mat   <- as.matrix(select(Data_val, N))
test_N_mat  <- as.matrix(select(Data_LGCP_test, N))

m_N  <- apply(train_N_mat, 2, mean)
sd_N <- apply(train_N_mat, 2, sd)
sd_N[sd_N == 0] <- 1

train_N <- scale(train_N_mat, center = m_N, scale = sd_N)
val_N   <- scale(val_N_mat,   center = m_N, scale = sd_N)
test_N  <- scale(test_N_mat,  center = m_N, scale = sd_N)

# --- Intensity features (8) ---
feat_I_cols <- c("quad_var", "quad_VMR", "quad_range_ratio",
                 "kde_var", "kde_skew", "kde_kurt", "kde_entropy", "kde_cv")

train_featI_mat <- as.matrix(select(Data_train,     all_of(feat_I_cols)))
val_featI_mat   <- as.matrix(select(Data_val,       all_of(feat_I_cols)))
test_featI_mat  <- as.matrix(select(Data_LGCP_test, all_of(feat_I_cols)))

m_featI  <- apply(train_featI_mat, 2, mean)
sd_featI <- apply(train_featI_mat, 2, sd)
sd_featI[sd_featI == 0] <- 1

train_featI <- scale(train_featI_mat, center = m_featI, scale = sd_featI)
val_featI   <- scale(val_featI_mat,   center = m_featI, scale = sd_featI)
test_featI  <- scale(test_featI_mat,  center = m_featI, scale = sd_featI)

# --- Targets ---
train_par_raw <- as.matrix(select(Data_train,     mu:scale))
val_par_raw   <- as.matrix(select(Data_val,       mu:scale))
test_par_raw  <- as.matrix(select(Data_LGCP_test, mu:scale))

m_par  <- apply(train_par_raw, 2, mean)
sd_par <- apply(train_par_raw, 2, sd)
sd_par[sd_par == 0] <- 1

train_par <- scale(train_par_raw, center = m_par, scale = sd_par)
val_par   <- scale(val_par_raw,   center = m_par, scale = sd_par)
test_par  <- scale(test_par_raw,  center = m_par, scale = sd_par)

true_par <- test_par_raw

# =============================================================================
# HELPER: CNN convolutional branch with batch normalization
# =============================================================================

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

# =============================================================================
# HELPER: training callbacks
# =============================================================================

make_callbacks <- function(patience_es = 15, patience_lr = 7) {
  list(
    callback_early_stopping(
      monitor              = "val_loss",
      patience             = patience_es,
      restore_best_weights = TRUE,
      mode                 = "min",
      verbose              = 1
    ),
    callback_reduce_lr_on_plateau(
      monitor  = "val_loss",
      factor   = 0.5,
      patience = patience_lr,
      min_lr   = 1e-6,
      verbose  = 1
    )
  )
}

EPOCHS     <- 200
BATCH_SIZE <- 64

# =============================================================================
# 4. MODEL 1 — BASELINE CNN (Vihrs replica + batch norm)
# =============================================================================

cat("\n===== Model 1: Baseline CNN =====\n")

main_in_1 <- layer_input(shape = dim(train_L)[-1], name = "main_input_1")
aux_in_1  <- layer_input(shape = c(1),             name = "aux_input_1")

conv_1 <- build_conv_branch(main_in_1, dim(train_L)[-1])

out_1 <- layer_concatenate(c(conv_1, aux_in_1)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = ncol(train_par), activation = "linear")

model1 <- keras_model(inputs = c(main_in_1, aux_in_1), outputs = out_1)
model1 %>% compile(loss = "mse", optimizer = optimizer_adam(learning_rate = 1e-3),
                   metrics = list("mae"))
summary(model1)

history1 <- model1 %>% fit(
  x               = list(train_L, train_N),
  y               = train_par,
  epochs          = EPOCHS,
  batch_size      = BATCH_SIZE,
  validation_data = list(list(val_L, val_N), val_par),
  callbacks       = make_callbacks(),
  verbose         = 2
)

pred1_std <- predict(model1, list(test_L, test_N))
colnames(pred1_std) <- c("mu", "var", "scale")
pred1 <- sweep(sweep(pred1_std, 2, sd_par, `*`), 2, m_par, `+`)
colnames(pred1) <- c("mu", "var", "scale")

# =============================================================================
# 5. MODEL 2 — CNN + 8 features intensidad (PROPUESTA)
# =============================================================================

cat("\n===== Model 2: CNN + 8 intensity features =====\n")

main_in_2 <- layer_input(shape = dim(train_L)[-1],      name = "main_input_2")
aux_in_2  <- layer_input(shape = c(1),                  name = "aux_input_2")
feat_in_2 <- layer_input(shape = c(length(feat_I_cols)), name = "feat_input_2")

conv_2 <- build_conv_branch(main_in_2, dim(train_L)[-1])

feat_branch_2 <- feat_in_2 %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_dense(units = 16, activation = "relu")

out_2 <- layer_concatenate(c(conv_2, aux_in_2, feat_branch_2)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = ncol(train_par), activation = "linear")

model2 <- keras_model(inputs = c(main_in_2, aux_in_2, feat_in_2), outputs = out_2)
model2 %>% compile(loss = "mse", optimizer = optimizer_adam(learning_rate = 1e-3),
                   metrics = list("mae"))
summary(model2)

history2 <- model2 %>% fit(
  x               = list(train_L, train_N, train_featI),
  y               = train_par,
  epochs          = EPOCHS,
  batch_size      = BATCH_SIZE,
  validation_data = list(list(val_L, val_N, val_featI), val_par),
  callbacks       = make_callbacks(),
  verbose         = 2
)

pred2_std <- predict(model2, list(test_L, test_N, test_featI))
colnames(pred2_std) <- c("mu", "var", "scale")
pred2 <- sweep(sweep(pred2_std, 2, sd_par, `*`), 2, m_par, `+`)
colnames(pred2) <- c("mu", "var", "scale")

# =============================================================================
# 6. METRICS
# =============================================================================

r2_ssr_sst <- function(true_vec, pred_vec) {
  ss_res <- sum((pred_vec - true_vec)^2)
  ss_tot <- sum((true_vec - mean(true_vec))^2)
  1 - ss_res / ss_tot
}

compute_metrics <- function(true_mat, pred_mat, label) {
  purrr::map_dfr(c("mu", "var", "scale"), function(p) {
    e <- pred_mat[, p] - true_mat[, p]
    tibble(
      model = label,
      param = p,
      R2    = round(r2_ssr_sst(true_mat[, p], pred_mat[, p]), 4),
      RMSE  = round(sqrt(mean(e^2)), 4),
      MAE   = round(mean(abs(e)), 4)
    )
  })
}

metrics1 <- compute_metrics(test_par, pred1_std, "CNN base")
metrics2 <- compute_metrics(test_par, pred2_std, "CNN + I-feat")
all_metrics <- bind_rows(metrics1, metrics2)
models_order <- c("CNN base", "CNN + I-feat")

cat("\n===== TEST METRICS (standardized scale) =====\n")
print(all_metrics, n = Inf)

# =============================================================================
# 7. LaTeX TABLE
# =============================================================================

param_latex <- c(mu = "$\\mu$", var = "$\\sigma^2$", scale = "$s$")

tex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{M\\'etricas de recuperaci\\'on de par\\'ametros en escala estandarizada",
  "sobre el conjunto de test — ventana Colombia.}",
  "\\label{tab:metrics_final}",
  "\\begin{tabular}{llrrr}",
  "\\hline",
  "Modelo & Par\\'ametro & $R^2$ & RMSE & MAE \\\\",
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

tex_out <- file.path(fig_dir, "tabla_metricas_final.tex")
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
      labs(x = "True", y = "Estimated", title = param_labels[[p]]) +
      theme_bw(base_size = 10) +
      theme(panel.grid  = element_blank(),
            plot.title  = element_text(hjust = 0.5, size = 11))
  }
  plots[["mu"]] | plots[["var"]] | plots[["scale"]]
}

p_s1 <- make_scatter_grid(true_par, pred1, "CNN base") +
  plot_annotation(title = "CNN base (Vihrs)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

p_s2 <- make_scatter_grid(true_par, pred2, "CNN + I-feat") +
  plot_annotation(title = "CNN + intensity features (propuesta)",
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
         Training = history1$metrics$loss,
         Validation = history1$metrics$val_loss,
         Model = "CNN base"),
  tibble(epoch = seq_along(history2$metrics$loss),
         Training = history2$metrics$loss,
         Validation = history2$metrics$val_loss,
         Model = "CNN + I-feat")
) %>%
  pivot_longer(c(Training, Validation), names_to = "Set", values_to = "value")

p_loss <- ggplot(loss_df, aes(epoch, value, colour = Set, linetype = Model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Training" = "steelblue", "Validation" = "tomato")) +
  labs(x = "Epoch", y = "MSE", colour = "Set", linetype = "Model") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "loss_final_combined.pdf"), p_loss,
       width = 7, height = 4.5, device = cairo_pdf)

model_cnn = subset(loss_df, Model == "CNN base")
p_loss_cnn <- ggplot(model_cnn, aes(epoch, value, colour = Set, linetype = Model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Training" = "steelblue", "Validation" = "tomato")) +
  labs(x = "Epoch", y = "MSE", colour = "Set", linetype = "Model") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "p_loss_cnn.pdf"), p_loss_cnn,
       width = 7, height = 4.5, device = cairo_pdf)

model_cnn_mejora = subset(loss_df, Model == "CNN + I-feat")
p_loss_cnn_mejora <- ggplot(model_cnn_mejora, aes(epoch, value, colour = Set, linetype = Model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Training" = "steelblue", "Validation" = "tomato")) +
  labs(x = "Epoch", y = "MSE", colour = "Set", linetype = "Model") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "p_loss_cnn_mejora.pdf"), p_loss_cnn_mejora,
       width = 7, height = 4.5, device = cairo_pdf)

mae_df <- bind_rows(
  tibble(epoch = seq_along(history1$metrics$mae),
         Training = history1$metrics$mae,
         Validation = history1$metrics$val_mae,
         Model = "CNN base"),
  tibble(epoch = seq_along(history2$metrics$mae),
         Training = history2$metrics$mae,
         Validation = history2$metrics$val_mae,
         Model = "CNN + I-feat")
) %>%
  pivot_longer(c(Training, Validation), names_to = "Set", values_to = "value")

p_mae <- ggplot(mae_df, aes(epoch, value, colour = Set, linetype = Model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Training" = "steelblue", "Validation" = "tomato")) +
  labs(x = "Epoch", y = "MAE", colour = "Set", linetype = "Model") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical")

ggsave(file.path(fig_dir, "mae_final_combined.pdf"), p_mae,
       width = 7, height = 4.5, device = cairo_pdf)

cat("Loss and MAE plots saved\n")

# =============================================================================
# 10. R² COMPARISON BAR CHART
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
  scale_fill_manual(values = c("CNN base"   = "grey65",
                               "CNN + I-feat" = "darkorange")) +
  labs(x = "Parameter", y = expression(R^2), fill = "Model") +
  ylim(0, 1.08) +
  theme_bw(base_size = 11) +
  theme(panel.grid      = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "r2_comparison_final.pdf"), p_r2,
       width = 5, height = 4, device = cairo_pdf)

cat("R-squared comparison saved\n")
cat("\nDone! All figures and table in:", fig_dir, "\n")
