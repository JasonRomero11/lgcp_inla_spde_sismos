# =============================================================================
# replot_figures.R  --  Regenera TODAS las figuras y tablas SIN reentrenar
# =============================================================================
#
# POR QUE EXISTE ESTE SCRIPT
#   El entrenamiento Keras/TF y las simulaciones rLGCP no son 100% reproducibles
#   entre corridas (ver README / seccion de semillas). Reentrenar solo para
#   cambiar una etiqueta CAMBIA los numeros de las figuras y de las tablas.
#   Este script regenera todo a partir de los objetos ya calculados, asi que
#   los numeros son EXACTAMENTE los de la corrida original.
#
# MODO DE USO
#   (A) Con la sesion de R viva (los objetos siguen en memoria):
#         source("scripts/ENTRENANDO_CNN/replot_figures.R")
#       -> congela los objetos en figures/figure_objects.rds y regenera todo.
#
#   (B) En una sesion nueva, despues de haber corrido (A) al menos una vez:
#         source("scripts/ENTRENANDO_CNN/replot_figures.R")
#       -> lee figures/figure_objects.rds y regenera todo. No necesita Keras,
#          ni GPU, ni los datos de simulacion.
#
#   El codigo de graficado es copia literal de las secciones 6-10 de
#   CNN_train_and_predict.R. Si editas etiquetas alli, vuelve a extraerlas aqui.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork)
})

setwd("/home/jasonromeroia/Documents/personal/paper_lgcp_features_computers-and-geosciences/lgcp-cnn-features")
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

RDS_OBJ <- file.path(fig_dir, "figure_objects.rds")

# Objetos minimos que necesitan las figuras (nada de modelos Keras ni sims_val)
NEEDED <- c(
  # secciones 6-10 (conjunto de prueba)
  "test_par", "pred1_std", "pred2_std", "true_par", "pred1", "pred2",
  "history1", "history2",
  # secciones 6a-6c (catalogo observado 2020)
  "r_obs", "Lc_obs", "L_med", "L_lo", "L_hi",
  "N_obs", "N_sims", "pred1_obs", "pred2_obs", "EN_base", "EN_feat"
)

in_memory <- vapply(NEEDED, exists, logical(1), envir = globalenv())

# --- PASO 1: congelar SIEMPRE lo que haya en memoria -------------------------
# Se hace ANTES de graficar: si algo falla despues, los objetos ya estan a salvo.
# Se fusiona con el rds previo, asi que una sesion parcial nunca borra lo ya
# congelado en una sesion anterior.
if (any(in_memory)) {
  obj <- mget(NEEDED[in_memory], envir = globalenv())
  # history: guardar solo $metrics (lista R pura, sin punteros de Python)
  if ("history1" %in% names(obj)) obj$history1 <- list(metrics = history1$metrics)
  if ("history2" %in% names(obj)) obj$history2 <- list(metrics = history2$metrics)
  if (file.exists(RDS_OBJ)) {
    prev <- readRDS(RDS_OBJ)
    prev[names(obj)] <- obj          # lo vivo gana sobre lo guardado
    obj <- prev
  }
  saveRDS(obj, RDS_OBJ)
  message("Congelados ", length(obj), "/", length(NEEDED), " objetos en ", RDS_OBJ,
          " (", round(file.size(RDS_OBJ) / 1024), " KB)")
} else {
  message("Ningun objeto en memoria; se usara ", RDS_OBJ)
}

# --- PASO 2: completar lo que falte desde el rds ------------------------------
if (!all(in_memory) && file.exists(RDS_OBJ)) {
  faltan <- NEEDED[!in_memory]
  saved  <- readRDS(RDS_OBJ)
  hay    <- intersect(faltan, names(saved))
  if (length(hay)) {
    message("Recuperando del rds: ", paste(hay, collapse = ", "))
    invisible(list2env(saved[hay], envir = globalenv()))
  }
}

# --- PASO 3: verificar antes de graficar --------------------------------------
missing <- NEEDED[!vapply(NEEDED, exists, logical(1), envir = globalenv())]
if (length(missing)) {
  stop("Faltan objetos para graficar: ", paste(missing, collapse = ", "), "\n",
       "  Los que si estaban ya quedaron congelados en ", RDS_OBJ, ".\n",
       "  Corre CNN_train_and_predict.R y vuelve a ejecutar este script ",
       "en la MISMA sesion.")
}
message("OK: los ", length(NEEDED), " objetos necesarios estan disponibles.\n")

# ############################################################################
# COPIA LITERAL: secciones 6-10 de CNN_train_and_predict.R
# ############################################################################

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

# ############################################################################
# COPIA LITERAL: seccion 6 (figuras del catalogo observado) de CNN_train_and_predict.R
# ############################################################################


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



cat("\n>>> Figuras regeneradas en:", normalizePath(fig_dir), "\n")
