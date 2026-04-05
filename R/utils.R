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
# Cargar el paquete
library(fmesher)
library(xtable)
if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable") # Instalar desde el repositorio de INLA
}
library(patchwork)




create_mesh = function(shapeZona, xy, 
                       param_cutoff = 10000, offset_param = c(100, 20000)){
  
  shapeZona <- st_make_valid(shapeZona)   
  shapeZona <- as_Spatial(shapeZona)   
  shapeZona <- st_as_sf(shapeZona)
  
  max.edge_params = max(c(diff(range(xy[,1])), diff(range(xy[,2]))))/30
  
  meshSismos = inla.mesh.2d(loc.domain = shapeZona,
                            loc = cbind(xy),
                            max.edge = c(1, 3) * max.edge_params,
                            cutoff = param_cutoff,
                            offset = offset_param,
                            crs = st_crs(shapeZona))
  
  return(meshSismos)
}


create_ppp = function(shapeZona, sismosdf, save_plot = FALSE, file_path = NULL){

  shapeOwin <- as.owin(shapeZona)
  xy <- st_coordinates(sismosdf)
  p <- ppp(xy[,1], xy[,2], window = shapeOwin)
  unitname(p) = c("meter", "meters")
  

  if (save_plot) {
    if (is.null(file_path)) {
      stop("You must specify 'file_path' when save_plot = TRUE.")
    }
    
    ext <- tools::file_ext(file_path)
    
    if (ext == "png") {
      png(file_path, width = 1200, height = 1200, res = 150)
    } else if (ext == "pdf") {
      pdf(file_path, width = 8, height = 8)
    } else if (ext == "jpeg" || ext == "jpg") {
      jpeg(file_path, width = 1200, height = 1200, res = 150)
    } else {
      stop("Unsupported file format. Use .png, .pdf, .jpg")
    }
  }
  
  par(mai = c(0,0,0,0))
  plot(p)
  title("Eventos Sísmicos Colombia")
  
  # Close device if saving
  if (save_plot) {
    dev.off()
    message("Plot saved to: ", file_path)
  }
  
  return(list(point_pattern = p, coordinates_eventos = xy))
}


calc_model_fit <- function(O, E) {
  # Seguridad numérica
  E <- pmax(E, 1e-10)
  
  # R2 clásico (basado en sumas de cuadrados)
  R2_model <- 1 - sum((O - E)^2, na.rm = TRUE) /
    sum((O - mean(O, na.rm = TRUE))^2, na.rm = TRUE)
  
  # Deviance del modelo
  D_model <- sum(2 * (O * log(pmax(O / E, 1e-10)) - (O - E)), na.rm = TRUE)
  
  # Deviance del modelo nulo
  mean_O <- mean(O, na.rm = TRUE)
  D_null <- sum(2 * (O * log(pmax(O / mean_O, 1e-10)) - (O - mean_O)), na.rm = TRUE)
  
  R2_deviance <- 1 - D_model / D_null
  
  # Log-verosimilitudes
  logLik_model <- sum(dpois(O, lambda = E, log = TRUE), na.rm = TRUE)
  logLik_null  <- sum(dpois(O, lambda = mean_O, log = TRUE), na.rm = TRUE)
  
  R2_McFadden <- 1 - logLik_model / logLik_null
  
  # R² tipo correlación
  R2_pearson <- cor(O, E, use = "complete.obs")^2
  
  # Retornar como lista
  return(list(
    R2_model     = R2_model,
    R2_deviance  = R2_deviance,
    R2_McFadden  = R2_McFadden,
    R2_pearson   = R2_pearson
  ))
}


plot_mean_spatial_effect <- function(result_model,
                                     output_path = NULL,
                                     width = 8,
                                     height = 5,
                                     dpi = 300) {
  
  preds <- dmesh %>%
    st_as_sf() %>%
    dplyr::mutate(
      posterior = result_model$summary.random$spatial.field$mean
    )
  
  p <- ggplot(preds) +
    theme_bw() +
    geom_sf(aes(fill = posterior), colour = NA) +
    
    # 🔹 MUY IMPORTANTE: elimina padding espacial
    coord_sf(expand = FALSE) +
    
    scale_fill_gradientn(
      name = "Estimación media posterior del campo espacial",
      colours = viridis(10),
      breaks = seq(-7, 7, 2),
      guide = guide_colourbar(
        nbin = 500,
        raster = TRUE,
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 1,
        barwidth = 20,
        barheight = 1,
        direction = "horizontal",
        title.position = "top",
        title.theme = element_text(hjust = 0.5)
      )
    ) +
    
    geom_sf(
      data = st_as_sf(shapeZona_sp),
      fill = NA,
      colour = "white",
      linewidth = 0.2
    ) +
    
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      
      # 🔹 Reduce espacio en blanco externo
      plot.margin = margin(2, 2, 2, 2)
    )
  
  if (!is.null(output_path)) {
    ggsave(
      filename = output_path,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(p)
}





results_model_to_tableLatex = function(result_model){
  table_results <- result_model$summary.fixed[, c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")]
  table_results <- round(table_results, 3)
  table_results <- cbind(Variable = rownames(table_results), table_results)
  tabla_latex <- xtable(table_results, 
                        caption = "Resumen de los efectos fijos del modelo",
                        label = "tab:fixed_effects",
                        digits = 3)
  print(tabla_latex, 
        include.rownames = FALSE, 
        caption.placement = "top",
        booktabs = TRUE,  
        sanitize.text.function = function(x){x})
}



plot_significance_result <- function(
    result_model,
    output_path = NULL,
    width = 8,
    height = 5,
    dpi = 300
) {
  
  summary_data <- result_model$summary.fixed
  summary_data <- summary_data[-1, ]  # quitar intercepto
  
  plot_data <- data.frame(
    Variable = rownames(summary_data),
    Mean = summary_data$mean,
    Lower = summary_data$`0.025quant`,
    Upper = summary_data$`0.975quant`
  )
  
  plot_data$Significativo <- ifelse(
    plot_data$Lower > 0 | plot_data$Upper < 0,
    "Significativo",
    "No significativo"
  )
  
  p <- ggplot(plot_data,
              aes(x = reorder(Variable, Mean), y = Mean)) +
    
    ## Intervalo “grueso” (mejor que línea fina)
    geom_linerange(
      aes(ymin = Lower, ymax = Upper, color = Significativo),
      linewidth = 2
    ) +
    
    ## Punto central
    geom_point(size = 3) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "red") +
    
    coord_flip() +
    
    labs(
      title = "",
      x = "Covariables",
      y = "Estimación posterior",
      color = "Significancia"
    ) +
    
    scale_color_manual(
      values = c(
        "Significativo" = "#1f78b4",
        "No significativo" = "gray60"
      )
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  ## Guardar si se especifica ruta
  if (!is.null(output_path)) {
    ggsave(
      filename = output_path,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}


plot_spatial_effects_2x2 <- function(resM0, resM1, resM2, resM3,
                                     titles = c("M0", "M1", "M2", "M3"),
                                     output_path = NULL,
                                     width = 12,
                                     height = 8,
                                     dpi = 300,
                                     breaks = seq(-7, 7, 2),
                                     limits = NULL,           # si NULL => calcula límites globales
                                     use_quantile_limits = FALSE,
                                     q = c(0.01, 0.99)) {
  
  stopifnot(length(titles) == 4)
  
  # ---- helper para extraer el vector de medias del campo espacial
  get_mean_field <- function(res) {
    v <- res$summary.random$spatial.field$mean
    if (is.null(v)) stop("No encuentro summary.random$spatial.field$mean en uno de los modelos.")
    as.numeric(v)
  }
  
  v0 <- get_mean_field(resM0)
  v1 <- get_mean_field(resM1)
  v2 <- get_mean_field(resM2)
  v3 <- get_mean_field(resM3)
  
  allv <- c(v0, v1, v2, v3)
  allv <- allv[is.finite(allv)]
  
  # ---- límites comunes (misma escala para todos)
  if (is.null(limits)) {
    if (use_quantile_limits) {
      lim <- as.numeric(stats::quantile(allv, probs = q, na.rm = TRUE))
      # asegura simetría opcional si te interesa (comenta si no)
      # m <- max(abs(lim)); lim <- c(-m, m)
      limits <- lim
    } else {
      limits <- range(allv, na.rm = TRUE)
    }
  }
  
  # Si el usuario da breaks tipo seq(-7,7,2), pero limits no cubre eso,
  # puedes dejar breaks tal cual o recalcularlos:
  # breaks <- pretty(limits, n = 7)
  
  # ---- helper: construir un ggplot para un modelo
  make_plot <- function(res, title_txt) {
    vals <- get_mean_field(res)
    
    preds <- dmesh |>
      sf::st_as_sf() |>
      dplyr::mutate(posterior = vals)
    
    ggplot2::ggplot(preds) +
      ggplot2::geom_sf(ggplot2::aes(fill = posterior), colour = NA) +
      ggplot2::geom_sf(
        data = sf::st_as_sf(shapeZona_sp),
        fill = NA,
        colour = "white",
        linewidth = 0.2
      ) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        name = "Estimación media posterior del campo espacial",
        colours = viridis::viridis(10),
        limits = limits,
        breaks = breaks,
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(
          nbin = 500,
          raster = TRUE,
          frame.colour = "black",
          ticks.colour = "black",
          frame.linewidth = 1,
          barwidth = 20,
          barheight = 1,
          direction = "horizontal",
          title.position = "top",
          title.theme = ggplot2::element_text(hjust = 0.5)
        )
      ) +
      ggplot2::labs(title = title_txt) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
  }
  
  p0 <- make_plot(resM0, titles[1])
  p1 <- make_plot(resM1, titles[2])
  p2 <- make_plot(resM2, titles[3])
  p3 <- make_plot(resM3, titles[4])
  
  # ---- armar 2x2 y dejar UNA sola leyenda común
  # install.packages("patchwork") si no lo tienes
  spatial_2x2 <- (p0 | p1) / (p2 | p3) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        legend.position = "bottom"
      )
    )
  
  
  if (!is.null(output_path)) {
    ggplot2::ggsave(
      filename = output_path,
      plot = spatial_2x2,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(list(
    plot = spatial_2x2,
    limits = limits
  ))
}

plot_spatial_effects_1x2 <- function(resA, resB,
                                     titles = c("Modelo A", "Modelo B"),
                                     output_path = NULL,
                                     width = 11,
                                     height = 5,
                                     dpi = 300,
                                     breaks = seq(-7, 7, 2),
                                     limits = NULL,
                                     use_quantile_limits = FALSE,
                                     q = c(0.01, 0.99),
                                     show_stats = TRUE,
                                     stats_position = "topleft",
                                     stats_size = 2.8) {
  
  stopifnot(length(titles) == 2)
  
  # ---- helper para extraer el vector de medias del campo espacial
  get_mean_field <- function(res) {
    v <- res$summary.random$spatial.field$mean
    if (is.null(v)) {
      stop("No encuentro summary.random$spatial.field$mean en uno de los modelos.")
    }
    as.numeric(v)
  }
  
  vA <- get_mean_field(resA)
  vB <- get_mean_field(resB)
  
  allv <- c(vA, vB)
  allv <- allv[is.finite(allv)]
  
  # ---- límites comunes (misma escala)
  if (is.null(limits)) {
    if (use_quantile_limits) {
      limits <- as.numeric(stats::quantile(allv, probs = q, na.rm = TRUE))
    } else {
      limits <- range(allv, na.rm = TRUE)
    }
  }
  
  # ---- helper: calcular estadísticos
  calc_stats <- function(v) {
    v <- v[is.finite(v)]
    data.frame(
      Min = min(v, na.rm = TRUE),
      Mediana = median(v, na.rm = TRUE),
      Media = mean(v, na.rm = TRUE),
      Max = max(v, na.rm = TRUE)
    )
  }
  
  # ---- helper: crear etiqueta de estadísticos
  make_stats_label <- function(v) {
    s <- calc_stats(v)
    paste0(
      "Min: ", sprintf("%.2f", s$Min), "\n",
      "Mediana: ", sprintf("%.2f", s$Mediana), "\n",
      "Media: ", sprintf("%.2f", s$Media), "\n",
      "Max: ", sprintf("%.2f", s$Max)
    )
  }
  
  # ---- helper: obtener posición del cuadro de estadísticos
  get_stats_position <- function(position, bbox) {
    x_range <- as.numeric(bbox["xmax"] - bbox["xmin"])
    y_range <- as.numeric(bbox["ymax"] - bbox["ymin"])
    
    positions <- list(
      "topright" = c(x = as.numeric(bbox["xmax"]) - 0.05 * x_range, 
                     y = as.numeric(bbox["ymax"]) - 0.05 * y_range,
                     hjust = 1, vjust = 1),
      "topleft" = c(x = as.numeric(bbox["xmin"]) + 0.05 * x_range, 
                    y = as.numeric(bbox["ymax"]) - 0.05 * y_range,
                    hjust = 0, vjust = 1),
      "bottomright" = c(x = as.numeric(bbox["xmax"]) - 0.05 * x_range, 
                        y = as.numeric(bbox["ymin"]) + 0.05 * y_range,
                        hjust = 1, vjust = 0),
      "bottomleft" = c(x = as.numeric(bbox["xmin"]) + 0.05 * x_range, 
                       y = as.numeric(bbox["ymin"]) + 0.05 * y_range,
                       hjust = 0, vjust = 0)
    )
    positions[[position]]
  }
  
  # ---- helper: construir un ggplot para un modelo
  make_plot <- function(res, title_txt, field_values) {
    
    preds <- dmesh |>
      sf::st_as_sf() |>
      dplyr::mutate(
        posterior = get_mean_field(res)
      )
    
    # Obtener bbox para posicionar estadísticos
    bbox <- sf::st_bbox(preds)
    pos <- get_stats_position(stats_position, bbox)
    
    p <- ggplot2::ggplot(preds) +
      ggplot2::geom_sf(ggplot2::aes(fill = posterior), colour = NA) +
      ggplot2::geom_sf(
        data = sf::st_as_sf(shapeZona_sp),
        fill = NA,
        colour = "white",
        linewidth = 0.2
      ) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        name = "Estimación media\nposterior del\ncampo espacial",
        colours = viridis::viridis(10),
        limits = limits,
        breaks = breaks,
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(
          nbin = 500,
          raster = TRUE,
          frame.colour = "black",
          ticks.colour = "black",
          frame.linewidth = 0.5,
          barwidth = 1.2,
          barheight = 15,
          direction = "vertical",
          title.position = "top",
          title.theme = ggplot2::element_text(hjust = 0.5, size = 9)
        )
      ) +
      ggplot2::labs(title = title_txt) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
    
    # Agregar cuadro de estadísticos si show_stats = TRUE
    if (show_stats) {
      stats_label <- make_stats_label(field_values)
      
      # Crear dataframe para geom_label (más control que annotate)
      stats_df <- data.frame(
        x = as.numeric(pos["x"]),
        y = as.numeric(pos["y"]),
        label = stats_label
      )
      
      p <- p +
        ggplot2::geom_label(
          data = stats_df,
          ggplot2::aes(x = x, y = y, label = label),
          hjust = as.numeric(pos["hjust"]),
          vjust = as.numeric(pos["vjust"]),
          size = stats_size,
          fill = ggplot2::alpha("white", 0.85),
          colour = "black",
          label.padding = ggplot2::unit(0.4, "lines"),
          label.r = ggplot2::unit(0.15, "lines"),
          family = "mono",
          inherit.aes = FALSE
        )
    }
    
    return(p)
  }
  
  pA <- make_plot(resA, titles[1], vA)
  pB <- make_plot(resB, titles[2], vB)
  
  # ---- armar 1x2 con una sola leyenda
  spatial_1x2 <- (pA | pB) +
    patchwork::plot_layout(guides = "collect")
  
  if (!is.null(output_path)) {
    ggplot2::ggsave(
      filename = output_path,
      plot = spatial_1x2,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  # ---- Crear tabla de estadísticos para retornar
  stats_table <- data.frame(
    Modelo = titles,
    Min = c(calc_stats(vA)$Min, calc_stats(vB)$Min),
    Mediana = c(calc_stats(vA)$Mediana, calc_stats(vB)$Mediana),
    Media = c(calc_stats(vA)$Media, calc_stats(vB)$Media),
    Max = c(calc_stats(vA)$Max, calc_stats(vB)$Max)
  )
  
  return(list(
    plot = spatial_1x2,
    limits = limits,
    stats = stats_table
  ))
}



tema_tesis <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "#f8fafc", color = NA),
      panel.background = element_rect(fill = "#ffffff", color = NA),
      panel.grid.major = element_line(color = "#e2e8f0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", color = "#1e293b", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "#64748b", hjust = 0.5),
      axis.title = element_text(size = 10, color = "#475569"),
      axis.text = element_text(size = 9, color = "#64748b"),
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      strip.text = element_text(size = 10, face = "bold", color = "#1e293b")
    )
}

# ============================================================================
# FUNCIÓN 1: Gráfico de barras comparativo Log-Score y LCPO
# ============================================================================

grafico_comparacion_metricas <- function(LS_M0, LS_M1, LS_M2, LS_M3,
                                         LCPO_M0, LCPO_M1, LCPO_M2, LCPO_M3,
                                         titulo = "") {
  
  # Crear dataframe
  datos <- data.frame(
    Modelo = factor(c("M0", "M1", "M2", "M3"), levels = c("M0", "M1", "M2", "M3")),
    Prior = c("Sin Prior", "Con Prior", "Sin Prior", "Con Prior"),
    Covariables = c("M0 - M1", "M0 - M1", "M2 - M3", "M2 - M3"),
    LogScore = c(LS_M0, LS_M1, LS_M2, LS_M3),
    LCPO = c(LCPO_M0, LCPO_M1, LCPO_M2, LCPO_M3)
  )
  
  # Colores académicos: gris para sin prior, azul moderado para con prior
  colores_prior <- c("Sin Prior" = "#7D8A96", "Con Prior" = "#4682B4")
  
  # Gráfico Log-Score
  p1 <- ggplot(datos, aes(x = Modelo, y = LogScore, fill = Prior)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.5, alpha = 0.85) +
    geom_text(aes(label = round(LogScore, 4)), 
              vjust = -0.5, size = 3, color = "#2d2d2d", fontface = "bold") +
    scale_fill_manual(values = colores_prior) +
    facet_wrap(~Covariables, scales = "free_x") +
    labs(
      title = "Log-Score",
      x = NULL,
      y = "Log-Score",
      fill = "Prior Informativo"
    ) +
    tema_tesis() +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(min(datos$LogScore) - 2, max(datos$LogScore) + 2))
  
  # Gráfico LCPO
  p2 <- ggplot(datos, aes(x = Modelo, y = LCPO/1000, fill = Prior)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.5, alpha = 0.85) +
    geom_text(aes(label = round(LCPO, 2)), 
              vjust = -0.5, size = 3, color = "#2d2d2d", fontface = "bold") +
    scale_fill_manual(values = colores_prior) +
    facet_wrap(~Covariables, scales = "free_x") +
    labs(
      title = "LCPO",
      x = NULL,
      y = "LCPO (miles)",
      fill = "Prior Informativo"
    ) +
    tema_tesis() +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(min(datos$LCPO/1000) - 2, max(datos$LCPO/1000) + 2))
  
  # Leyenda común
  leyenda <- ggplot(datos, aes(x = Modelo, y = LogScore, fill = Prior)) +
    geom_col() +
    scale_fill_manual(values = colores_prior) +
    tema_tesis() +
    theme(legend.position = "bottom")
  
  leyenda_grob <- cowplot::get_legend(leyenda)
  
  # Combinar
  grafico_final <- (p1 | p2) / 
    cowplot::ggdraw(leyenda_grob) +
    plot_layout(heights = c(10, 1)) +
    plot_annotation(
      title = titulo,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#1a1a1a")
      )
    )
  
  return(grafico_final)
}



grafico_bootstrap <- function(delta_vector, 
                              modelo_base = "M2", 
                              modelo_comparado = "M3",
                              titulo = NULL,
                              stats_position = "topright",
                              stats_size = 3,
                              output_path = NULL,
                              width = 10,
                              height = 6,
                              dpi = 300) {
  
  # Calcular estadísticos
  media <- mean(delta_vector)
  q025 <- quantile(delta_vector, 0.025)
  q975 <- quantile(delta_vector, 0.975)
  p_mejor <- mean(delta_vector > 0)
  B <- length(delta_vector)
  
  # Crear dataframe
  datos_boot <- data.frame(delta = delta_vector)
  
  # Colores consistentes con tema de tesis
  color_principal <- "#4682B4"
  color_gris <- "#7D8A96"
  color_ic <- "#5A8FAD"
  color_referencia <- "#2d2d2d"
  
  # Crear etiqueta de estadísticos
  stats_text <- paste0(
    "ESTADÍSTICOS BOOTSTRAP\n",
    "Media: ", sprintf("%.4f", media), "\n",
    "IC 95%: [", sprintf("%.4f", q025), ", ", sprintf("%.4f", q975), "]\n",
    "P(", modelo_comparado, " > ", modelo_base, "): ", sprintf("%.1f%%", p_mejor * 100), "\n",
    ifelse(p_mejor > 0.95, 
           paste0(modelo_comparado, " signif. mejor"),
           ifelse(p_mejor < 0.05,
                  paste0(modelo_base, " signif. mejor"),
                  "Sin diferencia signif."))
  )
  
  # Calcular posición del cuadro de estadísticos
  x_range <- range(delta_vector)
  y_dens <- density(delta_vector)
  y_max <- max(y_dens$y) * 1.15
  
  if (stats_position == "topright") {
    x_pos <- x_range[2] - 0.02 * diff(x_range)
    y_pos <- y_max * 0.95
    hjust_val <- 1
    vjust_val <- 1
  } else if (stats_position == "topleft") {
    x_pos <- x_range[1] + 0.02 * diff(x_range)
    y_pos <- y_max * 0.95
    hjust_val <- 0
    vjust_val <- 1
  } else if (stats_position == "bottomright") {
    x_pos <- x_range[2] - 0.02 * diff(x_range)
    y_pos <- y_max * 0.05
    hjust_val <- 1
    vjust_val <- 0
  } else {
    x_pos <- x_range[1] + 0.02 * diff(x_range)
    y_pos <- y_max * 0.05
    hjust_val <- 0
    vjust_val <- 0
  }
  
  # Gráfico principal (histograma + densidad + estadísticos)
  p_hist <- ggplot(datos_boot, aes(x = delta)) +
    # Área del IC
    annotate("rect", xmin = q025, xmax = q975, ymin = -Inf, ymax = Inf,
             fill = color_ic, alpha = 0.15) +
    # Histograma
    geom_histogram(aes(y = after_stat(density)), 
                   bins = 40, 
                   fill = color_principal, 
                   color = "white",
                   alpha = 0.7) +
    # Curva de densidad
    geom_density(color = color_principal, linewidth = 1.2, fill = color_principal, alpha = 0.2) +
    # Línea de referencia en 0
    geom_vline(xintercept = 0, color = color_referencia, linewidth = 1.2, linetype = "solid") +
    # Línea de la media
    geom_vline(xintercept = media, color = color_gris, linewidth = 0.8, linetype = "dashed", alpha = 0.8) +
    # IC 95%
    geom_vline(xintercept = q025, color = color_ic, linewidth = 1, linetype = "dotted") +
    geom_vline(xintercept = q975, color = color_ic, linewidth = 1, linetype = "dotted") +
    # Anotación de la media
    annotate("text", x = media, y = y_max * 0.98, 
             label = paste0("μ = ", round(media, 4)),
             hjust = 0.5, color = color_referencia, fontface = "bold", size = 3.5) +
    # Anotación del 0
    annotate("text", x = 0, y = y_max * 0.98, 
             label = "0",
             hjust = 0.5, color = color_referencia, fontface = "bold", size = 3.5) +
    # Cuadro de estadísticos
    annotate("label", 
             x = x_pos, 
             y = y_pos, 
             label = stats_text,
             hjust = hjust_val, 
             vjust = vjust_val,
             size = stats_size,
             fill = alpha("white", 0.9),
             colour = color_referencia,
             label.padding = unit(0.5, "lines"),
             label.r = unit(0.15, "lines"),
             family = "mono") +
    labs(
      title = titulo,
      x = expression(Delta ~ "Log-Score"),
      y = "Densidad"
    ) +
    tema_tesis() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = color_gris),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11)
    ) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")
  
  # Guardar si se especifica ruta
  if (!is.null(output_path)) {
    ggsave(
      filename = output_path,
      plot = p_hist,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    cat("Gráfico guardado en:", output_path, "\n")
  }
  
  # Retornar gráfico y estadísticos
  return(list(
    plot = p_hist,
    stats = data.frame(
      Media = media,
      IC_inferior = as.numeric(q025),
      IC_superior = as.numeric(q975),
      P_mejor = p_mejor,
      B = B
    )
  ))
}


plot_intensidades_modelos <- function(spdf_sf,
                                      cols_intensidad = c("expect_M0", "expect_M1", 
                                                          "expect_M2", "expect_M3"),
                                      modelos = c("M0", "M1", "M2", "M3"),
                                      shapeZona = NULL,
                                      output_path = NULL,
                                      width = 14,
                                      height = 5,
                                      dpi = 300,
                                      breaks = NULL,
                                      limits = NULL,
                                      use_quantile_limits = TRUE,
                                      q = c(0.01, 0.99),
                                      show_stats = TRUE,
                                      stats_position = "topleft",
                                      stats_size = 2.5,
                                      ncol = 4,
                                      palette = "viridis",
                                      legend_title = "Intensidad\nesperada (λ)") {
  
  stopifnot(length(cols_intensidad) == length(modelos))
  
  # Convertir a sf si es necesario
  if (!inherits(spdf_sf, "sf")) {
    spdf_sf <- sf::st_as_sf(spdf_sf)
  }
  
  # ---- Calcular límites comunes
  all_values <- unlist(lapply(cols_intensidad, function(col) spdf_sf[[col]]))
  all_values <- all_values[is.finite(all_values)]
  
  if (is.null(limits)) {
    if (use_quantile_limits) {
      limits <- as.numeric(quantile(all_values, probs = q, na.rm = TRUE))
    } else {
      limits <- range(all_values, na.rm = TRUE)
    }
  }
  
  if (is.null(breaks)) {
    breaks <- pretty(limits, n = 5)
  }
  
  # ---- Helper: calcular estadísticos
  calc_stats <- function(v) {
    v <- v[is.finite(v)]
    data.frame(
      Min = min(v, na.rm = TRUE),
      Mediana = median(v, na.rm = TRUE),
      Media = mean(v, na.rm = TRUE),
      Max = max(v, na.rm = TRUE)
    )
  }
  
  # ---- Helper: crear etiqueta de estadísticos
  make_stats_label <- function(v) {
    s <- calc_stats(v)
    paste0(
      "Min: ", sprintf("%.2f", s$Min), "\n",
      "Med: ", sprintf("%.2f", s$Mediana), "\n",
      "μ: ", sprintf("%.2f", s$Media), "\n",
      "Max: ", sprintf("%.2f", s$Max)
    )
  }
  
  # ---- Helper: obtener posición del cuadro de estadísticos
  get_stats_position <- function(position, bbox) {
    x_range <- as.numeric(bbox["xmax"] - bbox["xmin"])
    y_range <- as.numeric(bbox["ymax"] - bbox["ymin"])
    
    positions <- list(
      "topright" = c(x = as.numeric(bbox["xmax"]) - 0.05 * x_range, 
                     y = as.numeric(bbox["ymax"]) - 0.05 * y_range,
                     hjust = 1, vjust = 1),
      "topleft" = c(x = as.numeric(bbox["xmin"]) + 0.05 * x_range, 
                    y = as.numeric(bbox["ymax"]) - 0.05 * y_range,
                    hjust = 0, vjust = 1),
      "bottomright" = c(x = as.numeric(bbox["xmax"]) - 0.05 * x_range, 
                        y = as.numeric(bbox["ymin"]) + 0.15 * y_range,
                        hjust = 1, vjust = 0),
      "bottomleft" = c(x = as.numeric(bbox["xmin"]) + 0.05 * x_range, 
                       y = as.numeric(bbox["ymin"]) + 0.15 * y_range,
                       hjust = 0, vjust = 0)
    )
    positions[[position]]
  }
  
  # ---- Seleccionar paleta de colores
  color_palette <- switch(palette,
                          "viridis" = viridis::viridis(10),
                          "plasma" = viridis::plasma(10),
                          "inferno" = viridis::inferno(10),
                          "magma" = viridis::magma(10),
                          "cividis" = viridis::cividis(10),
                          viridis::viridis(10))
  
  # ---- Helper: construir un ggplot para un modelo
  make_plot <- function(col_name, title_txt) {
    
    preds <- spdf_sf %>%
      dplyr::mutate(intensidad = .data[[col_name]])
    
    # Obtener bbox para posicionar estadísticos
    bbox <- sf::st_bbox(preds)
    pos <- get_stats_position(stats_position, bbox)
    
    p <- ggplot2::ggplot(preds) +
      ggplot2::geom_sf(ggplot2::aes(fill = intensidad), colour = NA) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        name = legend_title,
        colours = color_palette,
        limits = limits,
        breaks = breaks,
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(
          nbin = 500,
          raster = TRUE,
          frame.colour = "black",
          ticks.colour = "black",
          frame.linewidth = 0.5,
          barwidth = 1.2,
          barheight = 12,
          direction = "vertical",
          title.position = "top",
          title.theme = ggplot2::element_text(hjust = 0.5, size = 9)
        )
      ) +
      ggplot2::labs(title = title_txt) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      )
    
    # Agregar shape de referencia si se proporciona
    if (!is.null(shapeZona)) {
      p <- p +
        ggplot2::geom_sf(
          data = sf::st_as_sf(shapeZona),
          fill = NA,
          colour = "white",
          linewidth = 0.3
        )
    }
    
    # Agregar cuadro de estadísticos si show_stats = TRUE
    if (show_stats) {
      field_values <- preds$intensidad
      stats_label <- make_stats_label(field_values)
      
      stats_df <- data.frame(
        x = as.numeric(pos["x"]),
        y = as.numeric(pos["y"]),
        label = stats_label
      )
      
      p <- p +
        ggplot2::geom_label(
          data = stats_df,
          ggplot2::aes(x = x, y = y, label = label),
          hjust = as.numeric(pos["hjust"]),
          vjust = as.numeric(pos["vjust"]),
          size = stats_size,
          fill = ggplot2::alpha("white", 0.85),
          colour = "black",
          label.padding = ggplot2::unit(0.3, "lines"),
          label.r = ggplot2::unit(0.15, "lines"),
          family = "mono",
          inherit.aes = FALSE
        )
    }
    
    return(p)
  }
  
  # ---- Crear lista de gráficos
  plots_list <- mapply(make_plot, 
                       cols_intensidad, 
                       modelos, 
                       SIMPLIFY = FALSE)
  
  # ---- Combinar gráficos
  if (ncol == 4) {
    combined_plot <- (plots_list[[1]] | plots_list[[2]] | plots_list[[3]] | plots_list[[4]]) +
      patchwork::plot_layout(guides = "collect")
  } else if (ncol == 2) {
    combined_plot <- (plots_list[[1]] | plots_list[[2]]) / 
      (plots_list[[3]] | plots_list[[4]]) +
      patchwork::plot_layout(guides = "collect")
  } else {
    combined_plot <- patchwork::wrap_plots(plots_list, ncol = ncol) +
      patchwork::plot_layout(guides = "collect")
  }
  
  # ---- Guardar si se especifica ruta
  if (!is.null(output_path)) {
    ggplot2::ggsave(
      filename = output_path,
      plot = combined_plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    cat("Gráfico guardado en:", output_path, "\n")
  }
  
  # ---- Crear tabla de estadísticos
  stats_table <- do.call(rbind, lapply(seq_along(modelos), function(i) {
    s <- calc_stats(spdf_sf[[cols_intensidad[i]]])
    data.frame(Modelo = modelos[i], s)
  }))
  
  return(list(
    plot = combined_plot,
    limits = limits,
    stats = stats_table
  ))
}

compute_oe_metrics <- function(observed, expected) {
  # Comprobaciones básicas
  if (length(observed) != length(expected)) {
    stop("observed y expected deben tener la misma longitud")
  }
  
  # Eliminar pares con NA
  valid <- is.finite(observed) & is.finite(expected)
  O <- observed[valid]
  E <- expected[valid]
  
  m <- length(O)
  if (m == 0) {
    stop("No hay observaciones válidas tras eliminar NA/Inf")
  }
  
  # Métricas de error
  rmse <- sqrt(mean((O - E)^2))
  mae  <- mean(abs(O - E))
  
  # R^2 exploratorio Observado–Esperado
  O_bar <- mean(O)
  r2 <- 1 - sum((O - E)^2) / sum((O - O_bar)^2)
  
  return(list(
    RMSE = rmse,
    MAE  = mae,
    R2   = r2,
    n    = m
  ))
}



#if(use_barrier){
#  tl = length(meshSismos$graph$tv[,1]) 
#  posTri = matrix(0, tl, 2) 
#  for (t in 1:tl){
#    temp = meshSismos$loc[meshSismos$graph$tv[t, ], ]
#    posTri[t,] = colMeans(temp)[c(1,2)] 
#  }
#  prj = "+proj=tmerc +lat_0=4.59620041666667 +lon_0=-74.0775079166667 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
#  posTri = SpatialPoints(posTri, proj4string = CRS(prj))
#  shapeZona_sp_ <- as(shapeZona_sp, "Spatial")
#  normal = unlist(over(shapeZona_sp_, posTri, returnList = T)) # check which mesh triangles are inside the normal area
#  barrier.triangles = setdiff(1:tl, normal)
#  spdesismos = inla.barrier.pcmatern(meshSismos,
#                                     barrier.triangles = barrier.triangles,
#                                     prior.range = c(range_simulated, 0.5),   
#                                     prior.sigma = c(scale_simulated, 0.5))
#  
#  n_spde = spdesismos$f$n
#}
