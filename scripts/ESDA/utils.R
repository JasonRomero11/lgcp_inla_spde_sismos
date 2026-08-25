




crear_mapa_cientifico_anual <- function(sismos_3857, shapeZona_3857, path_output) {
  
  bb <- st_bbox(shapeZona_3857)
  
  # Filtrar solo años 2005-2020 y ordenar
  sismos_filtrado <- sismos_3857 %>%
    filter(YEAR >= 2005 & YEAR <= 2020) %>%
    mutate(YEAR = factor(YEAR, levels = 2005:2020))
  
  # Crear clasificación de magnitud
  sismos_filtrado <- sismos_filtrado %>%
    mutate(
      mag_clase = cut(
        .data[[column_mag]],
        breaks = c(0, 3, 4, 5, 6, Inf),
        labels = c("< 3.0", "3.0 - 4.0", "4.0 - 5.0", "5.0 - 6.0", "> 6.0"),
        include.lowest = TRUE
      )
    )
  
  # Contar eventos por año para subtítulo
  conteo_anual <- sismos_filtrado %>%
    st_drop_geometry() %>%
    group_by(YEAR) %>%
    summarise(n = n(), .groups = "drop")
  
  total_eventos <- sum(conteo_anual$n)
  
  # Crear el gráfico
  p_cientifico_anual <- ggplot() +
    # Contorno de zona de estudio
    geom_sf(
      data = shapeZona_3857,
      fill = NA,
      color = "gray30",
      linewidth = 0.3
    ) +
    
    # Sismos con clasificación por magnitud
    geom_sf(
      data = sismos_filtrado,
      aes(fill = mag_clase),
      shape = 21,
      color = "gray40",
      size = 1.2,
      stroke = 0.15,
      alpha = 0.75
    ) +
    
    # Paleta de colores secuencial
    scale_fill_brewer(
      palette = "YlOrRd",
      name = "magnitud (Mw)",
      drop = FALSE,
      na.translate = FALSE
    ) +
    
    # Facetas por año (4 columnas para 16 años)
    facet_wrap(
      ~ YEAR,
      ncol = 4,
      labeller = labeller(YEAR = function(x) paste0(x))
    ) +
    
    # Coordenadas
    coord_sf(
      crs = st_crs(3857),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    
    # Tema científico minimalista
    theme_bw(base_size = 9) +
    theme(
      # Título y subtítulo
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5,
        color = "gray40",
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = 7,
        hjust = 1,
        color = "gray50",
        margin = margin(t = 10)
      ),
      
      # Paneles de facetas
      strip.text = element_text(
        face = "bold",
        size = 9,
        color = "gray20",
        margin = margin(t = 3, b = 3)
      ),
      strip.background = element_rect(
        fill = "gray95",
        color = "gray70",
        linewidth = 0.3
      ),
      
      # Panel
      panel.grid.major = element_line(color = "gray92", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray60", linewidth = 0.3),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(0.3, "lines"),
      
      # Ejes
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      # Leyenda
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.background = element_rect(fill = "white", color = "gray60", linewidth = 0.3),
      legend.margin = margin(5, 10, 5, 10),
      legend.box.margin = margin(t = 5),
      
      # Márgenes generales
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    
    # Guías de leyenda
    guides(
      fill = guide_legend(
        title.position = "left",
        title.vjust = 0.8,
        nrow = 1,
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    
    labs(
      caption = paste0(
        "Fuente: Servicio Geológico Colombiano | ",
        "Elaboración: [Institución] | "
      )
    )
  
  # Guardar en PNG alta resolución
  ggsave(
    filename = path_output,
    plot = p_cientifico_anual,
    width = 12,
    height = 14,
    dpi = 300,
    bg = "white"
  )
  
  # Retornar el gráfico
  return(p_cientifico_anual)
}


# Guardar figuras
save_fig <- function(p, name, w, h)
  ggsave(file.path(dir_out_st, name), plot = p, width = w, height = h, dpi = 300)



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



reproj_rast <- function(r, target_epsg = 3116, method = "bilinear") {
  
  # --- 1. Acepta tanto Raster* como SpatRaster ----------------------------
  is_raster <- inherits(r, "Raster")
  if (is_raster) r <- terra::rast(r)
  
  # --- 2. ¿Comparte CRS con el objetivo? ----------------------------------
  tgt_crs <- paste0("EPSG:", target_epsg)
  needs_proj <- !terra::same.crs(r, tgt_crs)   # FALSE si ya coincide
  
  # --- 3. Reproyecta solo cuando es necesario -----------------------------
  if (needs_proj) {
    message("   reproyectando a ", tgt_crs, " …")
    r <- terra::project(r, tgt_crs, method = method)
  } else {
    message("   CRS ya es ", tgt_crs)
  }
  
  # --- 4. Devuelve en la misma clase que entró ----------------------------
  if (is_raster) r <- raster::brick(r)
  return(r)
}
