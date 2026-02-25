library(INLA)
library(sp)
library(ggplot2)
library(sf)
library(spatial)
library(spData)
library(spdep)
library(maps)
library(gridExtra)

library(deldir)
library(raster)
library(viridis)
library("deldir")
library(terra)
library(lubridate)
# Cargar el paquete
library(fmesher)
library(sqldf)
library(dplyr)
library(leaflet)
library(scales)
library(knitr)
library(purrr)
library(ggspatial)
if (!requireNamespace("fmesher", quietly = TRUE)) {
  install.packages("fmesher", repos = "https://inla.r-inla-download.org/R/stable") # Instalar desde el repositorio de INLA
}
library(xtable)
library(forcats)
library(grid)
library(cowplot)
library(png)

library(terra)
library(tidyterra)
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)
library(inlabru)
library(xtable)
setwd("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/")

path_image_results = 'imagenes_doc'

path_file_seismic = "data_new/EventosColPointsPlanas31162005_2020_continental.gpkg"
files_rds = "covariables_rds"#"data_rds" #covariables_rds"
#sismosSp=st_read("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/eventos_sismicos1995_2024/eventos_declustering_2000_2020_3116.geojson")

sismosSp=st_read(path_file_seismic)


dept_sp = st_read("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/data_new/departamentos_col.geojson")


shapeZona_sp <- readRDS(paste0(files_rds,"/shapeZona_sp"))


shapeZona_sp = st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = T) #20000
names(sismosSp)

#sismosSp$fecha <- trimws(sismosSp$Date)
#sismosSp$fecha_convertida <- dmy_hms(sismosSp$fecha)
#sismosSp$fecha_convertida <- as.Date(sismosSp$fecha_convertida)
#sismosSp$YEAR <- format(sismosSp$fecha_convertida, "%Y")
sismosSp$X = st_coordinates(sismosSp)[,1]
sismosSp$Y = st_coordinates(sismosSp)[,2]




columns_view = c("mag")
column_mag = "mag"
columns_depth = "depth"



crs_zona <- st_crs(shapeZona_sp)  #EPSG:3116
if (is.na(crs_zona$epsg) || crs_zona$epsg != 3116) {
  stop("shapeZona_sp no está en EPSG:3116. Reprojéctalo a 3116 antes de continuar.")
}

# Transformar sismos (si vinieran en 4326, se pasan a 3116)
if (st_crs(sismosSp)$epsg != 3116) {
  sismos3116 <- st_transform(sismosSp, 3116)
} else {
  sismos3116 <- sismosSp
}

# 3) Geometrías válidas y recorte a la zona de estudio
shapeZona_sp   <- st_make_valid(shapeZona_sp)
sismos3116  <- st_make_valid(sismos3116)



# Recorte eficiente: quedarse solo con puntos dentro de shapeZona_sp
# (st_filter suele ser más rápido que st_intersection para puntos vs polígonos)
sismos_clip <- st_filter(sismos3116, st_union(shapeZona_sp))



unique(sismos_clip$YEAR)


cat(st_crs(shapeZona_sp)$epsg)
cat(st_crs(sismos_clip)$epsg)
cat(st_crs(dept_sp)$epsg)



shapeZona_sp_3857  <- st_transform(shapeZona_sp, 3857)
sismos_3857     <- st_transform(sismos_clip, 3857)
bb <- st_bbox(shapeZona_sp_3857)


names(sismos_3857)


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
        mag,
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
      name = "magnitud (ML)",
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


#####################################################################
###### Mapa de eventos sísmicos anual ###############################
#####################################################################
sismosSp_2005_2012 = subset(sismos_3857, YEAR >= 2005 & YEAR <= 2012)
sismosSp_2013_2020 = subset(sismos_3857, YEAR >= 2013 & YEAR <= 2020)
p1 <- crear_mapa_cientifico_anual(sismosSp_2005_2012, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2005_2012.png"))
p2 <- crear_mapa_cientifico_anual(sismosSp_2013_2020, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2013_2020.png"))
p3 <- crear_mapa_cientifico_anual(sismos_3857, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2005_2020.png"))

#####################################################################
###### Eventos sísmicos por anno      ###############################
#####################################################################

coords <- st_coordinates(sismos_3857)
df_pts <- sismos_3857 |>
  st_drop_geometry() |>
  mutate(X = coords[,1], Y = coords[,2]) |>
  filter(!is.na(YEAR))


conteo_year_porcentaje <- df_pts %>%
  count(YEAR) %>%
  mutate(
    porcentaje = scales::percent(n / sum(n), accuracy = 0.01, decimal.mark = ",")
  )


conteo_year_porcentaje <- df_pts %>%
  count(YEAR) %>%
  mutate(
    porcentaje = scales::percent(n / sum(n), accuracy = 0.01, decimal.mark = ",")
  )

latex_escape <- function(x) {
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x
}

tabla_conteo_year_porcentaje <- conteo_year_porcentaje %>%
  mutate(
    pct_fmt = latex_escape((porcentaje)),
  )


print(
  xtable(conteo_year_porcentaje),
  include.rownames = FALSE,
  sanitize.text.function = identity,         # keep our LaTeX
  sanitize.colnames.function = identity,     # keep \% in header
  comment = FALSE                            # remove timestamp comment line
)

conteo_year = conteo_year_porcentaje[,c("YEAR", "n")]
hist_x_anno = ggplot(conteo_year, aes(x = YEAR, y = n)) +
  geom_segment(aes(x = YEAR, xend = YEAR, y = 0, yend = n),
               color = "grey60", linewidth = 1) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(aes(label = n), vjust = -0.9, angle = 90, size = 3.5) +
  scale_x_continuous(breaks = seq(min(conteo_year$YEAR),
                                  max(conteo_year$YEAR), 1)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Año",
    y = "Número de registros"
  )


ggsave(paste0(path_image_results,"/histograma_x_anno.png"), 
       plot = hist_x_anno, 
       width = 12, height = 8, dpi = 300)


box_plot_anno = ggplot(df_pts, aes(x = factor(YEAR), y = .data[[column_mag]])) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(x = "Año", y = column_mag) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(path_image_results,"/box_plot_anno.png"), 
       plot = box_plot_anno, 
       width = 12, height = 8, dpi = 300)


########################################################
### NÚMERO DE EVENTOS POR DEPARTAMENTO #################
########################################################

conteo_dpto <- sismos_clip %>%
  st_join(dept_sp) %>%    
  st_drop_geometry() %>%  
  group_by(ADM1_ES) %>%   
  summarise(n_sismos = n(), .groups = "drop") %>%
  filter(!is.na(ADM1_ES))  # <- Agregar esta línea

conteo_dpto_global = ggplot(conteo_dpto,
                            aes(x = fct_reorder(ADM1_ES, n_sismos), y = n_sismos)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Departamento",
    y = "Número de sismos"
  ) +
  theme_minimal()

ggsave(paste0(path_image_results,"/conteo_dpto_global.png"), 
       plot = conteo_dpto_global, 
       width = 12, height = 8, dpi = 300)



conteo_dpto <- sismos_clip %>%
  st_join(dept_sp) %>%             
  st_drop_geometry() %>%            
  group_by(ADM1_ES, YEAR) %>%        
  summarise(n_sismos = n(), .groups = "drop")

totales <- conteo_dpto %>%
  group_by(ADM1_ES) %>%
  summarise(total_sismos = sum(n_sismos), .groups = "drop")


top10    <- totales %>% arrange(desc(total_sismos)) %>% slice_head(n = 10) %>% pull(ADM1_ES)


conteo_top10    <- conteo_dpto %>% filter(ADM1_ES %in% top10)



conteo_depto_anno_10 = ggplot(conteo_top10,
                              aes(x = YEAR, y = n_sismos, colour = ADM1_ES)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Año",
    y = "Número de sismos",
    colour = "Departamento"
  ) +
  theme_minimal()

ggsave(paste0(path_image_results,"/conteo_depto_anno_10.png"), 
       plot = conteo_depto_anno_10, 
       width = 12, height = 8, dpi = 300)



###############################################################
####### PROFUNDIDAD DE EVENTOS SISMICOS POR BREAKS ############
###############################################################

## Resumen estadístico por profundiad de eventos sísmicos

breaks_p <- c(-Inf, 100, 200, 300, 400, 500, Inf)
labels_p <- c(
  "p < 100",
  "100 ≤ p < 200",
  "200 ≤ p < 300",
  "300 ≤ p < 400",
  "400 ≤ p < 500",
  "500 ≤ p"
)


to_math_range <- function(x) {
  y <- gsub("≤", "\\leq", x, fixed = TRUE)
  paste0("$", y, "$")  # now < is valid as-is in math mode
}
names(df_pts)
tabla_depth <- df_pts %>%
  transmute(
    bin = cut(.data[[columns_depth]], breaks = breaks_p, labels = labels_p, right = FALSE)
  ) %>%
  count(bin, name = "Conteo") %>%
  mutate(`%` = Conteo / sum(Conteo)) %>%
  mutate(bin = factor(bin, levels = labels_p)) %>%
  arrange(bin)

tabla_depth_fmt <- tabla_depth %>%
  mutate(
    Conteo_fmt = number(Conteo, big.mark = ".", decimal.mark = ",", accuracy = 1),
    `%_fmt`    = percent(`%`, accuracy = 0.01, decimal.mark = ",")
  ) %>%
  select(Profundidad = bin, Conteo = Conteo_fmt, `%` = `%_fmt`)

tabla_depth_fmt <- tabla_depth %>%
  mutate(
    Conteo_fmt   = number(Conteo, big.mark = ".", decimal.mark = ",", accuracy = 1),
    pct_fmt      = latex_escape(percent(`%`, accuracy = 0.01, decimal.mark = ",")),
    Profundidad = to_math_range(bin)
  ) %>%
  select(Magnitud = Profundidad, Conteo = Conteo_fmt, `%` = pct_fmt)

xtable(tabla_depth_fmt)



box_plot_depth = ggplot(df_pts, aes(x = factor(YEAR), y = .data[[columns_depth]])) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(x = "Año", y = columns_depth) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(path_image_results,"/box_plot_depth_anno.png"), 
       plot = box_plot_depth, 
       width = 12, height = 8, dpi = 300)

#####################################################################################
############### Fallas Geológicas ###################################################
#####################################################################################
falla_rumbo_sinestral = st_read("data/shapefile/fallas/falla_rumbo_sinestral.shp")
## dextral
falla_rumbo_dextral = st_read("data/shapefile/fallas/falla_rumbo_dextral.shp")
## normal
falla_normal = st_read("data/shapefile/fallas/falla_normal.shp")
## inversas
falla_inversa = st_read("data/shapefile/fallas/falla_inversa.shp")

#lineamientos = st_read("data/shapefile/fallas/lineamiento.shp")

pliegues = st_read("data/shapefile/pliegues/pliegues.shp")

crs_target  = 3857
layers <- list(
  "Falla rumbo sinistral" = falla_rumbo_sinestral,
  "Falla rumbo dextral"   = falla_rumbo_dextral,
  "Falla normal"          = falla_normal,
  "Falla inversa"         = falla_inversa,
  #"Lineamientos"          = lineamientos,
  "Pliegues"              = pliegues
) |> 
  lapply(\(x) if (is.na(st_crs(x))) x else st_transform(x, crs_target))




# === 3) Unir en un solo objeto con tipo de estructura ===
capas_union <- bind_rows(
  lapply(names(layers), \(nm) mutate(layers[[nm]], tipo = nm))
)

capas_union_sum <- capas_union %>%
  mutate(long_km = as.numeric(st_length(geometry)) / 1000)

# --- resumen por tipo ---
resumen_longitudes <- capas_union_sum %>%
  group_by(tipo) %>%
  summarise(
    n_lineas   = n(),
    total_km   = sum(long_km, na.rm = TRUE),
    promedio_km = mean(long_km, na.rm = TRUE),
    max_km     = max(long_km, na.rm = TRUE),
    min_km     = min(long_km, na.rm = TRUE),
    .groups = "drop"
  )

print(resumen_longitudes)
capas_union <- st_cast(capas_union, "MULTILINESTRING", warn = FALSE)

# === 4) Estilos (colores y linetipos consistentes) ===
cols <- c(
  "Falla rumbo sinistral" = "#c51b7d",
  "Falla rumbo dextral"   = "#2c7fb8",
  "Falla normal"          = "#41ab5d",
  "Falla inversa"         = "#fec44f",
  #"Lineamientos"          = "#9e9ac8",
  "Pliegues"              = "#fb6a4a"
)

lts <- c(
  "Falla rumbo sinistral" = "dotdash",
  "Falla rumbo dextral"   = "longdash",
  "Falla normal"          = "solid",
  "Falla inversa"         = "twodash",
  #"Lineamientos"          = "dotted",
  "Pliegues"              = "dashed"
)


header_map_basic <- paste0(
  #"Fecha: ", format(Sys.Date(), "%Y-%m-%d"), "\n",
  "Universidad Distrital Francisco José de Caldas \n",
  "\n",
  "Maestría en Ciencias de la Información \n",
  "y Comunicaciones Énfasis en Geomática\n",
  "\n",
  "Elaborado por: Jason Romero\n\n",
  "Sistema de Referencia: EPSG:3857\n\n"
)

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
# --- 3) Cargar la imagen y crear el grob (sin magick) ---
logo <- read_png_rgba("imagenes_doc/logo_ud.png")
logo_grob <- rasterGrob(logo, interpolate = TRUE)

# Contar eventos totales para subtítulo
total_eventos <- nrow(sismos_3857)

# Mapa con estilo científico similar a crear_mapa_cientifico_anual
p_dist <- ggplot() +
  # Contorno de zona de estudio
  geom_sf(
    data = shapeZona_sp_3857,
    fill = NA,
    color = "gray30",
    linewidth = 0.3
  ) +

  # Sismos (puntos simples, sin clasificación)
  geom_sf(
    data = sismos_3857,
    color = "#3B528B",
    size = 0.5,
    alpha = 0.6
  ) +

  # Coordenadas
  coord_sf(
    crs = st_crs(3857),
    expand = FALSE
  ) +

  # Tema científico minimalista (igual que crear_mapa_cientifico_anual)
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

    # Panel
    panel.grid.major = element_line(color = "gray92", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray60", linewidth = 0.3),
    panel.background = element_rect(fill = "white"),

    # Ejes
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),

    # Márgenes generales
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA)
  ) +

  labs(
    title = "Distribución de Eventos Sísmicos en Colombia (2005-2020)",
    subtitle = paste0("Total de eventos: ", format(total_eventos, big.mark = ".", decimal.mark = ",")),
    caption = paste0(
      "Fuente: Servicio Geológico Colombiano | ",
      "Elaboración: Jason Romero"
    )
  )

# --- Composición final con logo y header (similar a p_topo_final y p_isostasia_final) ---
p_dist_final <- ggdraw() +
  draw_plot(p_dist) +
  draw_grob(logo_grob, x = 0.82, y = 0.90, width = 0.18, height = 0.18,
            hjust = 0.5, vjust = 1) +
  draw_text(header_map_basic,
            x = 0.82, y = 0.70,
            hjust = 0.5, vjust = 1,
            size = 9, lineheight = 1)

# Mostrar
print(p_dist_final)

ggsave(paste0(path_image_results,"/mapa_sismos.png"), 
       plot = p_dist_final, 
       width = 12, height = 8, dpi = 300)

header_con_estructuras <- paste0(header_map_basic, "Estructuras")
# --- 2) Mapa con la leyenda (incluye el header como title) ---
p <- ggplot() +
  annotation_map_tile(
    type = "cartolight",   # estilo
    zoom = 5,              # nivel de detalle
    cachedir = tempdir()   # para cache local
  ) +
  geom_sf(data = capas_union, aes(color = tipo, linetype = tipo), linewidth = 0.5) +
  scale_color_manual(values = cols, name = header_con_estructuras) +
  scale_linetype_manual(values = lts, name = header_con_estructuras) +
  guides(
    color    = guide_legend(override.aes = list(linewidth = 1), byrow = TRUE),
    linetype = guide_legend(byrow = TRUE)
  ) +
  annotation_scale(location = "bl", text_cex = 0.8, line_width = 0.6, height = unit(0.2, "cm")) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"), width = unit(1.1, "cm")
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, hjust = 0.5),   # centra el texto
    legend.title.align = 0.5,                             # centra todas las líneas del header
    legend.title.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey85")
  )

# --- 4) Componer: colocar el logo por encima de la zona de la leyenda ---
# Ajusta x/y/width/height para posicionarlo exactamente donde lo quieres
final_plot <- ggdraw(p) +
  draw_grob(logo_grob, x = 0.75, y = 0.9, width = 0.18, height = 0.18,
            hjust = 0.5, vjust = 1)


ggsave(paste0(path_image_results,"/fallas.png"), 
       plot = final_plot, 
       width = 12, height = 8, dpi = 300)



#####################################################################################
############### Topografia        ###################################################
#####################################################################################


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


topography <- terra::rast("data/raster/Topography_Colombia_2000m.tif")
topography = reproj_rast(topography, 3857)
topography[topography == 0] <- NA


p_topo <- ggplot() +
  annotation_map_tile(
    type = "cartolight",
    zoom = 5,
    cachedir = tempdir()
  ) +
  geom_spatraster(data = topography, alpha = 0.85) +  # transparencia para ver el basemap
  scale_fill_viridis_c(
    name = "Elevación (m)",
    option = "C",   # plasma
    direction = 1,
    na.value = NA
  ) +
  coord_sf(crs = st_crs(3857), expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    legend.position  = "right",
    panel.grid.major = element_line(color = "grey80", size = 0.2)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering
  )

# --- Extraer leyenda vertical (desde el mismo p_topo) ---
legend_topo <- get_legend(
  p_topo + theme(legend.position = "right")
)

# --- Composición final: mapa sin leyenda + logo + header + leyenda vertical ---
p_topo_final <- ggdraw() +
  draw_plot(p_topo + theme(legend.position = "none")) +
  draw_grob(logo_grob, x = 0.82, y = 0.90, width = 0.18, height = 0.18,
            hjust = 0.5, vjust = 1) +
  draw_text(header_map_basic,
            x = 0.82, y = 0.70,
            hjust = 0.5, vjust = 1,
            size = 9, lineheight = 1) +
  # Ajusta x/y/width/height para ubicar la leyenda donde te convenga
  draw_plot(legend_topo, x = 0.74, y = 0.27, width = 0.18, height = 0.28)

p_topo_final
ggsave(paste0(path_image_results,"/topografia.png"), 
       plot = p_topo_final, 
       width = 14, height = 8, dpi = 300)





#####################################################################################
############### Anomalia Isostasia        ###################################################
#####################################################################################



isostasia <- terra::rast("data/raster/mosaico_isostasia_cor.tif")
isostasia = reproj_rast(isostasia, 3857)
summary(isostasia$mosaico_isostasia_cor)

#vals <- values(isostasia, na.rm = TRUE)
#qs <- quantile(vals, probs = c(0.1, 0.9), na.rm = TRUE)

#isostasia <- clamp(isostasia, lower = qs[1], upper = qs[2], values = TRUE)


rminmax <- global(isostasia, fun = range, na.rm = TRUE)[1,]
vmin <- rminmax[1]; vmax <- rminmax[2]
brks_pretty <- pretty(c(vmin$X1, vmax$X2), n = 7)

#??pretty
p_isostasia <- ggplot() +
  annotation_map_tile(
    type = "cartolight",
    zoom = 5,
    cachedir = tempdir()
  ) +
  geom_spatraster(data = isostasia, alpha = 0.85) +  # transparencia para ver el basemap
  scale_fill_viridis_c(
    name = "Anomalías Isostáticas (mGal)",
    #limits = c(vmin, vmax),
    breaks = brks_pretty,
    labels = label_number(accuracy = 0.1),
    option = "C",
    direction = 1,
    na.value = NA
  ) +
  coord_sf(crs = st_crs(3857), expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    legend.position  = "right",
    panel.grid.major = element_line(color = "grey80", size = 0.2)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering
  )

# --- Extraer leyenda vertical (desde el mismo p_topo) ---
legend_isostasia <- get_legend(
  p_isostasia + theme(legend.position = "right")
)

# --- Composición final: mapa sin leyenda + logo + header + leyenda vertical ---
p_isostasia_final <- ggdraw() +
  draw_plot(p_isostasia + theme(legend.position = "none")) +
  draw_grob(logo_grob, x = 0.82, y = 0.90, width = 0.18, height = 0.18,
            hjust = 0.5, vjust = 1) +
  draw_text(header_map_basic,
            x = 0.82, y = 0.70,
            hjust = 0.5, vjust = 1,
            size = 9, lineheight = 1) +
  # Ajusta x/y/width/height para ubicar la leyenda donde te convenga
  draw_plot(legend_isostasia, x = 0.74, y = 0.27, width = 0.18, height = 0.28)

p_isostasia_final
ggsave(paste0(path_image_results,"/p_isostasia_final.png"), 
       plot = p_isostasia_final, 
       width = 14, height = 8, dpi = 300)




#####################################################################################
############### Análisis de Aleatoriedad Espacial Completa – CSR      ###################################################
#####################################################################################
SP <- as.owin(shapeZona_sp)
pSisAux <- as.ppp(sismos_clip)
pSismos <- ppp(pSisAux$x,pSisAux$y,window=as(SP,"owin"))
unitname(pSismos)=c("meter","meters")
plot(pSismos$window)
Cuad1 <- quadratcount(pSismos,ny=6,nx=6)
Cuad2 <- quadrat.test(pSismos,ny=6,nx=6)

png(paste0(path_image_results,"/mCuadCount.png",width = 1024, height = 1024, pointsize = 20))
plot(Cuad1,cex=1,main="",col="blue")
dev.off() ###


png(paste0(path_image_results,"/mCuadCountResiduales.png",width = 1024, height = 800, pointsize = 20))
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Gráfico 1: Observados vs Esperados
plot(Cuad2$observed, Cuad2$expected, 
     xlab = "Valores observados", 
     ylab = "Valores esperados",
     main = "",
     pch = 19, col = "blue", cex = 1.2)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2) # Línea y=x
grid()

# Gráfico 2: Observados vs Residuales
plot(Cuad2$observed, Cuad2$residuals, 
     xlab = "Valores observados", 
     ylab = "Valores residuales",
     main = "",
     pch = 19, col = "darkgreen", cex = 1.2)
abline(h = 0, col = "red", lwd = 2, lty = 2) # Línea en y=0
grid()

par(op)
dev.off() ###



################################################################################
####### 

#####################################################################################
############### Mapa de Sismos con Zoom a Zona de Interés ###########################
############### Sin API Key - Tiles Gratuitos #######################################
#####################################################################################

# Librerías necesarias
library(sf)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(png)
library(grid)

# === 1) Cargar zona de interés ===
zona_interes <- st_read("/home/jasonromeroia/Documents/Personal/TesisUDFJCMCIC/solucion2025/earthquakes_lgcp_inla/data_new/zona_interes.gpkg")

# === 2) Transformar a EPSG:3857 (Web Mercator) ===
crs_target <- 3857
zona_interes_3857 <- st_transform(zona_interes, crs_target)

# Asegurarse que sismos también esté en 3857
# (Asumiendo que ya tienes sismos_3857 cargado de tu código anterior)

# === 3) Obtener bbox de la zona de interés para el zoom ===
bbox_zona <- st_bbox(zona_interes_3857)

# Añadir un pequeño buffer al bbox (5% de margen)
margin_x <- (bbox_zona["xmax"] - bbox_zona["xmin"]) * 0.05
margin_y <- (bbox_zona["ymax"] - bbox_zona["ymin"]) * 0.05

xlim <- c(bbox_zona["xmin"] - margin_x, bbox_zona["xmax"] + margin_x)
ylim <- c(bbox_zona["ymin"] - margin_y, bbox_zona["ymax"] + margin_y)

# === 4) Header del mapa ===
header_map_basic <- paste0(
  "Universidad Distrital Francisco José de Caldas \n",
  "\n",
  "Maestría en Ciencias de la Información \n",
  "y Comunicaciones Énfasis en Geomática\n",
  "\n",
  "Elaborado por: Jason Romero\n\n",
  "Sistema de Referencia: EPSG:3857\n\n"
)

# === 5) Función para leer PNG ===
read_png_rgba <- function(path) {
  x <- png::readPNG(path)
  d <- dim(x)
  if (length(d) == 2) {
    x <- array(rep(x, 3), dim = c(d[1], d[2], 3))
  } else if (d[3] == 2) {
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

# Cargar logo
logo <- read_png_rgba("imagenes_doc/logo_ud.png")
logo_grob <- rasterGrob(logo, interpolate = TRUE)



p_zona <- ggplot() +
  annotation_map_tile(
    type = "osm",           
    zoom = 14,               
    cachedir = tempdir()
  ) +
  geom_sf(
    data = zona_interes_3857, 
    fill = NA, 
    color = "red", 
    linewidth = 1.5,
    linetype = "solid"
  ) +
  geom_sf(
    data = sismos_3857, 
    color = "#3B528B", 
    size = 1.5, 
    alpha = 0.01          # Alta transparencia
  ) +
  # Zoom a la zona de interés
  coord_sf(
    xlim = xlim, 
    ylim = ylim, 
    expand = FALSE
  ) +
  # Escala
  annotation_scale(
    location = "bl", 
    text_cex = 0.8, 
    line_width = 0.6, 
    height = unit(0.2, "cm")
  ) +
  # Norte
  annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"), 
    width = unit(1.1, "cm")
  ) +
  # Tema
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey85")
  )


p_zona_final <- ggdraw() +
  draw_plot(p_zona) +
  draw_grob(
    logo_grob, 
    x = 0.83, y = 0.90, 
    width = 0.18, height = 0.18,
    hjust = 0.5, vjust = 1
  ) +
  draw_text(
    header_map_basic, 
    x = 0.83, y = 0.7, 
    hjust = 0.5, vjust = 1, 
    size = 9
  )

# Mostrar
print(p_zona_final)

# === 8) Guardar ===
ggsave(
  paste0(path_image_results, "/mapa_sismos_zona_interes.png"), 
  plot = p_zona_final, 
  width = 12, height = 8, dpi = 300
)

