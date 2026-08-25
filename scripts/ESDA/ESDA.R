# =============================================================================
# SCRIPT: ESDA.R
# =============================================================================
# PROPÓSITO:
#   Análisis exploratorio de datos espaciales para los eventos sísmicos 
#   de Colombia 2005–2020
#
#
# AUTOR: Jason Mauricio Romero Ríos
# UNIVERSIDAD: Universidad Distrital Francisco José de Caldas
# TESIS: Maestría en Ciencias de la Información y Comunicaciones – Geomática
# =============================================================================

# =============================================================================
# SECCIÓN 1: CARGA DE PAQUETES
# =============================================================================



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
library(prettymapr)

# =============================================================================
# SECCIÓN 2: CONFIGURACIÓN DE RUTAS Y PARÁMETROS GLOBALES
# =============================================================================

setwd("/home/jasonromeroia/Documents/personal/Tesis_MCIC/lgcp_inla_spde_sismos/")
source("scripts/ESDA/utils.R")

path_file_seismic <- "Data/gdf_espacial_2005_2020.gpkg" 

# Covariables
files_rds = "covariables_rds"


logo <- read_png_rgba("imagenes_doc/logo_ud.png")
logo_grob <- rasterGrob(logo, interpolate = TRUE)

# Capa de Departamentos
dept_sp = st_read("Data/departamentos_col.gpkg")


shapeZona_sp <- readRDS(paste0(files_rds,"/shapeZona_sp"))

# Borde de Colombia
shapeZona_sp = st_simplify(shapeZona_sp, dTolerance = 5000, preserveTopology = T)



# =============================================================================
# SECCIÓN CARGA DEL CATÁLOGO SÍSMICO 
# =============================================================================


sismosSp=st_read(path_file_seismic)

sismosSp$X = st_coordinates(sismosSp)[,1]
sismosSp$Y = st_coordinates(sismosSp)[,2]



# Columnas de interés

sprintf(names(sismosSp))

columns_view = c("MAGNITUD_MW")
column_mag = "MAGNITUD_MW"
columns_depth = "PROFUNDIDAD_KM"

# Ruta de las imagenes de salida
path_image_results <- "imagenes_doc"
dir_out_st <- file.path(path_image_results, "ESDA")

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


# Validación Sistema de referencia 3116 
cat(st_crs(shapeZona_sp)$epsg)
cat(st_crs(sismos_clip)$epsg)
cat(st_crs(dept_sp)$epsg)



shapeZona_sp_3857  <- st_transform(shapeZona_sp, 3857)
sismos_3857     <- st_transform(sismos_clip, 3857)
bb <- st_bbox(shapeZona_sp_3857)




#####################################################################
###### Mapa de eventos sísmicos anual ###############################
#####################################################################
#sismosSp_2005_2012 = subset(sismos_3857, YEAR >= 2005 & YEAR <= 2012)
#sismosSp_2013_2020 = subset(sismos_3857, YEAR >= 2013 & YEAR <= 2020)
#p1 <- crear_mapa_cientifico_anual(sismosSp_2005_2012, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2005_2012.png"))
#p2 <- crear_mapa_cientifico_anual(sismosSp_2013_2020, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2013_2020.png"))
eventos_sismicos <- crear_mapa_cientifico_anual(sismos_3857, shapeZona_sp_3857, paste0(path_image_results, "/mapa_cientifico_anual2005_2020.png"))

save_fig(eventos_sismicos, "mapa_cientifico_anual2005_2020.png", 6.5, 6.5)

#####################################################################
###### Eventos sísmicos por anno      ###############################
#####################################################################

unique(sismos_3857$YEAR)
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


to_math_range <- function(x) {
  y <- gsub("≤", "\\leq", x, fixed = TRUE)
  paste0("$", y, "$")  # now < is valid as-is in math mode
}
conteo_year = conteo_year_porcentaje[,c("YEAR", "n")]

hist_x_anno <- ggplot(conteo_year, aes(x = YEAR, y = n)) +
  geom_segment(aes(xend = YEAR, y = 0, yend = n),
               color = "grey78", linewidth = 0.8) +
  geom_point(color = "#2C6E9B", size = 3.5) +
  geom_text(aes(label = comma(n)),
            angle = 90, hjust = -0.35, size = 3.1, color = "grey25") +
  scale_x_continuous(breaks = seq(min(conteo_year$YEAR),
                                  max(conteo_year$YEAR), by = 2)) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "",
    subtitle = "",
    x = NULL,
    y = "Número de registros",
    caption  = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "grey90"),
    plot.title.position = "plot",
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey35", margin = margin(b = 12)),
    plot.caption  = element_text(color = "grey55", size = 8),
    axis.title.y  = element_text(margin = margin(r = 8))
  )


save_fig(hist_x_anno, "histograma_x_anno.png", 7.5, 6.5)







box_plot_anno = ggplot(df_pts, aes(x = factor(YEAR), y = .data[[columns_view]])) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(x = "Año", y = "Magnitug Mw") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


save_fig(box_plot_anno, "box_plot_anno.png", 6.5, 4.5)

####### MAGNITUD DE EVENTOS SISMICOS POR BREAKS ###############
###############################################################

breaks_m <- c(-Inf, 1, 2, 3, 4, 5, 6, 7)
labels_m <- c(
  "m < 1",
  "1 ≤ m < 2",
  "2 ≤ m < 3",
  "3 ≤ m < 4",
  "4 ≤ m < 5",
  "5 ≤ m < 6",
  "6 ≤ m < 7"
)

tabla_mag <- df_pts %>%
  transmute(
    bin = cut(.data[[column_mag]], breaks = breaks_m, labels = labels_m, right = FALSE)
  ) %>%
  count(bin, name = "Conteo") %>%
  mutate(`%` = Conteo / sum(Conteo)) %>%
  mutate(bin = factor(bin, levels = labels_m)) %>%
  arrange(bin)

# Versión formateada para LaTeX (rangos en modo matemático + escape de %)
tabla_mag_fmt <- tabla_mag %>%
  mutate(
    Conteo_fmt = number(Conteo, big.mark = ".", decimal.mark = ",", accuracy = 1),
    pct_fmt    = latex_escape(percent(`%`, accuracy = 0.01, decimal.mark = ",")),
    Magnitud   = to_math_range(bin)
  ) %>%
  select(Magnitud, Conteo = Conteo_fmt, `%` = pct_fmt)

print(tabla_mag_fmt)          # tabla legible en consola
xtable(tabla_mag_fmt)     
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


save_fig(conteo_dpto_global, "conteo_dpto_global.png", 6.5, 6.5)


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


save_fig(conteo_depto_anno_10, "conteo_depto_anno_10.png", 8.5, 6.5)

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
  labs(x = "Año", y = "Profundidad (km)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


save_fig(box_plot_depth, "box_plot_depth_anno.png", 6.5, 6.5)

#####################################################################################
############### Fallas Geológicas ###################################################
#####################################################################################
falla_rumbo_sinestral = readRDS("Data/Covariables/falla_rumbo_sinestral.rds")
## dextral
falla_rumbo_dextral = readRDS("Data/Covariables/falla_rumbo_dextral.rds")
## normal
falla_normal = readRDS("Data/Covariables/falla_normal.rds")
## inversas
falla_inversa = readRDS("Data/Covariables/falla_inversa.rds")


pliegues = readRDS("Data/Covariables/pliegues.rds")

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
  "y Comunicaciones\n",
  "\n",
  "Elaborado por: Jason Romero\n\n",
  "Sistema de Referencia: EPSG:3857\n\n"
)





# Contar eventos totales para subtítulo
total_eventos <- nrow(sismos_3857)




#header_con_estructuras <- paste0(header_map_basic, "Estructuras")
titulo_leyenda <- "Estructuras geológicas"

# --- 2) Mapa con la leyenda (incluye el header como title) ---
p <- ggplot() +
  annotation_map_tile(
    type = "cartolight",   # estilo
    zoom = 5,              # nivel de detalle
    cachedir = tempdir()   # para cache local
  ) +
  geom_sf(data = capas_union, aes(color = tipo, linetype = tipo), linewidth = 0.5) +
  scale_color_manual(values = cols, name = titulo_leyenda) +
  scale_linetype_manual(values = lts, name = titulo_leyenda) +
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
  draw_grob(logo_grob, x = 0.82, y = 0.9, width = 0.18, height = 0.18,
            hjust = 0.5, vjust = 1)




save_fig(final_plot, "fallas.png", 8.5, 8.5)

#####################################################################################
############### Topografia        ###################################################
#####################################################################################





topography = readRDS("Data/Covariables/topography.rds")
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




save_fig(p_topo_final, "topografia.png", 14, 8)


#####################################################################################
############### Anomalia Isostasia        ###################################################
#####################################################################################


isostasia = readRDS("Data/Covariables/isostasia.rds")
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

save_fig(p_isostasia_final, "p_isostasia_final.png", 14, 8)

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

png(paste0(path_image_results,"/ESDA/mCuadCount.png",width = 1024, height = 1024, pointsize = 20))

plot(Cuad1,cex=1,main="",col="blue")
dev.off() ###


png(paste0(path_image_results,"/ESDA/mCuadCountResiduales.png",width = 1024, height = 800, pointsize = 20))
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



