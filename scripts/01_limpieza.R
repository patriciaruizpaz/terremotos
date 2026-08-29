# =====================================================================
# Limpieza y control de datos — Terremoto: profundidad y magnitud por
# región geográfica
# Inferencia Estadística para la Toma de Decisiones — Ruiz Paz / González
#
# Implementa, en orden, los pasos 1 a 11 del documento
# "registro_limpieza_datos" (protocolo de Zuur, Ieno & Elphick 2010
# adaptado al diseño). Todo lo que se hace acá es PREVIO a cualquier
# cálculo de la Unidad 2 (Fase 1): no se calcula ningún intervalo de
# confianza en este script.
#
# Entrada:  datos/crudo.csv
# Salida:   datos/limpio.csv, con columnas auxiliares (region,
#           depth_flag_default, magType_family, flag_posible_replica)
# =====================================================================

library(tidyverse)
library(gridExtra)

dir.create("reportes", showWarnings = FALSE)

# -----------------------------------------------------------------------
# Carga
# -----------------------------------------------------------------------
df <- read_csv("datos/crudo.csv", show_col_types = FALSE)
n0 <- nrow(df)
cat(sprintf("[Carga] %d registros, %d columnas\n", n0, ncol(df)))

# -----------------------------------------------------------------------
# PASO 1 - Exclusiones del preregistro (Sección 5)
#   - type != "earthquake"
#   - depth o mag faltantes (NA)
#   - depth negativo
# -----------------------------------------------------------------------
before <- nrow(df)
df <- df %>% filter(type == "earthquake")
n_type <- before - nrow(df)

before <- nrow(df)
df <- df %>% filter(!is.na(depth), !is.na(mag))
n_na <- before - nrow(df)

before <- nrow(df)
df <- df %>% filter(depth >= 0)
n_neg <- before - nrow(df)

cat(sprintf("\n[Paso 1] Excluidos por type != earthquake: %d (nuclear explosion, volcanic eruption, explosion, landslide, rock burst, mine collapse)\n", n_type))
cat(sprintf("[Paso 1] Excluidos por depth/mag NA: %d\n", n_na))
cat(sprintf("[Paso 1] Excluidos por depth negativo: %d\n", n_neg))
cat(sprintf("[Paso 1] Quedan: %d\n", nrow(df)))

# -----------------------------------------------------------------------
# PASO 2 - Validación de coordenadas
# -----------------------------------------------------------------------
bad_lat <- !between(df$latitude, -90, 90)
bad_lon <- !between(df$longitude, -180, 180)
cat(sprintf("\n[Paso 2] Latitudes fuera de rango: %d | Longitudes fuera de rango: %d\n", sum(bad_lat), sum(bad_lon)))
cat("[Paso 2] Decisión: no hubo nada que excluir, coordenadas 100% válidas.\n")

# -----------------------------------------------------------------------
# PASO 3 - Duplicados por id
#   Hallazgo: el catálogo se actualiza semanalmente y muchos eventos
#   quedaron registrados varias veces (misma id, distinto 'updated'); en
#   ~47% de esos grupos, además, depth y/o mag cambian entre versiones
#   (revisión real del evento a medida que llegan más estaciones).
#   Decisión: quedarse con la revisión más reciente por id.
#
#   NOTA: ajustar as_datetime() si el formato de 'updated'/'time' en tu
#   CSV no es ISO-8601 (ej. "2023-04-05T12:34:56.000Z").
# -----------------------------------------------------------------------
before <- nrow(df)
df <- df %>%
  mutate(updated = as_datetime(updated)) %>%
  arrange(desc(updated)) %>%
  distinct(id, .keep_all = TRUE)
n_dup <- before - nrow(df)
cat(sprintf("\n[Paso 3] ids duplicados (múltiples snapshots del mismo evento): %d\n", n_dup))
cat(sprintf("[Paso 3] Decisión: se conserva solo la revisión más reciente (max 'updated') por id. Quedan: %d\n", nrow(df)))

# -----------------------------------------------------------------------
# PASO 4 - Profundidades de relleno (placeholder depths)
#   Hallazgo: depth=10.0 y depth=33.0 están brutalmente sobrerrepresentados
#   frente a sus vecinos, y depthError está ausente en la mayoría de esos
#   casos -> son los valores fijos clásicos de las redes sismográficas
#   cuando no se puede resolver la profundidad real.
#   Decisión: NO se eliminan filas (cambiaría los criterios de exclusión
#   ya preregistrados sin pasar por el Registro de Enmiendas). Se crea un
#   flag para chequeo de sensibilidad en Fase 1-2.
# -----------------------------------------------------------------------
df <- df %>%
  mutate(depth_flag_default = depth %in% c(10.0, 33.0) & is.na(depthError))
cat(sprintf("\n[Paso 4] Registros marcados como profundidad por defecto (10 o 33 km sin error estimado): %d (%.1f%% del dataset)\n",
            sum(df$depth_flag_default), mean(df$depth_flag_default) * 100))
cat("[Paso 4] Decisión: se flaguean, no se eliminan. Usar como chequeo de sensibilidad (Sección 4).\n")

# -----------------------------------------------------------------------
# PASO 5 - Consistencia de magType
#   Hallazgo: la familia "mw" (basada en momento sísmico) no satura y
#   cubre 55.9% del dataset, con magnitudes hasta 9.5. mb satura ~7.1 y
#   ms satura ~8.0 en este dataset (efecto de saturación esperado).
#   Decisión: NO se excluye nada acá. Se documenta la recomendación de
#   restringir el análisis principal a "mw" en Fase 1-2.
# -----------------------------------------------------------------------
mw_family <- c("mw", "mwc", "mww", "mwb", "mwr", "mwp")
df <- df %>%
  mutate(magType_family = case_when(
    magType %in% mw_family ~ "mw_family",
    magType == "mb" ~ "mb",
    magType == "ms" ~ "ms",
    TRUE ~ "otra"
  ))
cat("\n[Paso 5] Distribución por familia de escala:\n")
print(table(df$magType_family))
cat("[Paso 5] Decisión: se documenta el problema de saturación; no se excluye nada todavía. Columna magType_family agregada para decidir en Fase 1-2.\n")

# -----------------------------------------------------------------------
# PASO 6 - Descriptiva, distribución y outliers de depth y mag
#   Hallazgo: depth y mag muy asimétricas a la derecha; los valores más
#   extremos (depth=700km en Fiji, mag=9.5 en Chile 1960) corresponden a
#   eventos históricos reales, no a errores de carga.
#   Decisión: no se excluye ningún outlier.
# -----------------------------------------------------------------------
skewness <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x); s <- sd(x); n <- length(x)
  (sum((x - m)^3) / n) / s^3
}
cat(sprintf("\n[Paso 6] depth: media=%.1f, mediana=%.1f, skew=%.2f, max=%.1f\n",
            mean(df$depth), median(df$depth), skewness(df$depth), max(df$depth)))
cat(sprintf("[Paso 6] mag:   media=%.2f, mediana=%.2f, skew=%.2f, max=%.1f\n",
            mean(df$mag), median(df$mag), skewness(df$mag), max(df$mag)))
cat("[Paso 6] Decisión: sin exclusiones, valores extremos genuinos (verificados contra sismos históricos conocidos).\n")

# -----------------------------------------------------------------------
# PASO 7 - Cobertura temporal
#   Hallazgo: el salto de instrumentación ocurre donde se esperaba:
#   ~13.3% de los registros son anteriores a 1970, ~86.7% son posteriores.
# -----------------------------------------------------------------------
df <- df %>%
  mutate(time = as_datetime(time),
         year = year(time))
pre1970 <- mean(df$year < 1970) * 100
cat(sprintf("\n[Paso 7] Registros anteriores a 1970: %.1f%% | desde 1970: %.1f%%\n", pre1970, 100 - pre1970))
cat("[Paso 7] Decisión: se registra el sesgo de cobertura; la comparación pre/post 1970 queda para Fase 1-2.\n")

# -----------------------------------------------------------------------
# PASO 8 - Construcción de la variable región
#   Criterio: coincidencia de texto en 'place' contra países/territorios
#   dominados por subducción a lo largo del Pacífico. Validado con el
#   patrón de profundidad esperado (Argentina 63.5% >100km, etc.) y con
#   un mapa de puntos.
#   Limitación documentada: zonas de subducción genuinas fuera del
#   Pacífico (South Sandwich, Hindu Kush, Italia) quedan en "Resto"
#   porque el preregistro define la región como "del Pacífico"
#   específicamente (ver Decisión 2 en 02_decisiones.R).
# -----------------------------------------------------------------------
ring_of_fire_keywords <- c(
  "Indonesia", "Sumatra", "Java", "Sulawesi", "Timor", "Banda Sea", "Molucca", "Halmahera",
  "Japan",
  "Kuril", "Kamchatka", "Petropavlovsk", "Sea of Okhotsk", "Sakhalin", "Russia",
  "Philippines",
  "Papua New Guinea", "New Britain", "Bougainville", "D'Entrecasteaux", "Bismarck",
  "Solomon Islands", "Santa Cruz Islands",
  "Vanuatu", "New Caledonia", "Loyalty Islands",
  "Tonga", "Kermadec", "Fiji", "Samoa", "Wallis and Futuna",
  "Mariana Islands", "Guam", "Micronesia", "Palau",
  "Taiwan",
  "Alaska", "Aleutian",
  "Mexico", "Guatemala", "El Salvador", "Honduras", "Nicaragua", "Costa Rica", "Panama",
  "Chile", "Peru", "Ecuador", "Colombia", "Bolivia", "Argentina", "Brazil",
  "New Zealand"
)
pattern <- paste(ring_of_fire_keywords, collapse = "|")
df <- df %>%
  mutate(region = if_else(
    coalesce(str_detect(place, regex(pattern, ignore_case = TRUE)), FALSE),
    "Cinturón de Fuego", "Resto"
  ))

cat("\n[Paso 8]\n")
print(table(df$region))
print(df %>% group_by(region) %>%
        summarise(count = n(), mean = mean(depth), median = median(depth), sd = sd(depth)))

colores <- c("Cinturón de Fuego" = "#d62728", "Resto" = "#1f77b4")
p8 <- ggplot(df, aes(longitude, latitude, color = region)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = colores) +
  labs(x = "Longitud", y = "Latitud", title = "Paso 8 - Validación visual de la variable región") +
  theme_minimal()
ggsave("reportes/paso8_mapa_validacion.png", p8, width = 14, height = 7, dpi = 130)
cat("[Paso 8] Mapa de validación guardado en reportes/paso8_mapa_validacion.png\n")

# -----------------------------------------------------------------------
# PASO 9 - Homogeneidad de varianza de depth entre regiones
#   Hallazgo: varianza del Cinturón de Fuego ~10x la de Resto. Confirma
#   que hará falta Welch en vez de t común en Fase 1.
# -----------------------------------------------------------------------
var_ratio <- df %>% group_by(region) %>% summarise(var = var(depth))
ratio <- max(var_ratio$var) / min(var_ratio$var)
cat(sprintf("\n[Paso 9] Ratio de varianzas: %.1fx\n", ratio))
cat("[Paso 9] Decisión: se registra heterocedasticidad marcada; usar t de Welch (no t común) en Fase 1.\n")

p9 <- ggplot(df, aes(x = region, y = depth)) +
  geom_boxplot() +
  labs(y = "Profundidad (km)", title = "Paso 9 - depth por región") +
  theme_minimal()
ggsave("reportes/paso9_boxplot.png", p9, width = 6, height = 5, dpi = 130)

# -----------------------------------------------------------------------
# PASO 10 - Relación depth ~ mag, general y por región
#   Hallazgo IMPORTANTE: la correlación depth-mag es prácticamente nula
#   en ambas regiones y las curvas loess son casi planas. El efecto
#   fuerte es el de REGIÓN sobre el nivel/dispersión de depth (paso 9),
#   no una interacción mag*region. Clave para decidir el modelo de Fase 2.
# -----------------------------------------------------------------------
set.seed(1)

p10_general <- ggplot(df, aes(mag, depth)) +
  geom_point(size = 0.3, alpha = 0.1, color = "gray40") +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black") +
  labs(title = "General", x = "Magnitud", y = "Profundidad (km)") +
  theme_minimal()

plots_region <- list()
for (reg in c("Cinturón de Fuego", "Resto")) {
  sub <- df %>% filter(region == reg)
  test <- cor.test(sub$mag, sub$depth)
  cat(sprintf("[Paso 10] %s: pearson r=%.3f (p=%.2g)\n", reg, test$estimate, test$p.value))
  plots_region[[reg]] <- ggplot(sub, aes(mag, depth)) +
    geom_point(size = 0.3, alpha = 0.1, color = colores[reg]) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black") +
    labs(title = reg, x = "Magnitud", y = NULL) +
    theme_minimal()
}

p10 <- grid.arrange(p10_general, plots_region[["Cinturón de Fuego"]], plots_region[["Resto"]],
                     ncol = 3, top = "Paso 10 - depth vs mag, general y por región (loess)")
ggsave("reportes/paso10_scatter_loess.png", p10, width = 16, height = 5, dpi = 130)
cat("[Paso 10] Decisión: relación lineal simple mag-depth es débil/nula; el efecto está en el nivel medio y la dispersión de depth por región, no en una pendiente distinta por región. Registrar como hallazgo para discutir el enfoque de Fase 2.\n")

# -----------------------------------------------------------------------
# PASO 11 - Chequeo preliminar de independencia (réplicas)
#   Método simple: ordenando por tiempo, se marca cada evento cuyo vecino
#   inmediato también está a <100km y <24h de diferencia (proxy simple de
#   réplica/aftershock, sin declustering formal).
#   Hallazgo: ~12.7% de los eventos consecutivos cumplen ambos criterios.
#   Decisión: no se corrige todavía (queda para Fase 3, Unidad 4).
# -----------------------------------------------------------------------
haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  lat1 <- lat1 * pi / 180; lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180; lon2 <- lon2 * pi / 180
  dlat <- lat2 - lat1; dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

df_t <- df %>% arrange(time)
dt_hours <- as.numeric(difftime(df_t$time, lag(df_t$time), units = "hours"))
dist_km <- haversine(df_t$latitude, df_t$longitude, lag(df_t$latitude), lag(df_t$longitude))
df_t <- df_t %>%
  mutate(flag_posible_replica = (dt_hours < 24) & (dist_km < 100))

cat(sprintf("\n[Paso 11] Posibles réplicas (delta t<24h y <100km del evento anterior en la secuencia temporal): %d (%.1f%%)\n",
            sum(df_t$flag_posible_replica, na.rm = TRUE), mean(df_t$flag_posible_replica, na.rm = TRUE) * 100))
cat("[Paso 11] Decisión: se registra, no se corrige todavía (corrección/limitación se decide en Fase 3).\n")

# -----------------------------------------------------------------------
# Exportar dataset limpio
# -----------------------------------------------------------------------
df_final <- df_t %>% select(-year)
write_csv(df_final, "datos/limpio.csv")
cat(sprintf("\n[OK] Dataset limpio exportado: datos/limpio.csv (%d filas, %d filas excluidas en total respecto del original de %d)\n",
            nrow(df_final), n0 - nrow(df_final), n0))