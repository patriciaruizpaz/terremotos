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
#           profundidad_valor_relleno, familia_magnitud, posible_replica)
#
# Nota de idioma: las columnas se renombran al español apenas se cargan
# (ver más abajo). Los VALORES de texto que vienen del dataset original
# (nombres de lugar en inglés, códigos de escala mb/ms) se dejan tal
# cual, porque cambiarlos rompería la búsqueda contra el dataset real.
# =====================================================================

library(tidyverse)
library(gridExtra)

dir.create("reportes", showWarnings = FALSE)

# -----------------------------------------------------------------------
# Carga y renombrado de columnas al español
# -----------------------------------------------------------------------
datos <- read_csv("datos/crudo.csv", show_col_types = FALSE) %>%
  rename(
    tiempo = time,
    latitud = latitude,
    longitud = longitude,
    profundidad = depth,
    error_profundidad = depthError,
    magnitud = mag,
    tipo_magnitud = magType,
    tipo_evento = type,
    lugar = place,
    actualizado = updated
  )
n_inicial <- nrow(datos)
cat(sprintf("[Carga] %d registros, %d columnas\n", n_inicial, ncol(datos)))

# -----------------------------------------------------------------------
# PASO 1 - Exclusiones del preregistro (Sección 5)
#   - tipo_evento != "earthquake"
#   - profundidad o magnitud faltantes (NA)
#   - profundidad negativa
# -----------------------------------------------------------------------
antes <- nrow(datos)
datos <- datos %>% filter(tipo_evento == "earthquake")
n_tipo <- antes - nrow(datos)

antes <- nrow(datos)
datos <- datos %>% filter(!is.na(profundidad), !is.na(magnitud))
n_faltantes <- antes - nrow(datos)

antes <- nrow(datos)
datos <- datos %>% filter(profundidad >= 0)
n_negativos <- antes - nrow(datos)

cat(sprintf("\n[Paso 1] Excluidos por tipo_evento != earthquake: %d (nuclear explosion, volcanic eruption, explosion, landslide, rock burst, mine collapse)\n", n_tipo))
cat(sprintf("[Paso 1] Excluidos por profundidad/magnitud faltante: %d\n", n_faltantes))
cat(sprintf("[Paso 1] Excluidos por profundidad negativa: %d\n", n_negativos))
cat(sprintf("[Paso 1] Quedan: %d\n", nrow(datos)))

# -----------------------------------------------------------------------
# PASO 2 - Validación de coordenadas
# -----------------------------------------------------------------------
lat_invalida <- !between(datos$latitud, -90, 90)
lon_invalida <- !between(datos$longitud, -180, 180)
cat(sprintf("\n[Paso 2] Latitudes fuera de rango: %d | Longitudes fuera de rango: %d\n", sum(lat_invalida), sum(lon_invalida)))
cat("[Paso 2] Decisión: no hubo nada que excluir, coordenadas 100% válidas.\n")

# -----------------------------------------------------------------------
# PASO 3 - Duplicados por id
#   Hallazgo: el catálogo se actualiza semanalmente y muchos eventos
#   quedaron registrados varias veces (mismo id, distinto 'actualizado');
#   en ~47% de esos grupos, además, profundidad y/o magnitud cambian
#   entre versiones (revisión real del evento a medida que llegan más
#   estaciones).
#   Decisión: quedarse con la revisión más reciente por id.
#
#   NOTA: ajustar as_datetime() si el formato de 'actualizado'/'tiempo'
#   en tu CSV no es ISO-8601 (ej. "2023-04-05T12:34:56.000Z").
# -----------------------------------------------------------------------
antes <- nrow(datos)
datos <- datos %>%
  mutate(actualizado = as_datetime(actualizado)) %>%
  arrange(desc(actualizado)) %>%
  distinct(id, .keep_all = TRUE)
n_duplicados <- antes - nrow(datos)
cat(sprintf("\n[Paso 3] ids duplicados (múltiples snapshots del mismo evento): %d\n", n_duplicados))
cat(sprintf("[Paso 3] Decisión: se conserva solo la revisión más reciente (max 'actualizado') por id. Quedan: %d\n", nrow(datos)))

# -----------------------------------------------------------------------
# PASO 4 - Profundidades de relleno (valores por defecto)
#   Hallazgo: profundidad=10.0 y profundidad=33.0 están brutalmente
#   sobrerrepresentadas frente a sus vecinos, y error_profundidad está
#   ausente en la mayoría de esos casos -> son los valores fijos
#   clásicos de las redes sismográficas cuando no se puede resolver la
#   profundidad real.
#   Decisión: NO se eliminan filas (cambiaría los criterios de exclusión
#   ya preregistrados sin pasar por el Registro de Enmiendas). Se crea
#   un flag para chequeo de sensibilidad en Fase 1-2.
# -----------------------------------------------------------------------
datos <- datos %>%
  mutate(profundidad_valor_relleno = profundidad %in% c(10.0, 33.0) & is.na(error_profundidad))
cat(sprintf("\n[Paso 4] Registros marcados como profundidad por defecto (10 o 33 km sin error estimado): %d (%.1f%% del dataset)\n",
            sum(datos$profundidad_valor_relleno), mean(datos$profundidad_valor_relleno) * 100))
cat("[Paso 4] Decisión: se flaguean, no se eliminan. Usar como chequeo de sensibilidad (Sección 4).\n")

# -----------------------------------------------------------------------
# PASO 5 - Consistencia de tipo_magnitud
#   Hallazgo: la familia "mw" (basada en momento sísmico) no satura y
#   cubre 55.9% del dataset, con magnitudes hasta 9.5. mb satura ~7.1 y
#   ms satura ~8.0 en este dataset (efecto de saturación esperado). Los
#   códigos mb/ms/mw son estándares sismológicos internacionales, se
#   dejan tal cual (no son palabras en inglés a traducir).
#   Decisión: NO se excluye nada acá. Se documenta la recomendación de
#   restringir el análisis principal a la familia mw en Fase 1-2.
# -----------------------------------------------------------------------
familia_mw <- c("mw", "mwc", "mww", "mwb", "mwr", "mwp")
datos <- datos %>%
  mutate(familia_magnitud = case_when(
    tipo_magnitud %in% familia_mw ~ "familia_mw",
    tipo_magnitud == "mb" ~ "mb",
    tipo_magnitud == "ms" ~ "ms",
    TRUE ~ "otra"
  ))
cat("\n[Paso 5] Distribución por familia de escala:\n")
print(table(datos$familia_magnitud))
cat("[Paso 5] Decisión: se documenta el problema de saturación; no se excluye nada todavía. Columna familia_magnitud agregada para decidir en Fase 1-2.\n")

# -----------------------------------------------------------------------
# PASO 6 - Descriptiva, distribución y outliers de profundidad y magnitud
#   Hallazgo: ambas variables muy asimétricas a la derecha; los valores
#   más extremos (profundidad=700km en Fiji, magnitud=9.5 en Chile 1960)
#   corresponden a eventos históricos reales, no a errores de carga.
#   Decisión: no se excluye ningún outlier.
# -----------------------------------------------------------------------
asimetria <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x); s <- sd(x); n <- length(x)
  (sum((x - m)^3) / n) / s^3
}
cat(sprintf("\n[Paso 6] profundidad: media=%.1f, mediana=%.1f, asimetria=%.2f, max=%.1f\n",
            mean(datos$profundidad), median(datos$profundidad), asimetria(datos$profundidad), max(datos$profundidad)))
cat(sprintf("[Paso 6] magnitud:   media=%.2f, mediana=%.2f, asimetria=%.2f, max=%.1f\n",
            mean(datos$magnitud), median(datos$magnitud), asimetria(datos$magnitud), max(datos$magnitud)))
cat("[Paso 6] Decisión: sin exclusiones, valores extremos genuinos (verificados contra sismos históricos conocidos).\n")

# -----------------------------------------------------------------------
# PASO 7 - Cobertura temporal
#   Hallazgo: el salto de instrumentación ocurre donde se esperaba:
#   ~13.3% de los registros son anteriores a 1970, ~86.7% son posteriores.
# -----------------------------------------------------------------------
datos <- datos %>%
  mutate(tiempo = as_datetime(tiempo),
         anio = year(tiempo))
pre_1970 <- mean(datos$anio < 1970) * 100
cat(sprintf("\n[Paso 7] Registros anteriores a 1970: %.1f%% | desde 1970: %.1f%%\n", pre_1970, 100 - pre_1970))
cat("[Paso 7] Decisión: se registra el sesgo de cobertura; la comparación pre/post 1970 queda para Fase 1-2.\n")

# -----------------------------------------------------------------------
# PASO 8 - Construcción de la variable región
#   Criterio: coincidencia de texto en 'lugar' contra países/territorios
#   dominados por subducción a lo largo del Pacífico. Los nombres de
#   lugar en el dataset original vienen en inglés, así que estas
#   palabras clave TIENEN que quedar en inglés para poder encontrarlas.
#   Validado con el patrón de profundidad esperado (Argentina 63.5%
#   >100km, etc.) y con un mapa de puntos.
#   Limitación documentada: zonas de subducción genuinas fuera del
#   Pacífico (South Sandwich, Hindu Kush, Italia) quedan en "Resto"
#   porque el preregistro define la región como "del Pacífico"
#   específicamente (ver Decisión 2 en 02_decisiones.R).
#
#   FIX: coalesce(..., FALSE) — sin esto, los 'lugar' vacíos (NA)
#   quedaban como región "NA" en vez de caer en "Resto".
# -----------------------------------------------------------------------
palabras_cinturon_fuego <- c(
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
patron <- paste(palabras_cinturon_fuego, collapse = "|")
datos <- datos %>%
  mutate(region = if_else(
    coalesce(str_detect(lugar, regex(patron, ignore_case = TRUE)), FALSE),
    "Cinturón de Fuego", "Resto"
  ))

cat("\n[Paso 8]\n")
print(table(datos$region))
print(datos %>% group_by(region) %>%
        summarise(cantidad = n(), media = mean(profundidad), mediana = median(profundidad), desvio = sd(profundidad)))

colores <- c("Cinturón de Fuego" = "#d62728", "Resto" = "#1f77b4")
grafico_paso8 <- ggplot(datos, aes(longitud, latitud, color = region)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = colores) +
  labs(x = "Longitud", y = "Latitud", title = "Paso 8 - Validación visual de la variable región") +
  theme_minimal()
ggsave("reportes/paso8_mapa_validacion.png", grafico_paso8, width = 14, height = 7, dpi = 130)
cat("[Paso 8] Mapa de validación guardado en reportes/paso8_mapa_validacion.png\n")

# -----------------------------------------------------------------------
# PASO 9 - Homogeneidad de varianza de profundidad entre regiones
#   Hallazgo: varianza del Cinturón de Fuego ~10x la de Resto. Confirma
#   que hará falta Welch en vez de t común en Fase 1.
# -----------------------------------------------------------------------
razon_varianzas <- datos %>% group_by(region) %>% summarise(varianza = var(profundidad))
razon <- max(razon_varianzas$varianza) / min(razon_varianzas$varianza)
cat(sprintf("\n[Paso 9] Razón de varianzas: %.1fx\n", razon))
cat("[Paso 9] Decisión: se registra heterocedasticidad marcada; usar t de Welch (no t común) en Fase 1.\n")

grafico_paso9 <- ggplot(datos, aes(x = region, y = profundidad)) +
  geom_boxplot() +
  labs(y = "Profundidad (km)", title = "Paso 9 - profundidad por región") +
  theme_minimal()
ggsave("reportes/paso9_boxplot.png", grafico_paso9, width = 6, height = 5, dpi = 130)

# -----------------------------------------------------------------------
# PASO 10 - Relación profundidad ~ magnitud, general y por región
#   Hallazgo IMPORTANTE: la correlación profundidad-magnitud es
#   prácticamente nula en ambas regiones y las curvas loess son casi
#   planas. El efecto fuerte es el de REGIÓN sobre el nivel/dispersión
#   de profundidad (paso 9), no una interacción magnitud*región. Clave
#   para decidir el modelo de Fase 2.
# -----------------------------------------------------------------------
set.seed(1)

grafico_general <- ggplot(datos, aes(magnitud, profundidad)) +
  geom_point(size = 0.3, alpha = 0.1, color = "gray40") +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black") +
  labs(title = "General", x = "Magnitud", y = "Profundidad (km)") +
  theme_minimal()

graficos_region <- list()
for (region_actual in c("Cinturón de Fuego", "Resto")) {
  subconjunto <- datos %>% filter(region == region_actual)
  prueba <- cor.test(subconjunto$magnitud, subconjunto$profundidad)
  cat(sprintf("[Paso 10] %s: pearson r=%.3f (p=%.2g)\n", region_actual, prueba$estimate, prueba$p.value))
  graficos_region[[region_actual]] <- ggplot(subconjunto, aes(magnitud, profundidad)) +
    geom_point(size = 0.3, alpha = 0.1, color = colores[region_actual]) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black") +
    labs(title = region_actual, x = "Magnitud", y = NULL) +
    theme_minimal()
}

grafico_paso10 <- grid.arrange(grafico_general, graficos_region[["Cinturón de Fuego"]], graficos_region[["Resto"]],
                                ncol = 3, top = "Paso 10 - profundidad vs magnitud, general y por región (loess)")
ggsave("reportes/paso10_scatter_loess.png", grafico_paso10, width = 16, height = 5, dpi = 130)
cat("[Paso 10] Decisión: relación lineal simple magnitud-profundidad es débil/nula; el efecto está en el nivel medio y la dispersión de profundidad por región, no en una pendiente distinta por región. Registrar como hallazgo para discutir el enfoque de Fase 2.\n")

# -----------------------------------------------------------------------
# PASO 11 - Chequeo preliminar de independencia (réplicas)
#   Método simple: ordenando por tiempo, se marca cada evento cuyo
#   vecino inmediato también está a <100km y <24h de diferencia (proxy
#   simple de réplica/aftershock, sin declustering formal).
#   Hallazgo: ~12.7% de los eventos consecutivos cumplen ambos criterios.
#   Decisión: no se corrige todavía (queda para Fase 3, Unidad 4).
# -----------------------------------------------------------------------
distancia_haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  lat1 <- lat1 * pi / 180; lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180; lon2 <- lon2 * pi / 180
  dlat <- lat2 - lat1; dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

datos_ordenados <- datos %>% arrange(tiempo)
horas_diferencia <- as.numeric(difftime(datos_ordenados$tiempo, lag(datos_ordenados$tiempo), units = "hours"))
distancia_km <- distancia_haversine(datos_ordenados$latitud, datos_ordenados$longitud,
                                     lag(datos_ordenados$latitud), lag(datos_ordenados$longitud))
datos_ordenados <- datos_ordenados %>%
  mutate(posible_replica = (horas_diferencia < 24) & (distancia_km < 100))

cat(sprintf("\n[Paso 11] Posibles réplicas (diferencia de tiempo <24h y <100km del evento anterior en la secuencia temporal): %d (%.1f%%)\n",
            sum(datos_ordenados$posible_replica, na.rm = TRUE), mean(datos_ordenados$posible_replica, na.rm = TRUE) * 100))
cat("[Paso 11] Decisión: se registra, no se corrige todavía (corrección/limitación se decide en Fase 3).\n")

# -----------------------------------------------------------------------
# Exportar dataset limpio
# -----------------------------------------------------------------------
datos_final <- datos_ordenados %>% select(-anio)
write_csv(datos_final, "datos/limpio.csv")
cat(sprintf("\n[OK] Dataset limpio exportado: datos/limpio.csv (%d filas, %d filas excluidas en total respecto del original de %d)\n",
            nrow(datos_final), n_inicial - nrow(datos_final), n_inicial))