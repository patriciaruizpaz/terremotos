# =====================================================================
# Decisiones abiertas del AE — Terremoto: profundidad y magnitud por
# región
# Continúa sobre datos/limpio.csv (salida de 01_limpieza.R)
# No calcula intervalos de confianza formales (eso es Fase 1, Unidad 2)
# — son chequeos de apoyo para las 3 decisiones ya discutidas y resueltas
# por el grupo (ver Registro de Enmiendas del preregistro).
# =====================================================================

library(tidyverse)

datos <- read_csv("datos/limpio.csv", show_col_types = FALSE)
cat(sprintf("Dataset limpio cargado: %d filas\n", nrow(datos)))

# ------------------------------------------------------------------
# DECISIÓN 1 — familia_magnitud: familia_mw como análisis principal,
# todas las escalas mezcladas como chequeo de sensibilidad.
# (la columna familia_magnitud ya viene calculada de 01_limpieza.R,
# no hace falta recalcularla)
# ------------------------------------------------------------------
datos_mw <- datos %>% filter(familia_magnitud == "familia_mw")

cat("\n=== DECISIÓN 1 — familia_magnitud ===\n")
cat(sprintf("n total: %d | n familia_mw: %d (%.1f%%)\n", nrow(datos), nrow(datos_mw), nrow(datos_mw) / nrow(datos) * 100))

reportar_correlacion <- function(datos_entrada, etiqueta) {
  prueba <- cor.test(datos_entrada$magnitud, datos_entrada$profundidad)
  cat(sprintf("\n%s: r(profundidad,magnitud)=%.3f (p=%.2g), n=%d\n", etiqueta, prueba$estimate, prueba$p.value, nrow(datos_entrada)))
  for (region_actual in unique(datos_entrada$region)) {
    subconjunto <- datos_entrada %>% filter(region == region_actual)
    prueba_region <- cor.test(subconjunto$magnitud, subconjunto$profundidad)
    cat(sprintf("   %s: r=%.3f (p=%.2g), n=%d\n", region_actual, prueba_region$estimate, prueba_region$p.value, nrow(subconjunto)))
  }
}
reportar_correlacion(datos, "TODAS las escalas (sensibilidad)")
reportar_correlacion(datos_mw, "Solo familia_mw (principal)")

# ------------------------------------------------------------------
# DECISIÓN 2 — región: mantener "Cinturón de Fuego del Pacífico" tal
# como está preregistrado. Chequeo de sensibilidad con una versión
# ampliada que suma otras zonas de subducción fuera del Pacífico.
# ------------------------------------------------------------------
zonas_subduccion_extra <- c("South Sandwich", "Hindu Kush", "Afghanistan", "Italy")
patron_extra <- paste(zonas_subduccion_extra, collapse = "|")

datos <- datos %>%
  mutate(region_ampliada = if_else(
    coalesce(str_detect(lugar, regex(patron_extra, ignore_case = TRUE)), FALSE) & region == "Resto",
    "Cinturón de Fuego", region
  ))

n_reclasificados <- sum(datos$region != datos$region_ampliada)
cat("\n=== DECISIÓN 2 — definición de región ===\n")
cat(sprintf("Eventos reclasificados si se amplía la definición: %d (%.2f%% del dataset)\n",
            n_reclasificados, n_reclasificados / nrow(datos) * 100))

for (columna in c("region", "region_ampliada")) {
  etiqueta <- if (columna == "region") "Definición preregistrada (Pacífico)" else "Definición ampliada (sensibilidad)"
  cat(sprintf("\n%s:\n", etiqueta))
  print(datos %>% group_by(.data[[columna]]) %>%
          summarise(cantidad = n(), media = mean(profundidad), mediana = median(profundidad), desvio = sd(profundidad)))
}

# ------------------------------------------------------------------
# DECISIÓN 3 — vista previa informal (NO reemplaza el cálculo formal de
# Fase 1-2, Unidad 2/3): sobre familia_mw, ¿cuánto pesa cada pieza de la
# hipótesis original — la relación r, o la diferencia d?
# ------------------------------------------------------------------
cat("\n=== DECISIÓN 3 — vista previa (r vs. d, sobre familia_mw) ===\n")

prueba_total <- cor.test(datos_mw$magnitud, datos_mw$profundidad)
cat(sprintf("Relación profundidad~magnitud (r de Pearson, familia_mw): r=%.3f (p=%.2g)\n",
            prueba_total$estimate, prueba_total$p.value))

grupo1 <- datos_mw %>% filter(region == "Cinturón de Fuego") %>% pull(profundidad)
grupo2 <- datos_mw %>% filter(region == "Resto") %>% pull(profundidad)
n_grupo1 <- length(grupo1); n_grupo2 <- length(grupo2)
desvio_combinado <- sqrt(((n_grupo1 - 1) * var(grupo1) + (n_grupo2 - 1) * var(grupo2)) / (n_grupo1 + n_grupo2 - 2))
d_cohen <- (mean(grupo1) - mean(grupo2)) / desvio_combinado

cat(sprintf("Diferencia de profundidad media (d de Cohen, familia_mw): d=%.3f\n", d_cohen))
cat(sprintf("   Media Cinturón de Fuego: %.1f km (n=%d)\n", mean(grupo1), n_grupo1))
cat(sprintf("   Media Resto: %.1f km (n=%d)\n", mean(grupo2), n_grupo2))
cat("(d sin corrección de Hedges todavía — esa corrección se calcula formalmente en Fase 2, Unidad 3)\n")