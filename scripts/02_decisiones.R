# =====================================================================
# Decisiones abiertas del AE — Terremoto: profundidad y magnitud por
# región
# Continúa sobre datos/limpio.csv (salida de 01_limpieza.R)
# No calcula intervalos de confianza formales (eso es Fase 1, Unidad 2)
# — son chequeos de apoyo para las 3 decisiones ya discutidas y resueltas
# por el grupo (ver Registro de Enmiendas del preregistro).
# =====================================================================

library(tidyverse)

df <- read_csv("datos/limpio.csv", show_col_types = FALSE)
cat(sprintf("Dataset limpio cargado: %d filas\n", nrow(df)))

# ------------------------------------------------------------------
# DECISIÓN 1 — magType: mw_family como análisis principal, todas las
# escalas mezcladas como chequeo de sensibilidad
# ------------------------------------------------------------------
mw_family <- c("mw", "mwc", "mww", "mwb", "mwr", "mwp")
df <- df %>% mutate(magType_family = if_else(magType %in% mw_family, "mw_family", "otra"))
df_mw <- df %>% filter(magType_family == "mw_family")

cat("\n=== DECISIÓN 1 — magType ===\n")
cat(sprintf("n total: %d | n mw_family: %d (%.1f%%)\n", nrow(df), nrow(df_mw), nrow(df_mw) / nrow(df) * 100))

reportar_correlacion <- function(data, etiqueta) {
  test <- cor.test(data$mag, data$depth)
  cat(sprintf("\n%s: r(depth,mag)=%.3f (p=%.2g), n=%d\n", etiqueta, test$estimate, test$p.value, nrow(data)))
  for (reg in unique(data$region)) {
    sub <- data %>% filter(region == reg)
    t <- cor.test(sub$mag, sub$depth)
    cat(sprintf("   %s: r=%.3f (p=%.2g), n=%d\n", reg, t$estimate, t$p.value, nrow(sub)))
  }
}
reportar_correlacion(df, "TODAS las escalas (sensibilidad)")
reportar_correlacion(df_mw, "Solo mw_family (principal)")

# ------------------------------------------------------------------
# DECISIÓN 2 — región: mantener "Cinturón de Fuego del Pacífico" tal
# como está preregistrado. Chequeo de sensibilidad con una versión
# ampliada que suma otras zonas de subducción fuera del Pacífico.
# ------------------------------------------------------------------
extra_subduccion <- c("South Sandwich", "Hindu Kush", "Afghanistan", "Italy")
pattern_extra <- paste(extra_subduccion, collapse = "|")

df <- df %>%
  mutate(region_ampliada = if_else(
    str_detect(place, regex(pattern_extra, ignore_case = TRUE)) & region == "Resto",
    "Cinturón de Fuego", region
  ))

n_reclasificados <- sum(df$region != df$region_ampliada)
cat("\n=== DECISIÓN 2 — definición de región ===\n")
cat(sprintf("Eventos reclasificados si se amplía la definición: %d (%.2f%% del dataset)\n",
            n_reclasificados, n_reclasificados / nrow(df) * 100))

for (col in c("region", "region_ampliada")) {
  etiqueta <- if (col == "region") "Definición preregistrada (Pacífico)" else "Definición ampliada (sensibilidad)"
  cat(sprintf("\n%s:\n", etiqueta))
  print(df %>% group_by(.data[[col]]) %>%
          summarise(count = n(), mean = mean(depth), median = median(depth), sd = sd(depth)))
}

# ------------------------------------------------------------------
# DECISIÓN 3 — vista previa informal (NO reemplaza el cálculo formal de
# Fase 1-2, Unidad 2/3): sobre mw_family, ¿cuánto pesa cada pieza de la
# hipótesis original — la relación r, o la diferencia d?
# ------------------------------------------------------------------
cat("\n=== DECISIÓN 3 — vista previa (r vs. d, sobre mw_family) ===\n")

test_total <- cor.test(df_mw$mag, df_mw$depth)
cat(sprintf("Relación depth~mag (r de Pearson, mw_family): r=%.3f (p=%.2g)\n",
            test_total$estimate, test_total$p.value))

g1 <- df_mw %>% filter(region == "Cinturón de Fuego") %>% pull(depth)
g2 <- df_mw %>% filter(region == "Resto") %>% pull(depth)
n1 <- length(g1); n2 <- length(g2)
pooled_sd <- sqrt(((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) / (n1 + n2 - 2))
cohen_d <- (mean(g1) - mean(g2)) / pooled_sd

cat(sprintf("Diferencia de profundidad media (Cohen's d, mw_family): d=%.3f\n", cohen_d))
cat(sprintf("   Media Cinturón de Fuego: %.1f km (n=%d)\n", mean(g1), n1))
cat(sprintf("   Media Resto: %.1f km (n=%d)\n", mean(g2), n2))
cat("(d sin corrección de Hedges todavía — esa corrección se calcula formalmente en Fase 2, Unidad 3)\n")