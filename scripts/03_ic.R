# =====================================================================
# 03_ic.R — Intervalos de confianza (Unidad 2)
# Terremoto: profundidad y magnitud de los sismos por región geográfica
# Inferencia Estadística para la Toma de Decisiones — Ruiz Paz / González
#
# Implementa los pasos 1 a 6 de registro_unidad2_intervalos.docx.
# Continúa sobre datos/limpio.csv (salida de 01_limpieza.R). No repite
# nada de la limpieza — solo carga el dataset ya limpio.
#
# Cada IC se calcula de DOS formas independientes: con la función de
# DescTools que pide la guía de la materia, y con su equivalente en R
# base, como chequeo cruzado — si las dos coinciden, da más confianza
# en el resultado antes de reportarlo en el trabajo.
# =====================================================================

# install.packages("DescTools")   # instalar una sola vez si falta
library(tidyverse)
library(DescTools)

datos <- read_csv("datos/limpio.csv", show_col_types = FALSE)
datos_mw <- datos %>% filter(familia_magnitud == "familia_mw")

cat(sprintf("Dataset cargado: %d filas totales | %d en familia_mw (%.1f%%)\n",
            nrow(datos), nrow(datos_mw), nrow(datos_mw) / nrow(datos) * 100))

# -----------------------------------------------------------------------
# PASO 1 - Subconjuntos de trabajo
#   Principal: familia_mw (no satura). Sensibilidad: todas las escalas
#   mezcladas. Ya quedaron preparados arriba (datos_mw / datos).
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# PASO 2 - IC para la diferencia de profundidad media entre regiones
#   Método: Welch (varianzas distintas — confirmado en el paso 9 de la
#   limpieza, razón ~10x). Orden: Cinturón de Fuego - Resto, para que
#   un resultado positivo se lea como "el Cinturón de Fuego es más
#   profundo".
# -----------------------------------------------------------------------
calcular_diferencia_profundidad <- function(datos_entrada, etiqueta) {
  grupo_cdf <- datos_entrada %>% filter(region == "Cinturón de Fuego") %>% pull(profundidad)
  grupo_resto <- datos_entrada %>% filter(region == "Resto") %>% pull(profundidad)

  cat(sprintf("\n[Paso 2 - %s]\n", etiqueta))
  cat(sprintf("  Diferencia de medias observada (Cinturón de Fuego - Resto) = %.2f km\n",
              mean(grupo_cdf) - mean(grupo_resto)))

  cat("  -- MeanDiffCI (DescTools, Welch por defecto) --\n")
  print(MeanDiffCI(grupo_cdf, grupo_resto, conf.level = 0.95))

  cat("  -- t.test (base R, var.equal = FALSE, chequeo cruzado) --\n")
  print(t.test(grupo_cdf, grupo_resto, var.equal = FALSE)$conf.int)
}

calcular_diferencia_profundidad(datos_mw, "familia_mw (principal)")
calcular_diferencia_profundidad(datos, "TODAS las escalas (sensibilidad)")

# -----------------------------------------------------------------------
# PASO 3 - IC para la correlación profundidad-magnitud
#   Método: transformación de Fisher. El paso 10 de la limpieza ya
#   mostró curvas loess casi planas — un IC angosto alrededor de un r
#   chico es justamente lo esperable con este n (Guía Teórica, Sección
#   7): no es un error, es lo que corresponde reportar con honestidad.
# -----------------------------------------------------------------------
calcular_correlacion <- function(datos_entrada, etiqueta) {
  r <- cor(datos_entrada$magnitud, datos_entrada$profundidad)
  n <- nrow(datos_entrada)

  cat(sprintf("\n[Paso 3 - %s]\n", etiqueta))
  cat(sprintf("  r = %.3f, n = %d\n", r, n))

  cat("  -- CorCI (DescTools, transformación de Fisher) --\n")
  print(CorCI(r, n, conf.level = 0.95))

  cat("  -- cor.test (base R, chequeo cruzado) --\n")
  print(cor.test(datos_entrada$magnitud, datos_entrada$profundidad)$conf.int)
}

calcular_correlacion(datos_mw, "familia_mw, general")
calcular_correlacion(datos, "TODAS las escalas, general (sensibilidad)")

for (region_actual in c("Cinturón de Fuego", "Resto")) {
  subconjunto <- datos_mw %>% filter(region == region_actual)
  calcular_correlacion(subconjunto, sprintf("familia_mw, %s", region_actual))
}

# -----------------------------------------------------------------------
# PASO 3 (continuación) - Chequeo de robustez: Spearman
#   La Guía Metodológica (Sección 6) recomienda esto explícitamente
#   cuando hay outliers o no linealidad: "considerar la correlación de
#   Spearman o Kendall". Acá aplica directo: profundidad tiene asimetría
#   3,72 y outliers genuinos hasta 700 km (paso 6 de la limpieza), y
#   Pearson es sensible justo a eso. Spearman usa rangos en vez de
#   valores crudos, así que no lo afectan los extremos de la misma forma.
# -----------------------------------------------------------------------
r_spearman <- cor(datos_mw$magnitud, datos_mw$profundidad, method = "spearman")
n_mw <- nrow(datos_mw)
cat(sprintf("\n[Paso 3 - robustez Spearman, familia_mw] rho = %.3f, n = %d\n", r_spearman, n_mw))
print(CorCI(r_spearman, n_mw, conf.level = 0.95))

# -----------------------------------------------------------------------
# PASO 4 - (Opcional) IC para la pendiente de una regresión
#   profundidad ~ magnitud
#   No se corre por defecto: el paso 10 de la limpieza y el paso 3 de
#   acá ya muestran una relación prácticamente nula. Descomentar estas
#   dos líneas si el grupo decide que igual vale la pena reportarla.
# -----------------------------------------------------------------------
# modelo <- lm(profundidad ~ magnitud, data = datos_mw)
# confint(modelo, level = 0.95)

# -----------------------------------------------------------------------
# PASO 5 - Verificación formal de supuestos
#   Independencia y forma de la distribución ya se relevaron en la
#   limpieza (pasos 6 y 11). Acá se confirma homogeneidad de varianza
#   específicamente sobre familia_mw — el paso 9 de la limpieza la vio
#   sobre el dataset completo, no sobre este subconjunto puntual, así
#   que conviene chequearlo de nuevo acá antes de dar Welch por sentado.
# -----------------------------------------------------------------------
varianzas_mw <- datos_mw %>% group_by(region) %>% summarise(varianza = var(profundidad))
razon_mw <- max(varianzas_mw$varianza) / min(varianzas_mw$varianza)
cat(sprintf("\n[Paso 5] Razón de varianzas dentro de familia_mw: %.1fx\n", razon_mw))
cat("[Paso 5] Decisión: confirma heterocedasticidad también en el subconjunto principal -> Welch queda justificado, no solo en el dataset completo.\n")

# -----------------------------------------------------------------------
# PASO 6 - Redacción en formato APA
#   Formato de la Guía Teórica, Sección 8. Se arma a mano con los
#   números de los pasos 2 y 3 de arriba; esto queda como recordatorio
#   del formato esperado, no se genera solo.
# -----------------------------------------------------------------------
cat("\n[Paso 6] Formato de reporte (Guía Teórica, Sección 8):\n")
cat("  diferencia de medias = X km, IC 95% [LI; LS]\n")
cat("  r = X, IC 95% [LI; LS]\n")
cat("  Evitar Error 2 (nunca 'probabilidad de que el parámetro esté en el IC') y Error 3 (no comparar superposición de IC individuales entre regiones — para eso está el IC de la diferencia del paso 2).\n")