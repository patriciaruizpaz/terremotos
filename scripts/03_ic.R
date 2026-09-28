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
#
# ADVERTENCIA GENERAL: todos los IC analíticos de este script suponen
# observaciones independientes. Los sismos no lo son (réplicas, paso 11
# de la limpieza), así que estos IC son, si acaso, angostos de más. La
# sensibilidad a esa dependencia se evalúa en 05_bootstrap.R.
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
# PASO 3 (continuación) - IC para la DIFERENCIA de correlaciones entre
#   regiones (Resto - Cinturón de Fuego), sobre familia_mw.
#   Por qué hace falta: la pregunta del proyecto es si la relación
#   profundidad-magnitud CAMBIA según la región. Tener r = 0,007 (CdF,
#   IC que incluye el 0) y r = 0,067 (Resto, IC que excluye el 0) NO
#   responde eso: comparar si cada IC incluye o no el cero es el Error 3
#   de la Guía Teórica (comparar ICs individuales). Hay que estimar la
#   diferencia directamente.
#   Método: IC de Zou (2007) para la diferencia de dos correlaciones
#   independientes, armado con los IC de Fisher de cada una. Se agrega la
#   prueba z clásica de Fisher solo como referencia.
#   Limitaciones: (1) supone dos muestras independientes (las secuencias
#   sísmicas cruzan regiones en muy pocos casos, pero la dependencia
#   entre réplicas sigue estando); (2) una diferencia de r chica con este
#   n puede ser "distinguible de cero" y aun así trivial: interpretar el
#   tamaño, no solo si excluye el 0. La versión bootstrap está en
#   05_bootstrap.R.
# -----------------------------------------------------------------------
ic_diferencia_correlaciones <- function(r1, n1, r2, n2, nivel = 0.95) {
  ic1 <- CorCI(r1, n1, conf.level = nivel)
  ic2 <- CorCI(r2, n2, conf.level = nivel)
  l1 <- ic1[["lwr.ci"]]; u1 <- ic1[["upr.ci"]]
  l2 <- ic2[["lwr.ci"]]; u2 <- ic2[["upr.ci"]]
  dif <- r1 - r2
  c(diferencia = dif,
    lwr.ci = dif - sqrt((r1 - l1)^2 + (u2 - r2)^2),
    upr.ci = dif + sqrt((u1 - r1)^2 + (r2 - l2)^2))
}

datos_mw_cdf <- datos_mw %>% filter(region == "Cinturón de Fuego")
datos_mw_resto <- datos_mw %>% filter(region == "Resto")
r_cdf <- cor(datos_mw_cdf$magnitud, datos_mw_cdf$profundidad); n_cdf <- nrow(datos_mw_cdf)
r_resto <- cor(datos_mw_resto$magnitud, datos_mw_resto$profundidad); n_resto <- nrow(datos_mw_resto)

ic_dif_r <- ic_diferencia_correlaciones(r_resto, n_resto, r_cdf, n_cdf)
z_dif <- (atanh(r_resto) - atanh(r_cdf)) / sqrt(1 / (n_resto - 3) + 1 / (n_cdf - 3))
cat(sprintf("\n[Paso 3 - diferencia de r entre regiones, familia_mw]\n"))
cat(sprintf("  r Resto = %.3f (n = %d) | r Cinturón de Fuego = %.3f (n = %d)\n", r_resto, n_resto, r_cdf, n_cdf))
cat(sprintf("  Diferencia (Resto - CdF) = %.3f, IC 95%% [%.3f; %.3f] (Zou)\n",
            ic_dif_r[["diferencia"]], ic_dif_r[["lwr.ci"]], ic_dif_r[["upr.ci"]]))
cat(sprintf("  Referencia: z de Fisher = %.2f, p = %.2g\n", z_dif, 2 * pnorm(-abs(z_dif))))

# -----------------------------------------------------------------------
# PASO 3 (continuación) - Chequeo de robustez: Spearman
#   La Guía Metodológica (Sección 6) recomienda esto explícitamente
#   cuando hay outliers o no linealidad: "considerar la correlación de
#   Spearman o Kendall". Acá aplica directo: profundidad tiene asimetría
#   3,72 y outliers genuinos hasta 700 km (paso 6 de la limpieza), y
#   Pearson es sensible justo a eso. Spearman usa rangos en vez de
#   valores crudos, así que no lo afectan los extremos de la misma forma.
#   Nota: CorCI(rho, n) aplica la fórmula de Fisher pensada para Pearson;
#   para Spearman el error estándar real es algo mayor (Bonett y Wright,
#   2000), así que este IC es aproximado y levemente angosto. Con rho ~0,01
#   la diferencia es despreciable, pero conviene decirlo si se reporta.
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
# PASO 5 - Verificación de supuestos
#   Homogeneidad de varianza: se confirma acá específicamente sobre
#   familia_mw — el paso 9 de la limpieza la vio sobre el dataset
#   completo, no sobre este subconjunto puntual. (La razón de varianzas
#   es una descripción, no una prueba formal; alcanza para justificar
#   Welch, que de todos modos no exige varianzas iguales.)
#   Forma de la distribución: profundidad es muy asimétrica (asimetría
#   3,72); con n de decenas de miles la media muestral es aproximadamente
#   normal igual (TCL), pero declararlo.
#   INDEPENDENCIA: NO se cumple. ~9% de los sismos de familia_mw son
#   posibles réplicas (paso 11 de la limpieza), y los IC de este script
#   la suponen. No se corrige acá: se declara como limitación y se evalúa
#   con el bootstrap por secuencias de 05_bootstrap.R.
# -----------------------------------------------------------------------
varianzas_mw <- datos_mw %>% group_by(region) %>% summarise(varianza = var(profundidad))
razon_mw <- max(varianzas_mw$varianza) / min(varianzas_mw$varianza)
cat(sprintf("\n[Paso 5] Razón de varianzas dentro de familia_mw: %.1fx\n", razon_mw))
cat("[Paso 5] Decisión: confirma heterocedasticidad también en el subconjunto principal -> Welch queda justificado, no solo en el dataset completo.\n")
cat("[Paso 5] Independencia: no se cumple del todo (réplicas). Los IC de este script deben leerse junto con la sensibilidad de 05_bootstrap.R.\n")

# -----------------------------------------------------------------------
# PASO 6 - Redacción en formato APA
#   Formato de la Guía Teórica, Sección 8. Se arma a mano con los
#   números de los pasos 2 y 3 de arriba; esto queda como recordatorio
#   del formato esperado, no se genera solo.
# -----------------------------------------------------------------------
cat("\n[Paso 6] Formato de reporte (Guía Teórica, Sección 8):\n")
cat("  diferencia de medias = X km, IC 95% [LI; LS]\n")
cat("  r = X, IC 95% [LI; LS]\n")
cat("  diferencia de r = X, IC 95% [LI; LS]\n")
cat("  Evitar Error 2 (nunca 'probabilidad de que el parámetro esté en el IC') y Error 3 (no comparar superposición de IC individuales entre regiones — para eso están el IC de la diferencia del paso 2 y el de la diferencia de r del paso 3).\n")