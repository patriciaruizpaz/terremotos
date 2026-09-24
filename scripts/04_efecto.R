# =====================================================================
# 04_efecto.R — Tamaño del efecto (Unidad 3)
# Terremoto: profundidad y magnitud de los sismos por región geográfica
# Inferencia Estadística para la Toma de Decisiones — Ruiz Paz / González
#
# Continúa sobre datos/limpio.csv (salida de 01_limpieza.R). No repite
# la limpieza ni los IC de la Unidad 2 (03_ic.R): acá la pregunta ya no
# es "¿hay diferencia / relación?" sino "¿qué tan grande es?".
#
# Igual que en 03_ic.R, la medida principal se calcula de DOS formas:
# con la función de DescTools y a mano con la fórmula, como chequeo
# cruzado antes de reportarla en el trabajo.
#
# Salida: reportes/04_tamanos_efecto.csv (tabla resumen de todo lo
#         calculado, para pasar al trabajo escrito)
# =====================================================================

# install.packages("DescTools")   # instalar una sola vez si falta
library(tidyverse)
library(DescTools)

dir.create("reportes", showWarnings = FALSE)

datos <- read_csv("datos/limpio.csv", show_col_types = FALSE)
datos_mw <- datos %>% filter(familia_magnitud == "familia_mw")

cat(sprintf("Dataset cargado: %d filas totales | %d en familia_mw (%.1f%%)\n",
            nrow(datos), nrow(datos_mw), nrow(datos_mw) / nrow(datos) * 100))

# Etiquetas orientativas de Cohen (1988). Son una referencia, no un
# veredicto: si un efecto es "importante" depende del contexto
# sismológico, y eso hay que justificarlo en el trabajo.
etiqueta_d <- function(d) {
  a <- abs(d)
  if (a < 0.2) "trivial" else if (a < 0.5) "pequeño" else if (a < 0.8) "mediano" else "grande"
}
etiqueta_r <- function(r) {
  a <- abs(r)
  if (a < 0.1) "trivial" else if (a < 0.3) "pequeño" else if (a < 0.5) "mediano" else "grande"
}

separar_grupos <- function(datos_entrada) {
  list(
    cdf   = datos_entrada %>% filter(region == "Cinturón de Fuego") %>% pull(profundidad),
    resto = datos_entrada %>% filter(region == "Resto") %>% pull(profundidad)
  )
}

# -----------------------------------------------------------------------
# PASO 1 - g de Hedges para la diferencia de profundidad media entre
#   regiones, con IC 95%.
#   Es la versión formal de la d de Cohen que 02_decisiones.R dejó como
#   vista previa (Decisión 3): la misma d multiplicada por el factor de
#   corrección J = 1 - 3 / (4(n1 + n2) - 9), que corrige el leve sesgo
#   hacia arriba de d en muestras chicas. Con este n, J es prácticamente
#   1, así que g y d van a coincidir casi exacto — es lo esperable.
#   Orden: Cinturón de Fuego - Resto (positivo = CdF más profundo),
#   igual que el IC de la diferencia en 03_ic.R.
# -----------------------------------------------------------------------
calcular_g_hedges <- function(datos_entrada, etiqueta) {
  g_ <- separar_grupos(datos_entrada)
  n1 <- length(g_$cdf); n2 <- length(g_$resto)

  cat(sprintf("\n[Paso 1 - %s] n CdF = %d | n Resto = %d\n", etiqueta, n1, n2))

  cat("  -- CohenD (DescTools, correct = TRUE -> g de Hedges) --\n")
  print(CohenD(g_$cdf, g_$resto, pooled = TRUE, correct = TRUE, conf.level = 0.95))

  # Cálculo manual (chequeo cruzado)
  desvio_combinado <- sqrt(((n1 - 1) * var(g_$cdf) + (n2 - 1) * var(g_$resto)) / (n1 + n2 - 2))
  d <- (mean(g_$cdf) - mean(g_$resto)) / desvio_combinado
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  ee_g <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
  ic <- g + c(-1, 1) * qnorm(0.975) * ee_g

  cat(sprintf("  -- Manual: d = %.3f | J = %.5f | g = %.3f, IC 95%% [%.3f; %.3f] -> %s\n",
              d, J, g, ic[1], ic[2], etiqueta_d(g)))
  cat("     (si difieren en el 3er decimal es por el método del IC que usa cada uno; con este n es despreciable)\n")

  tibble(analisis = etiqueta, medida = "g de Hedges (desvío combinado)",
         n = n1 + n2, valor = g, li = ic[1], ls = ic[2], magnitud = etiqueta_d(g))
}

res_g_mw    <- calcular_g_hedges(datos_mw, "familia_mw (principal)")
res_g_todas <- calcular_g_hedges(datos, "TODAS las escalas (sensibilidad)")

# -----------------------------------------------------------------------
# PASO 2 - Alternativas a la d "clásica" por la heterocedasticidad
#   Problema: la g del paso 1 usa el desvío COMBINADO, que supone
#   varianzas iguales. Acá la varianza del CdF es ~10x la de Resto
#   (paso 9 de la limpieza, paso 5 de 03_ic.R) y los n son distintos,
#   así que el desvío combinado queda tirado hacia el grupo más grande.
#   Es la misma razón por la que en la Unidad 2 se usó Welch.
#   Por eso se reportan al lado dos estandarizaciones que no suponen
#   varianzas iguales:
#     - d_av: divide por la raíz del promedio SIMPLE de las dos
#       varianzas, así ningún grupo pesa más por tener más casos.
#     - Delta de Glass: estandariza con el desvío de un solo grupo, el de
#       referencia (Resto). Se lee como "a cuántos desvíos típicos de
#       Resto está corrida la media del CdF".
#   Si las tres medidas cuentan la misma historia (mismo signo y
#   categoría parecida), la conclusión no depende de esta elección.
# -----------------------------------------------------------------------
calcular_alternativas_d <- function(datos_entrada, etiqueta) {
  g_ <- separar_grupos(datos_entrada)
  diferencia <- mean(g_$cdf) - mean(g_$resto)

  d_av  <- diferencia / sqrt((var(g_$cdf) + var(g_$resto)) / 2)
  glass <- diferencia / sd(g_$resto)

  cat(sprintf("\n[Paso 2 - %s]\n", etiqueta))
  cat(sprintf("  Desvío CdF = %.1f km | Desvío Resto = %.1f km\n", sd(g_$cdf), sd(g_$resto)))
  cat(sprintf("  d_av (promedio de varianzas) = %.3f -> %s\n", d_av, etiqueta_d(d_av)))
  cat(sprintf("  Delta de Glass (desvío de Resto) = %.3f -> %s\n", glass, etiqueta_d(glass)))

  n <- length(g_$cdf) + length(g_$resto)
  bind_rows(
    tibble(analisis = etiqueta, medida = "d_av (promedio de varianzas)", n = n,
           valor = d_av, li = NA_real_, ls = NA_real_, magnitud = etiqueta_d(d_av)),
    tibble(analisis = etiqueta, medida = "Delta de Glass (ref. Resto)", n = n,
           valor = glass, li = NA_real_, ls = NA_real_, magnitud = etiqueta_d(glass))
  )
}

res_alt_mw <- calcular_alternativas_d(datos_mw, "familia_mw (principal)")
cat("  (IC para d_av y Glass: se pueden sacar por bootstrap en 05_bootstrap.R)\n")

# -----------------------------------------------------------------------
# PASO 3 - Tamaño del efecto no paramétrico: probabilidad de
#   superioridad (A) y delta de Cliff
#   Profundidad tiene asimetría 3,72 y outliers genuinos hasta 700 km
#   (paso 6 de la limpieza). Todas las d se basan en medias y desvíos,
#   que son sensibles justo a eso. Esta medida usa solo el orden
#   (misma lógica que el chequeo con Spearman en 03_ic.R):
#     A = probabilidad de que un sismo del CdF elegido al azar sea más
#         profundo que uno de Resto elegido al azar (empates cuentan 0,5).
#     delta de Cliff = 2A - 1 (va de -1 a 1; 0 = sin diferencia).
#   Es además la más fácil de explicar con palabras en el trabajo.
#   Sale directo del estadístico W de Wilcoxon: A = W / (n1 * n2).
#   Umbrales de Romano et al. (2006) para |delta|: 0,147 / 0,33 / 0,474.
# -----------------------------------------------------------------------
calcular_superioridad <- function(datos_entrada, etiqueta) {
  g_ <- separar_grupos(datos_entrada)
  n1 <- length(g_$cdf); n2 <- length(g_$resto)

  W <- wilcox.test(g_$cdf, g_$resto, exact = FALSE)$statistic
  A <- as.numeric(W) / (n1 * n2)
  delta_cliff <- 2 * A - 1
  a <- abs(delta_cliff)
  etiqueta_cliff <- if (a < 0.147) "trivial" else if (a < 0.33) "pequeño" else if (a < 0.474) "mediano" else "grande"

  cat(sprintf("\n[Paso 3 - %s]\n", etiqueta))
  cat(sprintf("  Probabilidad de superioridad A = %.3f\n", A))
  cat(sprintf("  -> En el %.1f%% de los pares al azar (CdF, Resto), el sismo del CdF es más profundo.\n", A * 100))
  cat(sprintf("  Delta de Cliff = %.3f -> %s (umbrales de Romano et al.)\n", delta_cliff, etiqueta_cliff))

  bind_rows(
    tibble(analisis = etiqueta, medida = "Probabilidad de superioridad (A)", n = n1 + n2,
           valor = A, li = NA_real_, ls = NA_real_, magnitud = NA_character_),
    tibble(analisis = etiqueta, medida = "Delta de Cliff", n = n1 + n2,
           valor = delta_cliff, li = NA_real_, ls = NA_real_, magnitud = etiqueta_cliff)
  )
}

res_sup_mw <- calcular_superioridad(datos_mw, "familia_mw (principal)")

# -----------------------------------------------------------------------
# PASO 4 - La correlación como tamaño del efecto (r y r²)
#   El r de Pearson YA es un tamaño del efecto; su IC se calculó en el
#   paso 3 de 03_ic.R. Acá se agrega lo que falta para interpretarlo:
#   r² = proporción de la variabilidad de profundidad que acompaña a la
#   magnitud. Con un r chico, r² deja claro que casi toda la variación
#   de profundidad se explica por otra cosa — lo que respalda el
#   hallazgo del paso 10 de la limpieza.
#   Ojo: que el IC de r excluya al 0 (con este n, muy probable) NO lo
#   vuelve un efecto relevante. Esa es justamente la diferencia entre
#   significación y tamaño del efecto que trabaja la unidad.
# -----------------------------------------------------------------------
resumir_correlacion <- function(datos_entrada, etiqueta, metodo = "pearson") {
  r <- cor(datos_entrada$magnitud, datos_entrada$profundidad, method = metodo)
  n <- nrow(datos_entrada)
  ic <- CorCI(r, n, conf.level = 0.95)
  nombre <- if (metodo == "pearson") "r de Pearson" else "rho de Spearman"

  cat(sprintf("\n[Paso 4 - %s] %s = %.3f, IC 95%% [%.3f; %.3f] | r² = %.4f (%.2f%% de la variabilidad) -> %s\n",
              etiqueta, nombre, r, ic[["lwr.ci"]], ic[["upr.ci"]], r^2, r^2 * 100, etiqueta_r(r)))

  tibble(analisis = etiqueta, medida = nombre, n = n,
         valor = r, li = ic[["lwr.ci"]], ls = ic[["upr.ci"]], magnitud = etiqueta_r(r))
}

res_r <- bind_rows(
  resumir_correlacion(datos_mw, "familia_mw, general"),
  resumir_correlacion(datos_mw %>% filter(region == "Cinturón de Fuego"), "familia_mw, Cinturón de Fuego"),
  resumir_correlacion(datos_mw %>% filter(region == "Resto"), "familia_mw, Resto"),
  resumir_correlacion(datos_mw, "familia_mw, general", metodo = "spearman"),
  resumir_correlacion(datos, "TODAS las escalas, general (sensibilidad)")
)

# -----------------------------------------------------------------------
# PASO 5 - Sensibilidad: profundidades de relleno (10 / 33 km)
#   El paso 4 de la limpieza flagueó esos registros sin eliminarlos y
#   dejó anotado usarlos como chequeo de sensibilidad. Se repite la g de
#   Hedges sin ellos: si cambia poco, el tamaño del efecto no depende de
#   esos valores fijos que ponen las redes sismográficas.
# -----------------------------------------------------------------------
datos_mw_sin_relleno <- datos_mw %>% filter(!profundidad_valor_relleno)
cat(sprintf("\n[Paso 5] familia_mw sin profundidades de relleno: %d filas (se sacan %d)\n",
            nrow(datos_mw_sin_relleno), nrow(datos_mw) - nrow(datos_mw_sin_relleno)))
res_g_sin_relleno <- calcular_g_hedges(datos_mw_sin_relleno, "familia_mw sin relleno (sensibilidad)")

# -----------------------------------------------------------------------
# PASO 6 - Tabla resumen y formato de reporte
# -----------------------------------------------------------------------
tabla_efectos <- bind_rows(res_g_mw, res_alt_mw, res_sup_mw, res_g_todas, res_g_sin_relleno, res_r) %>%
  mutate(across(c(valor, li, ls), ~ round(.x, 3)))

cat("\n[Paso 6] Tabla resumen de tamaños del efecto:\n")
print(tabla_efectos, n = Inf, width = Inf)
write_csv(tabla_efectos, "reportes/04_tamanos_efecto.csv")
cat("[Paso 6] Tabla guardada en reportes/04_tamanos_efecto.csv\n")

cat("\n[Paso 6] Formato de reporte:\n")
cat("  g = X, IC 95% [LI; LS]  (junto con la diferencia en km del paso 2 de 03_ic.R)\n")
cat("  r = X, IC 95% [LI; LS], r² = X\n")
cat("  Reportar siempre el efecto también en unidades originales (km): 'X desvíos' es abstracto,\n")
cat("  'X km más profundo en promedio' se entiende solo.\n")
cat("  Las etiquetas de Cohen son orientativas: justificar la relevancia desde el contexto, no solo por el umbral.\n")