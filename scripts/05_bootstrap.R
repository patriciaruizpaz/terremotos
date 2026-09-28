# =====================================================================
# 05_bootstrap.R — Bootstrap y permutaciones (Unidad 4)
# Terremoto: profundidad y magnitud de los sismos por región geográfica
# Inferencia Estadística para la Toma de Decisiones — Ruiz Paz / González
#
# VERSIÓN RÁPIDA (v2). Qué hace y qué cambió respecto de la versión anterior:
#   1) Sin BCa: boot.ci(type = "bca") hacía un jackknife de n = 56.993
#      evaluaciones y no terminaba. Se usa IC percentil.
#   2) Bootstrap ingenuo y por secuencia con sumas suficientes por fila /
#      por secuencia: respeta las repeticiones de secuencias (el %in% de
#      la primera versión no las repetía) y es rápido.
#   3) NUEVO: el remuestreo es propio (rmultinom) en vez de boot(): boot()
#      arma una matriz R x n de índices (5000 x 57 mil = más de 1 GB de
#      RAM). Es el mismo bootstrap ordinario, con memoria constante.
#   4) NUEVO: el mismo remuestreo devuelve TODAS las medidas de una vez:
#      r general y por región, diferencia de r entre regiones (la parte
#      "¿la relación cambia según la región?" de la pregunta), diferencia
#      de medias en km, d de Cohen, g de Hedges, d_av y Delta de Glass
#      (los IC de d_av y Glass que 04_efecto.R dejaba pendientes).
#   5) NUEVO: la síntesis compara ingenuo vs. por secuencia con la razón
#      de anchos (cuánto se ensancha el IC al respetar las réplicas) y si
#      ambos IC excluyen o no el 0.
#   6) IC por inversión con el efecto observado ajustado por delta
#      (corregido en la versión anterior; se mantiene).
#   7) p permutacional: nunca se reporta p = 0, sino "menor que 1/(B+1)".
#   8) CACHÉ: la primera vez calcula y guarda en reportes/*_v2.rds; las
#      siguientes solo carga y reporta. El nombre cambió (_v2) para no
#      cargar por error la caché de la versión anterior, que tiene otra
#      estructura.
#
# Ya no necesita el paquete effectsize: d y g se calculan con sumas.
# =====================================================================

library(tidyverse)   # único paquete necesario (boot y effectsize ya no se usan)

# ---------------- Parámetros ----------------
MODO_PRUEBA <- FALSE   # TRUE = corrida corta para probar que anda
RECALCULAR  <- FALSE   # TRUE = ignora la caché y vuelve a calcular todo

R_BOOT <- if (MODO_PRUEBA) 200   else 5000    # réplicas bootstrap
B_PERM <- if (MODO_PRUEBA) 500   else 10000   # permutaciones (libre y restringida)
B_INV  <- if (MODO_PRUEBA) 100   else 1000    # permutaciones por punto de grilla (inversión)
ARCHIVO_CACHE <- if (MODO_PRUEBA) "reportes/05_bootstrap_prueba_v2.rds" else
  "reportes/05_bootstrap_resultados_v2.rds"
dir.create("reportes", showWarnings = FALSE)

set.seed(2026)  # fijada una sola vez, al principio

datos <- read_csv("datos/limpio.csv", show_col_types = FALSE)
datos_mw <- datos %>% filter(familia_magnitud == "familia_mw")

cat(sprintf("Dataset cargado: %d filas totales | %d en familia_mw (%.1f%%)\n",
            nrow(datos), nrow(datos_mw), nrow(datos_mw) / nrow(datos) * 100))

# -----------------------------------------------------------------------
# Haversine (misma fórmula que en 01_limpieza.R, para que corra solo)
# -----------------------------------------------------------------------
distancia_haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  lat1 <- lat1 * pi / 180; lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180; lon2 <- lon2 * pi / 180
  dlat <- lat2 - lat1; dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

# -----------------------------------------------------------------------
# PASO 1 - Estructura de secuencias sísmicas (recalculada dentro de familia_mw)
#   Limitación (declarar en el reporte): cada sismo se compara solo con el
#   inmediato anterior en el tiempo GLOBAL, así que dos secuencias
#   simultáneas en lugares distintos se cortan entre sí, y una secuencia
#   puede mezclar regiones. Es un proxy, no un declustering formal.
# -----------------------------------------------------------------------
construir_secuencias <- function(df) {
  df <- df %>% arrange(tiempo)
  horas_diferencia <- as.numeric(difftime(df$tiempo, lag(df$tiempo), units = "hours"))
  distancia_km <- distancia_haversine(df$latitud, df$longitud, lag(df$latitud), lag(df$longitud))
  es_replica <- (horas_diferencia < 24) & (distancia_km < 100)
  es_replica[is.na(es_replica)] <- FALSE
  df$secuencia_id <- cumsum(!es_replica)
  df
}

datos_mw <- construir_secuencias(datos_mw)
n_secuencias <- n_distinct(datos_mw$secuencia_id)
cat(sprintf("\n[Paso 1] familia_mw: %d sismos agrupados en %d secuencias (%.2f sismos por secuencia en promedio)\n",
            nrow(datos_mw), n_secuencias, nrow(datos_mw) / n_secuencias))
tam_sec <- table(datos_mw$secuencia_id)
cat(sprintf("[Paso 1] Secuencias de 1 sismo: %d (%.1f%%) | de 2 o más: %d | tamaño máximo: %d\n",
            sum(tam_sec == 1), mean(tam_sec == 1) * 100, sum(tam_sec >= 2), max(tam_sec)))

# -----------------------------------------------------------------------
# Datos auxiliares (vectores, mucho más rápidos que data frames)
# -----------------------------------------------------------------------
prof   <- datos_mw$profundidad
es_cdf <- datos_mw$region == "Cinturón de Fuego"
N  <- length(prof); n1 <- sum(es_cdf); S <- sum(prof)
dif_obs <- mean(prof[es_cdf]) - mean(prof[!es_cdf])

# Sumas suficientes por sismo: n, Σx, Σy, Σxy, Σx², Σy² para (a) todos,
# (b) Cinturón de Fuego (_c) y (c) Resto (_r). x e y se centran con la
# media global por estabilidad numérica (r y d no cambian al centrar).
x <- datos_mw$magnitud - mean(datos_mw$magnitud)
y <- prof - mean(prof)
sumas_por_grupo <- function(mascara, sufijo) {
  m <- cbind(n  = as.numeric(mascara),
             x  = ifelse(mascara, x, 0),      y  = ifelse(mascara, y, 0),
             xy = ifelse(mascara, x * y, 0),
             xx = ifelse(mascara, x^2, 0),    yy = ifelse(mascara, y^2, 0))
  colnames(m) <- paste0(colnames(m), sufijo)
  m
}
M_sismo <- cbind(sumas_por_grupo(rep(TRUE, N), ""),
                 sumas_por_grupo(es_cdf, "_c"),
                 sumas_por_grupo(!es_cdf, "_r"))
M_secuencia <- rowsum(M_sismo, datos_mw$secuencia_id)   # una fila por secuencia

r_desde_sumas <- function(n, sx, sy, sxy, sxx, syy) {
  (sxy / n - (sx / n) * (sy / n)) /
    sqrt((sxx / n - (sx / n)^2) * (syy / n - (sy / n)^2))
}

# Todas las medidas a partir de sumas ponderadas por las frecuencias f.
# Convención de signo: diferencias como Cinturón de Fuego - Resto (positivo =
# CdF más profundo), igual que en 03_ic.R y 04_efecto.R. La diferencia de r
# se define al revés (Resto - CdF) porque r_Resto > r_CdF.
estadisticos <- function(M, f) {
  s <- as.numeric(f %*% M); names(s) <- colnames(M)
  r_gen <- r_desde_sumas(s[["n"]],   s[["x"]],   s[["y"]],   s[["xy"]],   s[["xx"]],   s[["yy"]])
  r_c   <- r_desde_sumas(s[["n_c"]], s[["x_c"]], s[["y_c"]], s[["xy_c"]], s[["xx_c"]], s[["yy_c"]])
  r_r   <- r_desde_sumas(s[["n_r"]], s[["x_r"]], s[["y_r"]], s[["xy_r"]], s[["xx_r"]], s[["yy_r"]])
  m_c <- s[["y_c"]] / s[["n_c"]];  m_r <- s[["y_r"]] / s[["n_r"]]
  v_c <- (s[["yy_c"]] - s[["n_c"]] * m_c^2) / (s[["n_c"]] - 1)
  v_r <- (s[["yy_r"]] - s[["n_r"]] * m_r^2) / (s[["n_r"]] - 1)
  n_c <- s[["n_c"]]; n_r <- s[["n_r"]]
  sp   <- sqrt(((n_c - 1) * v_c + (n_r - 1) * v_r) / (n_c + n_r - 2))
  dif  <- m_c - m_r
  d    <- dif / sp
  J    <- 1 - 3 / (4 * (n_c + n_r) - 9)
  c(r = r_gen, r_cdf = r_c, r_resto = r_r, dif_r = r_r - r_c,
    dif_km = dif, d = d, g = J * d,
    d_av = dif / sqrt((v_c + v_r) / 2), glass = dif / sqrt(v_r))
}

etiquetas_medidas <- c(r = "r general", r_cdf = "r Cinturón de Fuego", r_resto = "r Resto",
                       dif_r = "diferencia de r (Resto - CdF)", dif_km = "diferencia de medias (km)",
                       d = "d de Cohen", g = "g de Hedges",
                       d_av = "d_av (promedio de varianzas)", glass = "Delta de Glass (ref. Resto)")

# Bootstrap ordinario con frecuencias: cada réplica reparte K extracciones
# entre las K filas de M (multinomial) y recalcula todas las medidas.
# Con M_sismo remuestrea sismos sueltos (ingenuo); con M_secuencia
# remuestrea secuencias completas (respeta la dependencia por réplicas).
bootstrap_pesos <- function(M, R) {
  K <- nrow(M)
  t0 <- estadisticos(M, rep(1, K))
  out <- matrix(NA_real_, R, length(t0), dimnames = list(NULL, names(t0)))
  for (b in seq_len(R)) out[b, ] <- estadisticos(M, as.numeric(rmultinom(1, K, rep(1, K))))
  out
}

# Chequeo: tiene que reproducir lo ya calculado en 03_ic.R y 04_efecto.R
chk <- estadisticos(M_sismo, rep(1, nrow(M_sismo)))
cat(sprintf("\n[Chequeo] r = %.3f (esperado 0.028) | g = %.3f (esperado 0.487) | dif = %.2f km (esperado 48.89)\n",
            chk[["r"]], chk[["g"]], chk[["dif_km"]]))
cat(sprintf("[Chequeo] d_av = %.3f (esperado 0.564) | Glass = %.3f (esperado 1.295) | r CdF = %.3f (0.007) | r Resto = %.3f (0.067)\n",
            chk[["d_av"]], chk[["glass"]], chk[["r_cdf"]], chk[["r_resto"]]))

# =======================================================================
# CÁLCULO PESADO (solo si no hay caché)
# =======================================================================
if (!RECALCULAR && file.exists(ARCHIVO_CACHE)) {
  cat(sprintf("\n[Caché] Cargando resultados guardados: %s (para recalcular: RECALCULAR <- TRUE)\n", ARCHIVO_CACHE))
  res <- readRDS(ARCHIVO_CACHE)
} else {
  cat("\n[Cálculo] No hay caché: calculando (una sola vez)...\n")
  t0 <- Sys.time()
  seg <- function() as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # --- PASOS 2 y 4: bootstrap de todas las medidas, ingenuo vs. por secuencia ---
  b_ingenuo   <- bootstrap_pesos(M_sismo,     R_BOOT)
  b_secuencia <- bootstrap_pesos(M_secuencia, R_BOOT)
  cat(sprintf("  bootstrap listo (%.0f s)\n", seg()))

  # --- PASO 3: permutación libre (a nivel sismo) ---
  perm_libre <- replicate(B_PERM, {
    s1 <- sum(prof[sample.int(N, n1)])
    s1 / n1 - (S - s1) / (N - n1)
  })

  # --- PASO 3: permutación restringida (a nivel secuencia) ---
  seq_tab <- datos_mw %>%
    group_by(secuencia_id) %>%
    summarise(n = n(), s = sum(profundidad),
              cdf = first(region) == "Cinturón de Fuego",
              mixta = n_distinct(region) > 1, .groups = "drop")
  perm_restringida <- replicate(B_PERM, {
    lab <- sample(seq_tab$cdf)
    n_c <- sum(seq_tab$n[lab]); s_c <- sum(seq_tab$s[lab])
    s_c / n_c - (S - s_c) / (N - n_c)
  })
  cat(sprintf("  permutaciones listas (%.0f s)\n", seg()))

  # --- PASO 3: IC por inversión del test (versión libre) ---
  # OJO: la permutación supone que, bajo H0, las dos regiones tienen la MISMA
  # distribución; con varianzas ~10x distintas y n distintos eso no es
  # cierto para un test de medias. Por eso este IC es un complemento, y el
  # IC bootstrap de la diferencia de medias (más abajo) es el de referencia.
  grilla <- seq(dif_obs - 3, dif_obs + 3, by = 0.25)
  p_grilla <- sapply(grilla, function(delta) {
    prof_aj <- prof
    prof_aj[es_cdf] <- prof_aj[es_cdf] - delta   # H0: diferencia = delta
    S_aj <- sum(prof_aj)
    dif_aj <- dif_obs - delta                     # efecto observado bajo esa H0
    dist_nula <- replicate(B_INV, {
      s1 <- sum(prof_aj[sample.int(N, n1)])
      s1 / n1 - (S_aj - s1) / (N - n1)
    })
    (1 + sum(abs(dist_nula) >= abs(dif_aj))) / (B_INV + 1)
  })
  cat(sprintf("  inversión lista (%.0f s)\n", seg()))

  res <- list(t0 = chk, b_ingenuo = b_ingenuo, b_secuencia = b_secuencia,
              perm_libre = perm_libre, perm_restringida = perm_restringida,
              n_secuencias_mixtas = sum(seq_tab$mixta),
              grilla = grilla, p_grilla = p_grilla,
              parametros = list(R_BOOT = R_BOOT, B_PERM = B_PERM, B_INV = B_INV))
  saveRDS(res, ARCHIVO_CACHE)
  cat(sprintf("[Cálculo] Guardado en %s\n", ARCHIVO_CACHE))
}

# =======================================================================
# REPORTE (rápido; corre siempre)
# =======================================================================
ic_perc <- function(t, nivel = 0.95) unname(quantile(t, c((1 - nivel) / 2, 1 - (1 - nivel) / 2)))

tabla <- map_dfr(c("ingenuo", "secuencia"), function(nm) {
  b <- res[[paste0("b_", nm)]]
  map_dfr(colnames(b), function(m) {
    ic <- ic_perc(b[, m])
    tibble(esquema = nm, medida = etiquetas_medidas[[m]], estimado = res$t0[[m]],
           ee = sd(b[, m]), li = ic[1], ls = ic[2], ancho = ic[2] - ic[1],
           excluye_cero = (ic[1] > 0) | (ic[2] < 0))
  })
})

mostrar <- function(claves) {
  for (m in claves) {
    for (nm in c("ingenuo", "secuencia")) {
      fila <- tabla %>% filter(medida == etiquetas_medidas[[m]], esquema == nm)
      cat(sprintf("  %-28s | %-9s: %8.4f | EE = %.4f | IC 95%% [%.4f; %.4f]\n",
                  etiquetas_medidas[[m]], nm, fila$estimado, fila$ee, fila$li, fila$ls))
    }
  }
}
fmt_p <- function(p, B) if (p <= 1 / (B + 1)) sprintf("< %.5f (ninguna de las %d permutaciones alcanzó la diferencia observada)", 1 / (B + 1), B) else sprintf("= %.4f", p)

# --- PASO 2: correlación (general, por región y diferencia entre regiones) ---
cat("\n[Paso 2 - correlación profundidad-magnitud] (IC percentil, R =", res$parametros$R_BOOT, ")\n")
mostrar(c("r", "r_cdf", "r_resto", "dif_r"))
cat("  Analítico (Fisher, 03_ic.R): r = 0.028, IC 95% [0.020; 0.036] | CdF 0.007 [-0.003; 0.017] | Resto 0.067 [0.052; 0.082]\n")
cat("  Diferencia de r (Resto - CdF): comparar con el IC de Zou de 03_ic.R\n")

# --- PASO 3: permutaciones y diferencia de medias ---
B <- res$parametros$B_PERM
p_libre <- (1 + sum(abs(res$perm_libre) >= abs(dif_obs))) / (B + 1)
p_restr <- (1 + sum(abs(res$perm_restringida) >= abs(dif_obs))) / (B + 1)
cat(sprintf("\n[Paso 3] Diferencia observada (Cinturón de Fuego - Resto) = %.2f km\n", dif_obs))
cat(sprintf("  Permutación libre:       p %s\n", fmt_p(p_libre, B)))
cat(sprintf("  Permutación restringida: p %s (secuencias con sismos de ambas regiones: %d)\n",
            fmt_p(p_restr, B), res$n_secuencias_mixtas))
cat(sprintf("  Diferencia más extrema en las permutaciones: libre %.2f km | restringida %.2f km\n",
            max(abs(res$perm_libre)), max(abs(res$perm_restringida))))
ok <- res$p_grilla > 0.05
if (any(ok)) {
  ic_inv <- range(res$grilla[ok])
  cat(sprintf("  IC 95%% por inversión (libre): [%.2f; %.2f] (resolución de grilla 0.25 km)\n", ic_inv[1], ic_inv[2]))
  if (ic_inv[1] == min(res$grilla) || ic_inv[2] == max(res$grilla))
    cat("  AVISO: el IC toca el borde de la grilla; ampliar la grilla.\n")
} else cat("  AVISO: ningún punto de la grilla con p > 0.05; ampliar la grilla.\n")
cat("  IC bootstrap de la diferencia de medias (no supone varianzas iguales):\n")
mostrar("dif_km")
cat("  Analítico (Welch, 03_ic.R): 48.89 km, IC 95% [47.62; 50.17]\n")

# --- PASO 4: tamaños de efecto ---
cat("\n[Paso 4 - tamaños de efecto estandarizados] (IC percentil)\n")
mostrar(c("d", "g", "d_av", "glass"))
cat("  Analítico (04_efecto.R): g = 0.487, IC 95% [0.469; 0.505]. d_av y Glass no tenían IC.\n")
cat("  OJO: el IC analítico de g supone varianzas iguales. Acá el grupo grande (CdF) es también el más variable,\n")
cat("  así que esa fórmula sobreestima el error estándar (compará EE bootstrap vs. ee_g de 04_efecto.R).\n")
cat("  Para el reporte, preferir el IC bootstrap de g y explicar por qué difiere del analítico.\n")
cat("  Convención: g es la medida principal (la prevista en el plan); d_av y Glass son contraste por la heterocedasticidad.\n")

# --- PASO 5: síntesis ---
razon <- tabla %>%
  select(esquema, medida, ancho, excluye_cero) %>%
  pivot_wider(names_from = esquema, values_from = c(ancho, excluye_cero)) %>%
  mutate(razon_ancho = ancho_secuencia / ancho_ingenuo,
         misma_lectura = excluye_cero_ingenuo == excluye_cero_secuencia)
cat("\n[Paso 5] Síntesis: ¿el IC ingenuo y el por secuencia llevan a la misma conclusión?\n")
cat("  razon_ancho = ancho del IC por secuencia / ancho del IC ingenuo (>1: respetar las réplicas ensancha el IC)\n")
print(as.data.frame(razon %>% select(medida, ancho_ingenuo, ancho_secuencia, razon_ancho, misma_lectura)), digits = 4)
write_csv(tabla, "reportes/05_bootstrap_tabla.csv")
cat("[Paso 5] Tabla guardada en reportes/05_bootstrap_tabla.csv\n")

# --- PASO 6 ---
cat("\n[Paso 6] Semilla fijada: set.seed(2026), al inicio de este script.\n")
cat("[Paso 6] Limitación declarada: la dependencia entre sismos (réplicas) se maneja vía remuestreo/permutación por secuencia, no vía declustering real (fuera del alcance del trabajo). Las secuencias son un proxy (ver comentario del Paso 1). Ver paso 11 de la limpieza.\n")