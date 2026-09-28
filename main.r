# =====================================================================
# main.R — corre el pipeline completo del proyecto, en orden.
# Correr este archivo entero (Ctrl+Shift+Enter) reproduce todo el
# análisis desde cero a partir de datos/crudo.csv.
#
# Antes de correr: abrir el proyecto en la carpeta raíz (la que contiene
# datos/, scripts/ y reportes/), porque los scripts usan rutas relativas.
# =====================================================================

# Paquetes usados en todo el proyecto: instalar una sola vez si falta
# alguno.
#   tidyverse  -> todos los scripts
#   gridExtra  -> 01_limpieza.R (gráficos del paso 10)
#   DescTools  -> 03_ic.R y 04_efecto.R (IC y tamaños de efecto)
# 05_bootstrap.R usa solo tidyverse (ya no necesita boot ni effectsize).
# install.packages(c("tidyverse", "gridExtra", "DescTools"))

source("scripts/01_limpieza.R")
source("scripts/02_decisiones.R")
source("scripts/03_ic.R")
source("scripts/04_efecto.R")

# 05_bootstrap.R: la primera vez tarda unos 3 minutos (bootstrap, permutaciones
# e inversión) y guarda una caché en reportes/*_v2.rds; las siguientes
# corridas la cargan y reportan en segundos. Para recalcular todo, poner
# RECALCULAR <- TRUE dentro de ese script (o borrar el .rds).
source("scripts/05_bootstrap.R")