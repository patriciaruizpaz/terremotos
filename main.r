# =====================================================================
# main.R — corre el pipeline completo del proyecto, en orden.
# Correr este archivo entero (Ctrl+Shift+Enter) reproduce todo el
# análisis desde cero a partir de datos/crudo.csv.
# =====================================================================

# Paquetes usados en todo el proyecto: instalar una sola vez si falta
# alguno.
# install.packages(c("tidyverse", "gridExtra"))

source("scripts/01_limpieza.R")
source("scripts/02_decisiones.R")
source("scripts/03_ic.R")
source("scripts/04_efecto.R")
# source("scripts/05_bootstrap.R")