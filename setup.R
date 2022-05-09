if(!require("BiocManager")){
  install.packages("BiocManager")
  library("BiocManager")
}

if(!require("Rsubread")){
  BiocManager::install("Rsubread")
  library("Rsubread")
}

if(!require("shiny")){
  install.packages("shiny")
  library("shiny")
}

if(!require("shinyWidgets")){
  install.packages("shinyWidgets")
  library("shinyWidgets")
}

if(!require("bslib")){
  install.packages("bslib")
  library("bslib")
}

if(!require("reactable")){
  install.packages("reactable")
  library("reactable")
}

if(!require("dplyr")){
  install.packages("dplyr")
  library("dplyr")
}

if(!require("purrr")){
  install.packages("purrr")
  library("purrr")
}

if(!require("knitr")){
  install.packages("knitr")
  library("knitr")
}

if(!require("stringi")){
  install.packages("stringi")
  library("stringi")
}

if(!require("stringdist")){
  install.packages("stringdist")
  library("stringdist")
}

source("functions/get_quadruplotype.R")
source("functions/read_sam.R")
source("functions/get_common_string.R")
source("ui.R")
source("server.R")

options(shiny.maxRequestSize = 1*1024^3)