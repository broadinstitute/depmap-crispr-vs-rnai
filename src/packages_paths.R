library(tidyverse)
library(magrittr)
library(data.table)

data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")

source(file.path("src","id_utility.R"))
source(file.path("src","rank_utility.R"))
source(file.path("src","BDP_helper_functions.R"))
source(file.path("src","load_data_utility.R"))