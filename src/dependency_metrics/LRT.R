
library(purrr)
library(tidyverse)
library(magrittr)
library(MASS) 
library(sn)
library(data.table)
library(readr)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
start_index <- as.numeric(args[2])
end_index <- as.numeric(args[3])
out_file <- args[4]

data <- fread(data_file) %>% column_to_rownames(.,var="Row.name")

if (is.na(start_index)){
  start_index <- 0
} 

if (is.na(end_index)){
  end_index <- ncol(data)
}

suffix <- paste0("-",start_index,"-",end_index)
prefix <- file_path_sans_ext(out_file)
out_file <- paste0(prefix,suffix,".",file_ext(out_file))

data <- data[,(start_index+1):end_index]

is.error <- function(x) inherits(x, "try-error")

LRT_init <- function(vec,starting_nu){
  init_mod <- data.frame(data = vec) %>%
    selm(data ~ 1, family = "ST", fixed.param = list(nu=starting_nu), data = .)
  st_LL <- data.frame(data = vec) %>%
    selm(data ~ 1, family = "ST", data = ., start = c(coef(init_mod, param.type = 'DP'), list(nu = starting_nu))) %>%
    logLik %>%
    as.numeric()
  return(st_LL)
}

LRT_test <- function(vec) {
  min_length <- 10
  stopifnot(is.vector(vec))
  vec <- vec[!is.na(vec)]
  if (length(vec) < min_length) {return(NA)}
  st_LL <- tryCatch({
    data.frame(data = vec) %>%
      selm(data ~ 1, family = "ST", data = .) %>%
      logLik %>%
      as.numeric()
  }, error = function(e) {
    
    print('Default fit failed, trying set nu values')
    nu_range <- c(2, 5, 10, 25, 50, 100, 250, 500, 1000)
    i = 1; unsolved = T
    while (i <= length(nu_range) & unsolved){
      st_LL <- try(LRT_init(vec,nu_range[i]), silent=T)
      unsolved <- is.error(st_LL)
      i = i + 1
    }
    if (unsolved){
      return(NA)
    } else {
      return(st_LL)
    }
  })
  if (!is.na(st_LL)){
    n_LL <- fitdistr(vec, 'normal')$loglik
    return(2*(st_LL - n_LL))
  } else {
    return(NA)
  }
}

res <- map_dbl(data,LRT_test)

df <- data.frame(`Row.name`=names(res),LRT=res)

write_csv(df,path=out_file)
