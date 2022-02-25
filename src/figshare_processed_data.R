
library(httr)

data_dir <- file.path("data","processed")
zipfile <- file.path(data_dir,"collection.zip")
if (!file.exists(data_dir)) {dir.create(data_dir)}

url <- "https://figshare.com/ndownloader/articles/17948639/versions/1"
GET(url, progress(), write_disk(zipfile, overwrite=TRUE))

unzip(zipfile,exdir=data_dir)
file.remove(zipfile)
