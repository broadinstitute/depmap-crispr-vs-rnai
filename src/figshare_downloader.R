
library(httr)

data_dir <- file.path("data","raw")
zipfile <- file.path(data_dir,"collection.zip")
dir.create(data_dir)

url <- "https://figshare.com/ndownloader/articles/16735132?private_link=24131e4e0b894826036e"
GET(url, progress(), write_disk(zipfile, overwrite=TRUE))

unzip(zipfile,exdir=data_dir)
file.remove(zipfile)
