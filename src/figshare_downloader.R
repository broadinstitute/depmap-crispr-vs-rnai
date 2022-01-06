
data_dir <- file.path("data","raw")
zipfile <- file.path(data_dir,"collection.zip")

options(timeout=200)
download.file(url ="https://figshare.com/ndownloader/articles/16735132?private_link=24131e4e0b894826036e",
              destfile=zipfile)

unzip(zipfile,exdir=data_dir)
