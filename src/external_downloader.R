
destination_dir <- file.path("data","raw")

external_files <- list("avana-public-20q2-gene-effect-unscaled.csv"="https://ndownloader.figshare.com/files/22629077")

download.file(url = "https://ndownloader.figshare.com/files/22629077",
              destfile = file.path(destination_dir,"avana-public-20q2-gene-effect-unscaled.csv"))

for (fname in names(external_files)){
  download.file(url =figshare_files[[fname]],
                destfile=file.path(destination_dir,fname))
}
