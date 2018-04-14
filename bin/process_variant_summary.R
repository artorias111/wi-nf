# This script takes parsed variant summary json files and converts them to TSVs
library(tidyverse)

flist <- list.files(pattern="*summary.json")

for(f in flist) {
    df <- jsonlite::read_json(f) %>%
    dplyr::bind_rows(.) %>%
    replace(is.na(.), 0)
    df <- df[,names(df) != ""]
    df %>%
    readr::write_tsv(gsub("json", "tsv", f))
}
