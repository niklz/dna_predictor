# build_shareable_script.R
current_hash <- substr(gert::git_log(max = 1)$commit, 1, 7)

script_text <- readLines("pbix_r_script.R")
script_text <- gsub(r"(model_ver <- "DEV")", sprintf('model_ver <- "%s"', current_hash), script_text)

writeLines(script_text, stringr::str_c("shared/pbix_r_script_", current_hash, ".R"))
