#' Data structuring
#'
#' @description Takes Bio-TraDIS output .csvs and restructures data for further tradisanalyseR scripts. All .csv files in the folder will be obtained.
#'
#' @param csvpath The path to your Bio-TraDIS csv output files.
#'
#' @return
#' @importFrom dplyr %>% full_join select contains
#' @export
#'
structure_csv <- function(csvpath = ""){
  myfiles <- lapply(list.files(path = csvpath, pattern = "*.csv", full.names = TRUE), read.delim, sep = ",")

  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  info <- joined[,c(1:3)]
  colnames(info) <- c("locus_tag", "gene", "function")

  filenames <- list.files(path = csvpath, pattern = "*.csv") %>% gsub(pattern = ".csv", replacement = "")
  replace <- joined[,-c(1:3)]
  replace <- replace %>% select(-contains(c("CPM", "PValue", "gene_name", "function")))

  names <- character(0)
  for (i in 1:length(filenames)){
    logfc <- paste0(filenames[i], "_logFC")
    names <- append(names, logfc)
    q <- paste0(filenames[i], "_qvalue")
    names <- append(names, q)
  }

  colnames(replace) <- names
  out <- cbind(info, replace)
  return(out)
}
