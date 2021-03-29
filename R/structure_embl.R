#' Data structuring
#'
#' @description Takes Bio-TraDIS output .csvs and corresponding EMBL file and restructures data for further tradisanalyseR scripts. All .csv files and .embl files will be uploaded.
#'
#' @param csvpath The path to your Bio-TraDIS csv output files.
#' @param emblpath The path to the folder containing the corresponding organism .embl file (optional - leave blank if the embl file is in the same directory as your input files).
#'
#' @return
#' @importFrom dplyr %>% full_join select contains
#' @export
#'
#' @examples
#'
#' @usage x <- structure_embl(csvpath = "~/path/to/csvs/", emblpath = "~/path/to/embl/")
#'
structure_embl <- function(csvpath = "", emblpath = ""){
  wd <- getwd()
  if(missing(emblpath)){emblpath = csvpath}
  setwd(emblpath)
  embl <- read.csv2(list.files(pattern = "*.embl"), sep = " ")
  embl <- as.data.frame(embl)
  embl_filter <- embl[, sapply(embl, function (x) any(grepl('locus_tag', x)))]
  tag <- stringr::str_match(embl_filter, pattern = 'locus_tag.*')
  tag <- tag[!is.na(tag)][1]
  tag <- gsub("locus_tag=", "", tag)
  tag <- gsub("_.*", "", tag)

  locus_tags <- stringr::str_extract(embl_filter, paste0(tag, "_[0-9]+"))
  locus_tags <- as.data.frame(unique(locus_tags[!is.na(locus_tags)]))
  colnames(locus_tags) <- "locus_tag"

  setwd(csvpath)
  myfiles <- lapply(list.files(pattern = "*.csv"), read.delim)

  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  all_locus <- list(joined, locus_tags) %>% purrr::reduce(full_join, by = "locus_tag")

  info <- all_locus[,c(1:3)]
  colnames(info) <- c("locus_tag", "gene", "function")

  filenames <- list.files(pattern = "*.csv") %>% gsub(pattern = ".csv", replacement = "")
  setwd(wd)
  replace <- all_locus[,-c(1:3)]
  replace2 <- replace %>% select(-contains(c("CPM", "PValue", "gene_name", "function")))

  names <- character(0)
  for (i in 1:length(filenames)){
    logfc <- paste0(filenames[i], "_logFC")
    names <- append(names, logfc)
    q <- paste0(filenames[i], "_qvalue")
    names <- append(names, q)
  }

  colnames(replace2) <- names
  out <- cbind(info, replace2)
  return(out)
}
