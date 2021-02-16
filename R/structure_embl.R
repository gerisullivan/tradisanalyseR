#' Data structuring
#'
#' @description Takes Bio-TraDIS output .csvs and corresponding EMBL file and restructures data for further tradisanalyseR scripts. Working directory must be set to folder containing the correct files. All .csv files and .embl files will be uploaded.
#'
#' @return
#' @importFrom dplyr %>% full_join
#' @export
#'
#' @examples
#'
structure_embl <- function(){
  print(paste0("Your working directory is set to: ", getwd()))
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

  myfiles <- lapply(list.files(pattern = "*.csv"), read.delim)

  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  all_locus <- list(joined, locus_tags) %>% purrr::reduce(full_join, by = "locus_tag")

  info <- all_locus[,c(1:3)]
  colnames(info) <- c("locus_tag", "gene", "function")

  filenames <- list.files(pattern = "*.csv") %>% gsub(pattern = ".csv", replacement = "")
  replace <- all_locus[,-c(1:3)]
  replace2 <- replace %>% select(-contains(c("CPM", "PValue", "gene_name", "function")))

  names <- character(0)
  for (i in 1:length(filenames)){
    logfc <- paste0(filenames[i], "_logFC")
    names <- append(names, logfc)
    q <- paste0(filenames[i], "_qvalue")
    names <- append(names, logfc)
  }

  colnames(replace2) <- names

  x <- cbind(info, replace2)

  assign(x = "logFCs_all", value = x)
}
