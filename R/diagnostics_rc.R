#' Diagnostics plots for read counts
#'
#' @description Generates diagnostic plots using read counts as input.
#'
#' @param path Path to file with *tradis_gene_insert_sites.csv files
#'
#' @export
diagnostics_rc <- function(path){
  myfiles <- lapply(list.files(path = path, pattern = "*sites.csv", full.names = TRUE), read.delim)
  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag") #join together by locus tag
  filenames <- list.files(path = path, pattern = "*sites.csv") %>%
    gsub(pattern = ".tradis_gene_insert_sites.csv", replacement = "") #make file names to rename columns later
  rc <- joined %>% select(contains(c("locus_tag", "read_count"))) #extract only read counts
  #rc <- rc[rowSums(rc[, -1])>0, ]
  locus_tags <- rc$locus_tag # handy for later

  names <- character(0) # rename columns
  for (i in 1:length(filenames)){
    readcount <- paste0(filenames[i])
    names <- append(names, readcount)
    rm(readcount)
  }

  colnames(rc)[2:ncol(rc)] <- names #rename columns
  rownames(rc) <- rc[,1]
  rc <- rc[,-1]
  rc <- rc[rowSums(rc[, -1])>0, ]

  set <- EDASeq::newSeqExpressionSet(as.matrix(rc))
  EDASeq::plotPCA(set)

  png(paste0(path, "PCAplot_all.png"), height = 800, width = 1300, units = "px")
  EDASeq::plotPCA(set)
  dev.off()
}
