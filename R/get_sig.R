#' Extract significant genes
#'
#' @description Takes the logFC and significance values as outputted by tradis_comparison.R and makes new files. Will output only gene name and logFC (perfect for input into EcoCyc/BioCyc SmartTables)
#'
#' @param x output from structure_embl()
#' @param save_path Path to directory to save output files. If not provided, will save in current working directory.
#' @param sig Significant level for cutoff (default 0.05)
#' @param abslogFC Cutoff level for the absolute value of the logFC (default 0)
#' @param subset Either "gene" or "all". "Gene" will output gene and logFC information only - great for input into EcoCyc. "all" will output locus tags, gene names and function with logFC data. (default= "all)
#'
get_sig <- function(x, save_path, sig = 0.05, abslogFC = 0, subset){
  if(missing(subset)){subset = "all"}
  if (abslogFC < 0) stop("'abslogFC' must be greater than 0")
  if (sig > 1) stop("'sig' must be a proportion between 0 and 1")
  if (sig < 0) stop("'sig' must be a proportion between 0 and 1")
  for (i in seq(4, ncol(x), by = 2)){
    data <- x[,c(1:3, i, i+1)]
    sub <- subset(data, data[,5] < sig)
    if (abslogFC > 0){
      sub <- subset(sub, abs(sub[,4]) > abslogFC)
    }
    if (subset == "all"){
      write.table(sub, file = paste0(save_path, gsub(pattern = "_logFC", replacement = ".tsv", colnames(data)[4])),
                  sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE)
      print(paste0("Writing all data to file ", save_path, gsub(pattern = "_logFC", replacement = ".tsv", colnames(data)[4])))
    }
    if (subset == "gene"){
      write.table(sub[,c(2,4)], file = paste0(save_path, gsub(pattern = "_logFC", replacement = "_gene_FC.tsv", colnames(data)[4])),
                  sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE)
      print(paste0("Writing gene name and logFC info to file ", save_path, gsub(pattern = "_logFC", replacement = "_gene_FC.tsv", colnames(data)[4])))
    }
  }
}
