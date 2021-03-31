get_sig <- function(x, save_path, sig = 0.05, abslogFC = 0, subset){
  x = logfcs
  subset = "gene"
  save_path = "~/phd/testing/"
  sig = 0.05
  abslogFC = 0
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
