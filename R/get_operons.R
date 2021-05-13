#' Plot transcription units
#'
#' @description Takes a tab separated file with lists of genes in transcription units (made with BioCyc in mind) and output from structure_csv() to show transcription unit fold changes.
#'
#' @param logfcs Output from structure_csv()
#' @param x Tab separated file of BioCyc transcription unit output.
#' @param abslogfc LogFC cutoff for output. Defaults to 0 to obtain all significant genes.
#' @param sig Significance level cutoff for output. Defaults to 0.05.
#'
#' @export
get_operons <- function(logfcs, x, abslogfc = 0, sig = 0.05){
#logfcs <- read.delim("~/phd/data/TraDIS/analysis/corrected_logfcs.txt", sep = "\t")
#ectable <- read.delim("~/phd/data/TraDIS/analysis/allgenes_operons_pathways.txt", sep = "\t")
ectable <- x[,c(1,3)]
colnames(ectable) <- c("gene", "operon")
operons <- tidyr::separate_rows(ectable, operon, sep = " // ", convert = TRUE)

for (i in seq(4, ncol(logfcs), by = 2)){
  data <- logfcs[,c(1:3, i, i+1)]
  sub <- subset(data, data[,5] < sig)
  sub <- subset(sub, abs(sub[,4]) > abslogfc)

  merged <- merge(sub, operons, by = "gene", all.x = TRUE)
  genes <- as.data.frame(merged[,c(1,6)])
  colnames(genes) <- c("sig_gene", "gene")
  mergeFC <- merge(genes, data, by = "gene", all.x = TRUE)
  mergeFC$sig <- ifelse(mergeFC[,6] >= 0.05, "Insignificant",
                              ifelse(mergeFC[,5] > 0, "Significant +ve", "Significant -ve"))
  paste0("Preparing plot for ", gsub("_logFC", replacement = "", colnames(data)[4]))
  p <- ggplot2::ggplot(mergeFC, aes(x = gene, y = mergeFC[,5], group = gene, fill = sig)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("azure4", brewer.pal(2, "Set2"))) +
    ggplot2::facet_wrap(~sig_gene, scales = "free") +
    ggplot2::theme(axis.text.x = element_text(angle = 25)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::labs(x = "Gene", title = gsub("_logFC", replacement = "", colnames(data)[4]), y = "Log2 Fold Change")
  print(p)
}
}
