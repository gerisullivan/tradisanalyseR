#' Plot transcription units
#'
#' @description Takes a tab separated file with lists of genes in transcription units (made with BioCyc in mind) and output from structure_csv() to show transcription unit fold changes.
#'
#' @param logfcs output from structure_csv()
#' @param x tab separated file of BioCyc transcription unit output.
#' @param abslogfc absolute log fold change cutoff for output. Defaults to 0.
#' @param sig significance level cutoff for output. Defaults to 0.05.
#' @param n maximum number of facets to output in one graph. If there are too many significant genes to plot, will show the top n plots based on most significant. Defaults to 36. If you have too many significant genes, decrease your significance cutoff.
#'
#' @export
get_operons <- function(logfcs, x, abslogfc = 0, sig = 0.05, n = 36){
  ectable <- x[,c(1,3)]
  colnames(ectable) <- c("gene", "operon")
  operons <- tidyr::separate_rows(ectable, operon, sep = " // ", convert = TRUE)

  for (i in seq(4, ncol(logfcs), by = 2)){
  data <- logfcs[,c(1:3, i, i+1)]
  sub <- subset(data, data[,5] < sig)
  sub <- subset(sub, abs(sub[,4]) > abslogfc)
  sort <- sub[,c(2,5)]
  sort <- sort[order(sort[,2]),]
  siggenes <- sort$gene

  merged <- merge(sub, operons, by = "gene", all.x = TRUE)
  genes <- as.data.frame(merged[,c(1,6)])
  colnames(genes) <- c("sig_gene", "gene")
  mergeFC <- merge(genes, data, by = "gene", all.x = TRUE)
  mergeFC$sig <- ifelse(mergeFC[,6] >= sig, "Insignificant",
                              ifelse(mergeFC[,5] > 0, "Significant +ve", "Significant -ve"))
  mergeFC <- mergeFC %>% dplyr::group_by(sig_gene) %>% dplyr::filter(n()>1)
  mergeFC <- as.data.frame(mergeFC[!is.na(mergeFC[,5]),])
  mergeFC <- subset(mergeFC, !(as.character(mergeFC$gene) == as.character(mergeFC$sig_gene) & mergeFC$sig == "Insignificant"))
  mergeFC <- mergeFC[!Biobase::isUnique(mergeFC$sig_gene),]

  genes <- intersect(siggenes, mergeFC$sig_gene)[1:n]
  input <- mergeFC[mergeFC$sig_gene %in% genes, ]
  input$sig_gene <- factor(input$sig_gene, levels = genes)

  print(paste0("Preparing plot for ", gsub("_logFC", replacement = "", colnames(data)[4])))
  p <- ggplot2::ggplot(input, aes(x = gene, y = input[,5], group = gene, fill = sig)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("gray20", "green4", "red")) +
    ggplot2::facet_wrap(~sig_gene, scales = "free") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.background = ggplot2::element_rect(colour = "black", fill = "gray90"),
          strip.text = element_text(color = "black", size = 12)) +
    ggplot2::labs(x = "Gene", title = gsub("_logFC", replacement = "", colnames(data)[4]), y = "Log2 Fold Change") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = paste0("Significance (", sig, ")"))) +
    ggplot2::geom_abline(intercept = 0, slope = 0, size = 0.5)
  print(p)
  }
}
