#' Plot fold changes
#'
#' @description Creates fold change plots for a specified gene list.
#'
#' @param logfcs The output from structure_csv()
#' @param gene_list A character vector of genes for which you want to plot the log2 fold changes for.
#' @param sig Significance cut-off level. Default = 0.05.
#' @param pathway_name Optional. A character of the pathway name, or what you want to prefix your plots (e.g. "Amino acid biosynthesis pathway").
#' @param plot_type Either "gene" or "condition". Chooses which variable to facet by.
#' @param save_plot TRUE or FALSE. If TRUE, will save to your current directory instead of outputting to the plot window. Default is FALSE.
#'
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_manual theme element_text labs
#' @importFrom stats complete.cases
#' @importFrom grDevices dev.off pdf
#' @importFrom utils read.csv2 read.delim
#' @export
#'
logfc_plots <- function(logfcs, gene_list, sig = 0.05, plot_title, plot_type, save_plot = FALSE)
{
  if(missing(plot_title)){plot_title = "Genes of Interest"}
  subset <- logfcs[logfcs$gene %in% gene_list, ]
  subset2 <- as.data.frame(subset)
  rownames(subset2) <- subset2$gene
  subset2 <- subset2[-c(1:3)]

  k <- 2
  nr = nrow(subset2)
  nc <- ncol(subset2)
  unames <- c("logFC", "qvalue")

  a <- array(as.matrix(subset2), c(nr, k, nc/k))
  m <- matrix(aperm(a, c(1,3,2)),, k, dimnames = list(NULL, unames))
  data <- as.data.frame(m, stringsAsFactors = FALSE)

  data$gene <- paste(rownames(subset2))
  char <- gsub("_logFC", "", colnames(logfcs)[grep("*logFC", colnames(logfcs))])
  data$condition <- rep(char, each = length(unique(subset$gene)))
  data$Significance <- ifelse(data$qvalue >= sig, "Insignificant",
                              ifelse(data$logFC > 0, "Significant +ve", "Significant -ve"))
  data <- data[complete.cases(data),]
  data$logFC <- as.numeric(data$logFC)

  if (plot_type == "condition"){
    # plot facet by condition
    plot <- ggplot2::ggplot(data, aes(x = gene, y = logFC, group = gene, fill = Significance)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~condition) +
      ggplot2::scale_fill_manual(values = c("gray20", "green4", "red")) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
            plot.title = element_text(hjust = 0.5, size = 20),
            strip.background = ggplot2::element_rect(colour = "black", fill = "gray90"),
            strip.text = element_text(color = "black", size = 12)) +
      ggplot2::labs(x = "Gene", title = plot_title, y = "Log2 Fold Change") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = paste0("Significance (", sig, ")")))
    print(plot)
  }
  if (plot_type == "gene"){
    plot <- ggplot2::ggplot(data, aes(x = condition, y = logFC, group = condition, fill = Significance)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~gene) +
      ggplot2::scale_fill_manual(values = c("gray20", "green4", "red")) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
            plot.title = element_text(hjust = 0.5, size = 20),
            strip.background = ggplot2::element_rect(colour = "black", fill = "gray90"),
            strip.text = element_text(color = "black", size = 12)) +
      ggplot2::labs(x = "Gene", title = plot_title, y = "Log2 Fold Change") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = paste0("Significance (", sig, ")")))
    print(plot)
  }
  if (save_plot == TRUE){
    print(paste0("Saving to Directory: ", getwd()))
    pdf(file = paste0(plot_title, "_", paste0(plot_type), ".pdf"), width = 10, height = 10)
    print(plot)
    dev.off()
  }
}
