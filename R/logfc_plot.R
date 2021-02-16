logfc_plots <- function(x, gene_list, pathway_name, plot_type, save_plot)
{
  if(missing(save_plot)){save_plot = FALSE}
  subset <- x[x$gene %in% gene_list, ]
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
  char <- gsub("_logFC", "", colnames(x)[grep("*logFC", colnames(x))])
  data$condition <- rep(char, each = length(unique(subset$gene)))
  data$Significance <- ifelse(data$qvalue >= 0.05, "Insignificant",
                              ifelse(data$logFC > 0, "Significant +ve", "Significant -ve"))
  data <- data[complete.cases(data),]
  data$logFC <- as.numeric(data$logFC)

  if (plot_type == "condition"){
    # plot facet by condition
    plot <- ggplot(data, aes(x = gene, y = logFC, group = gene, fill = Significance)) +
      geom_bar(stat = "identity") +
      facet_wrap(~condition) +
      scale_fill_manual(values = c("gray20", "green4", "red")) +
      theme(axis.text.x = element_text(angle = 45)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(x = "Gene", title = pathway_name, y = "Log2 Fold Change")
    plot
  }
  if (plot_type == "gene"){
    plot <- ggplot(data, aes(x = condition, y = logFC, group = condition, fill = Significance)) +
      geom_bar(stat = "identity") +
      facet_wrap(~gene) +
      scale_fill_manual(values = c("gray20", "green4", "red")) +
      theme(axis.text.x = element_text(angle = 45)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(x = "Gene", title = pathway_name, y = "Log2 Fold Change")
    plot
  }
  if (save_plot == TRUE){
    print(paste0("Saving to Directory: ", getwd()))
    pdf(file = paste0(pathway_name, "_", paste0(plot_type), ".pdf"), width = 10, height = 10)
    print(plot)
    dev.off()
  }
}
