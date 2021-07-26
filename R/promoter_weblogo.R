#' Promoter sequence weblogo
#'
#' @description Matches the significant genes with their promoter regions and creates a weblogo of consensus regions. This currently doesn't work for rRNA.
#'
#' @param logfcs The output from structure_csv()
#' @param promoters A data frame of locus tags and their associated upstream sequences. Can be generated using Artemis (Select --> All CDS features. File --> Write --> Upstream bases of selected features --> FASTA format). However, will need to edit for locus tags instead of gene name. Need the first column to be named locus_tag.
#' @param save_path Location to save the promoter weblogo. Default is the current working directory.
#' @param save_fasta TRUE or FALSE. If TRUE, will save to same directory as weblogo. (default = FALSE)
#'
#' @export
promoter_weblogo <- function(logfcs, promoters, save_plot = FALSE, save_path, save_fasta = FALSE){
  for (i in seq(4, ncol(x), by=2)){
    wd <- getwd()
    if(missing(save_plot)){save_path = FALSE}
    if(missing(save_path)){save_path = wd}
    data <- logfcs[,c(1:3,i,i+1)]
    data <- subset(data, data[,5]<0.05)
    data2 <- merge(data, promoters, by = "locus_tag", all.y = FALSE, all.x = TRUE)
    data2$ins <- ifelse(data2[,4]>0, "inc", "dec")
    data2 <- data2[order(data2[,4], decreasing = TRUE),]
    data2$newname <- paste0(data2$gene, "_", data2$ins)
    data2 <- data2[!is.na(data2$promoter),]
    name <- gsub(pattern = "_logFC", replacement = "", colnames(data)[4])

    p1 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data2$promoter), 1, 25),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE) +
      theme(axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 40))
    p1$scales$scales[[1]] <- scale_x_continuous(breaks= seq(1,25,by=5),labels=c("-50", "-45", "-40", "-35", "-30"))
    p2 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data2$promoter), 26, 50),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE) +
      labs(x = "Position Relative to TSS") +
      theme(axis.title = element_text(size = 20),
            plot.title = element_text(size = 40))
    p2$scales$scales[[1]] <- scale_x_continuous(breaks= c(seq(1,25,by=5), 25),labels=c("-25", "-20", "-15", "-10", "-5", "-1"))
    combined <- p1 + p2

    final <- combined & patchwork::plot_layout(nrow = 2, guides = "collect") & patchwork::plot_annotation(
      title = "Promoter Sequence Logo 50bp Upstream\nof Significant Genes",
      subtitle = paste0("Condition: ", name),
      theme = theme(plot.title = element_text(size = 40, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5)))

    print(final)

    if (save_plot == TRUE){
    png(file = paste0("weblogo_", name, ".png"), height = 600, width = 1050, units = "px")
    combined & patchwork::plot_layout(nrow = 2, guides = "collect") & patchwork::plot_annotation(
      title = "Promoter Sequence Logo 50bp Upstream\nof Significant Genes",
      theme = theme(plot.title = element_text(size = 40, hjust = 0.5)))
    dev.off()
    }
    if (save_fasta == TRUE){
    write.fasta(as.list(data2$promoter), as.vector(data2$newname), file.out = paste0(save_path, name, ".fa"), as.string = FALSE, )
    }
  }
}
