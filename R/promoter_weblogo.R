#' Promoter sequence weblogo
#'
#' @description Matches the significant genes with their promoter regions and creates a weblogo of consensus regions. This currently doesn't work for rRNA.
#'
#' @param x The output from structure_csv()
#' @param promoters A data frame of locus tags and their associated upstream sequences. Can be generated using Artemis (Select --> All CDS features. File --> Write --> Upstream bases of selected features --> FASTA format). However, will need to edit for locus tags instead of gene name. Need the first column to be named locus_tag.
#' @param save_path Location to save the promoter weblogo
#' @param save_fasta TRUE or FALSE. If TRUE, will save to same directory as weblogo. (default = FALSE)
#'
#' @export
promoter_weblogo <- function(x, promoters, save_path, save_fasta = FALSE){
  for (i in seq(4, ncol(x), by=2)){
    wd <- getwd()
    if(missing(save_path)){save_path = wd}
    data <- x[,c(1:3,i,i+1)]
    data <- subset(data, data[,5]<0.05)
    data2 <- merge(data, promoters, by = "locus_tag", all.y = FALSE, all.x = TRUE)
    data2$ins <- ifelse(data2[,4]>0, "inc", "dec")
    data2 <- data2[order(data2[,4], decreasing = TRUE),]
    data2$newname <- paste0(data2$gene, "_", data2$ins)
    data2 <- data2[!is.na(data2$promoter),]
    data_inc <- subset(data2, data2$ins == "inc")
    data_dec <- subset(data2, data2$ins == "dec")
    name <- gsub(pattern = "_logFC", replacement = "", colnames(data)[4])
    p1 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data2$promoter), 1, 25),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE)
    p1$scales$scales[[1]] <- scale_x_continuous(breaks= seq(1,25,by=5),labels=c("-50", "-45", "-40", "-35", "-30"))
    p2 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data2$promoter), 26, 50),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE) +
      labs(x = "Position Relative to TSS")
    p2$scales$scales[[1]] <- scale_x_continuous(breaks= c(seq(1,25,by=5), 25),labels=c("-25", "-20", "-15", "-10", "-5", "-1"))
    p3 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data_inc$promoter), 1, 25),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE)
    p3$scales$scales[[1]] <- scale_x_continuous(breaks= seq(1,25,by=5),labels=c("-50", "-45", "-40", "-35", "-30"))
    p4 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data_inc$promoter), 26, 50),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE) +
      labs(x = "Position Relative to TSS")
    p4$scales$scales[[1]] <- scale_x_continuous(breaks= c(seq(1,25,by=5), 25),labels=c("-25", "-20", "-15", "-10", "-5", "-1"))
    p5 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data_dec$promoter), 1, 25),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE)
    p5$scales$scales[[1]] <- scale_x_continuous(breaks= seq(1,25,by=5),labels=c("-50", "-45", "-40", "-35", "-30"))
    p6 <- ggplot() +
      ggseqlogo::geom_logo(substring(toupper(data_dec$promoter), 26, 50),
                           seq_type = "dna", method = "probability", rev_stack_order = TRUE) +
      labs(x = "Position Relative to TSS")
    p6$scales$scales[[1]] <- scale_x_continuous(breaks= c(seq(1,25,by=5), 25),labels=c("-25", "-20", "-15", "-10", "-5", "-1"))
    gridExtra::grid.arrange(p1,p3,p5,p2,p4,p6, nrow = 2, top = name)
    gridExtra::grid.arrange(p1,p2, nrow = 2, top = name)
    if (save_fasta == TRUE){
    write.fasta(as.list(data2$promoter), as.vector(data2$newname), file.out = paste0(save_path, name, ".fa"), as.string = FALSE, )
    }
  }
}
