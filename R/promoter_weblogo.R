#' Promoter sequence weblogo
#'
#' @description Matches the significant genes with their promoter regions and creates a weblogo of consensus regions. This currently doesn't work for rRNA.
#'
#' @param x The output from structure_embl()
#' @param promoters A data frame of locus tags and their associated upstream sequences. Can be generated using Artemis (Select --> All CDS features. File --> Write --> Upstream bases of selected features --> FASTA format). However, will need to edit for locus tags instead of gene name
#' @param save_path Location to save the promoter weblogo
#' @param save_fasta TRUE or FALSE. If TRUE, will save to same directory as weblogo. (default = FALSE)
#'
#' @export
promoter_weblogo <- function(x, promoters, save_path, save_fasta = FALSE){
  for (i in seq(4, ncol(x), by=2)){
    wd <- getwd()
    data <- x[,c(1:3,i,i+1)]
    data <- subset(data, data[,5]<0.05)
    data2 <- merge(data, promoters, by = "locus_tag", all.y = FALSE, all.x = TRUE)
    data2$ins <- ifelse(data2[,4]>0, "inc", "dec")
    data2 <- data2[order(data2[,4], decreasing = TRUE),]
    data2$newname <- paste0(data2$gene, "_", data2$ins)
    data2 <- data2[!is.na(data2$promoter),]
    name <- gsub(pattern = "_logFC", replacement = "", colnames(data)[4])
    weblogo(as.character(data2$promoter), file.out = paste0(save_path, name, "_weblogo.png"), open = FALSE, format = "png",
            title = name, annotate = c(-100,"","","","",-95,"","","","",-90,"","","","",-85,"","","","",-80,
                                       "","","","",-75,"","","","",-70,"","","","",-65,"","","","",-60,"","",
                                       "","",-55,"","","","",-50,"","","","",-45,"","","","",-40,"","","",
                                       "",-35,"","","","",-30,"","","","",-25,"","","","",-20,"","","","",-15,
                                       "","","","",-10,"","","","",-5,"","","",-1), xlabel = "position relative to TSS")
    if (save_fasta == TRUE){
    write.fasta(as.list(data2$promoter), as.vector(data2$newname), file.out = paste0(save_path, name, ".fa"), as.string = FALSE, )
    }
  }
}
