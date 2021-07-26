#' Diagnostics plots for log fold changes
#'
#' @description Generates diagnostic plots of output from tradis_comparison.R.
#'
#' @param path Path to file with *.csv files
#'
#' @export
diagnostics_fc <- function(path){

  path = "~/Desktop/logFC/"

  myfiles <- lapply(list.files(path = path, pattern = "*.csv", full.names = TRUE), read.delim, sep = ",")
  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  filenames <- list.files(path = path, pattern = "*.csv") %>%
    gsub(pattern = ".csv", replacement = "")
  info <- joined[,c(1:3)]
  colnames(info) <- c("locus_tag", "gene", "function")
  replace <- joined[,-c(1:3)]
  replace2 <- replace %>% select(-contains(c("CPM", "q.value", "gene_name", "function")))

  names <- character(0)
  for (i in 1:length(filenames)){
    logfc <- paste0(filenames[i], "_logFC")
    names <- append(names, logfc)
    p <- paste0(filenames[i], "_pvalue")
    names <- append(names, p)
  }

  colnames(replace2) <- names
  replace2 <- replace2[,c(1:(2*length(filenames)))]
  logfc <- replace2 %>% select(-contains(c("_pvalue")))
  logfc$ob <- 1:nrow(logfc)
  meltfc <- reshape2::melt(logfc, id.vars = "ob")
  meltfc$cond <- gsub(pattern = "_logFC", replacement = "", meltfc$variable)

  fc <- ggplot2::ggplot(meltfc, aes(x = ob, y = value)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::facet_wrap(~cond) +
    ggplot2::ylab(label = "Log2 Fold Change") +
    ggplot2::xlab(label = "Locus") +
    ggplot2::scale_x_continuous(breaks = c(0, (round(nrow(replace2), digits = -3))/2, round(nrow(replace2), digits = -3))) +
    ggplot2::theme_bw()

  pvalue <- replace2 %>% select(-contains(c("logFC")))
  pvalue$ob <- 1:nrow(pvalue)
  pmelt <- reshape2::melt(pvalue, id.vars = "ob")
  pmelt$cond <- gsub(pattern = "_pvalue", replacement = "", pmelt$variable)

  pv <- ggplot2::ggplot(pmelt, aes(x = value)) +
    ggplot2::geom_histogram(binwidth = 0.05, fill = "gray65") +
    ggplot2::facet_wrap(~cond, scales = "free") +
    ggplot2::ylab(label = "Frequency") +
    ggplot2::xlab(label = "P-Value") +
    ggplot2::scale_x_continuous(breaks = c(0,0.5,1)) +
    ggplot2::theme_bw()

    together <- gridExtra::grid.arrange(fc, pv, ncol = 2)
    plot(together)
    png(paste0(path, "diagnostics_FC_P_all.png"), height = 570, width = 1050, units = "px")
    plot(together)
    dev.off()
    paste0("Your diagnostics file has been saved to ", path, "diagnostics_FC_P_all.png")
}



