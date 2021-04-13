#' Diagnostics plots for log fold changes
#'
#' @description Generates diagnostic plots of output from tradis_comparison.R.
#'
#' @param path Path to file with *.csv files
#'
#' @export
diagnostics_fc <- function(path){
  myfiles <- lapply(list.files(path = path, pattern = "*.csv", full.names = TRUE), read.delim)
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
  logfc <- replace2 %>% select(-contains(c("_pvalue")))
  logfc$ob <- 1:nrow(logfc)
  meltfc <- reshape2::melt(logfc, id.vars = "ob")
  meltfc$cond <- gsub(pattern = "_logFC", replacement = "", meltfc$variable)

  fc <- ggplot2::ggplot(meltfc, aes(x = ob, y = value)) +
    geom_point(size = 0.1) +
    facet_wrap(~cond) +
    ylab(label = "Log2 Fold Change") +
    xlab(label = "Locus") +
    scale_x_continuous(breaks = c(0, (round(nrow(replace2), digits = -3))/2, round(nrow(replace2), digits = -3))) +
    theme_bw()

  pvalue <- replace2 %>% select(-contains(c("logFC")))
  pvalue$ob <- 1:nrow(pvalue)
  pmelt <- reshape2::melt(pvalue, id.vars = "ob")
  pmelt$cond <- gsub(pattern = "_pvalue", replacement = "", pmelt$variable)

  pv <- ggplot2::ggplot(pmelt, aes(x = value)) +
    geom_histogram(binwidth = 0.05, fill = "gray65") +
    facet_wrap(~cond, scales = "free") +
    ylab(label = "Frequency") +
    xlab(label = "P-Value") +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    theme_bw()

    together <- grid.arrange(fc, pv, ncol = 2)
    plot(together)
    png(paste0(path, "diagnostics_FC_P_all.png"), height = 570, width = 1050, units = "px")
    plot(together)
    dev.off()
    paste0("Your diagnostics file has been saved to ", path, "diagnostics_FC_P_all.png")
}



