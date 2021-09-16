#' Get info from CSVs
#'
#' @description Gives the number of significant genes at log fold change thresholds for all conditions in a file path.
#'
#' @param path path to folder containing csv files
#' @param abslogfc absolute log fold change cutoff for output. Defaults to 0.
#' @param sig significance level cutoff for output. Defaults to 0.05.
#'
#' @export
sig_info <- function(path, abslogfc = 1, sig = 0.05){
  myfiles <- lapply(list.files(path = path, pattern = "*.csv", full.names = TRUE), read.csv)
  filenames <- list.files(path = path, pattern = "*.csv") %>%
    gsub(pattern = ".csv", replacement = "")
stats <- data.frame(condition = character(),
                    numsig = numeric(),
                    numfc = numeric())
  for (i in 1:length(filenames)){
  data <- myfiles[[i]]
  dsig <- subset(data, data$q.value < sig)
  nsig <- as.numeric(nrow(dsig))
  nfc <- as.numeric(nrow(subset(dsig, abs(dsig$logFC) > abslogfc)))
  stats[i,] <- c(filenames[i], nsig, nfc)
  }
print(stats)
}
