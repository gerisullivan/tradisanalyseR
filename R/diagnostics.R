#' Diagnostics plots for read counts
#'
#' @description Generates diagnostic plots using read counts as input. File names are expected to follow \*_1.tradis_gene_insert_sites.csv, \*_2.tradis_gene_insert_sites.csv. Conditions will be named based on input names before the replicate number. All conditions and controls are expected to have two replicates.
#'
#' @param path Path to file with *tradis_gene_insert_sites.csv files
#' @param controls A string (whole or partial) to identify the controls (e.g. "control", or "time0")
#'
#' @export
diagnostics <- function(path, controls){
  print("Merging files together")
  myfiles <- lapply(list.files(path = path, pattern = "*sites.csv", full.names = TRUE), read.delim)
  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag") #join together by locus tag
  filenames <- list.files(path = path, pattern = "*sites.csv") %>%
    gsub(pattern = ".tradis_gene_insert_sites.csv", replacement = "") #make file names to rename columns later
  rc <- joined %>% select(contains(c("locus_tag", "read_count"))) #extract only read counts
  #rc <- rc[rowSums(rc[, -1])>0, ]
  locus_tags <- rc$locus_tag # handy for later

  names <- character(0) # rename columns
  for (i in 1:length(filenames)){
    readcount <- paste0(filenames[i])
    names <- append(names, readcount)
    rm(readcount)
  }

  colnames(rc)[2:ncol(rc)] <- names #rename columns
  rownames(rc) <- rc[,1]
  rc <- rc[,-1]
  rc <- rc[rowSums(rc[, -1])>0, ]

  conds <- unique(gsub("_[0-9]", replacement = "", x = filenames))
  conds <- gsub(paste0(".*", controls, ".*"), replacement = "control", conds)
  group <- as.factor(rep(c(conds), each = 2))
  #group <- relevel(group, ref="control")

  col <- as.factor(ifelse(group == "control", "control", "condition"))
  set <- EDASeq::newSeqExpressionSet(as.matrix(rc))

  print("Running comparisons")
  design <- model.matrix(~0+group) # create design
  y <- edgeR::DGEList(counts = as.matrix(rc), group=group) # design DGElist for input
  y <- edgeR::estimateGLMCommonDisp(y, design) # calculate negative binomial dispersion parameter
  y <- edgeR::estimateGLMTagwiseDisp(y, design) # calculate Bayes estimate of negative binomial dispersion parameters
  names <- levels(y$samples$group)
  names <- names[-1]

  conts <- character()
  for (i in 1:length(conds)){
    x <- paste0("group",conds[i], " - groupcontrol")
    conts <- append(conts, x)
  }

  conts <- conts[!conts %in% "groupcontrol - groupcontrol"]

  fit <- edgeR::glmFit(y, design, robust=TRUE)

  tags <- list()
  for (i in 1:length(conts)){
    contrast <- limma::makeContrasts(contrasts = conts[i], levels = design)
    lrt <- edgeR::glmLRT(fit, contrast = contrast)
    tags[[i]] <- lrt$table
  }

  # pvalues <- as.data.frame(rownames(rc))
  # for (i in 1:length(conts)){ #get all p values out
  #   pvalues[,i+1] <- tags[[i]]$PValue
  # }
  #
  # colnames(pvalues) <- c("locus_tag", names) #rename
  # melt_p <- reshape2::melt(pvalues, id.vars = "locus_tag") #melt for formatting for ggplot
  # ggplot(melt_p, aes(x = value)) + #plot p values
  #   geom_histogram(binwidth = 0.05) +
  #   facet_wrap(~variable, scales = "free")

  logfc <- as.data.frame(1:length(rownames(rc)))
  for (i in 1:length(conts)){ # get logFC for all conditions
    logfc[,i+1] <- tags[[i]]$logFC
  }

  colnames(logfc) <- c("observation", names) #rename
  melt_fc <- reshape2::melt(logfc, id.vars = "observation") # melt for formatting for ggplot

  png(paste0(path, "diagnostics.png"), height = 800, width = 1300, units = "px")
  par(mfrow=c(1,2), mar = c(5,5,5,2)) #(b,l,t,r)
  EDASeq::plotPCA(set, col=rep(c("gray60",2)[col]), cex = 1.5)
  plot.new()
  vps <- gridBase::baseViewports()
  grid::pushViewport(vps$figure)
  vp1 <- grid::plotViewport(c(1.3,1,4.5,1.5))
  p2 <- ggplot(melt_fc, aes(x = observation, y = value)) + #plot logFC for all conditions
    geom_point(size = 0.1) +
    facet_wrap(~variable) +
    labs(x = "Locus", y = "Log2 Fold Change") +
    theme(strip.background = NULL, strip.text = element_text(size = 12),
          axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
    scale_x_continuous(breaks = c(0, (round(nrow(rc), digits = -3))/2, round(nrow(rc), digits = -3)))
  print(p2, vp = vp1)
  mtext(expression(paste(bold("Diagnostics"), " - Principle Component Analysis and Locus vs Log" [2], " Fold Change")), side = 3, line = -4, outer = TRUE, cex = 2.2)
  dev.off()

  paste0("Your diagnostics file has been saved to ", path, "diagnostics.png")

}
