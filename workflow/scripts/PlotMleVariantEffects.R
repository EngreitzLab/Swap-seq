# Jesse Engreitz
# 5/30/21
# 09/03/26 SC: Rscript to plot MLE effect size estimates


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--mleEffects", type="character", help="MLE effect sizes"),
  make_option("--effectColumn", type="character", default="effect_size"),
  make_option("--outfile", type="character", help="Output plot filename")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


mle <- read.delim(opt$mleEffects, check.names=F, stringsAsFactors=F)
mle <- subset(mle, VariantID != "")  ## Only plot results for desired variants
mle <- subset(mle, sum1 >= 1000) ## Only plot the variants with at least 1000 cells

stopifnot(opt$effectColumn %in% colnames(mle))

## x axis = pegRNA index (trailing integer of the VariantID); reference and controls have none
pegIndex <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", mle$VariantID)))
fallback <- is.na(pegIndex)
mle$Index <- ifelse(fallback,
                    ifelse(grepl("InferredReference", mle$VariantID), "ref",
                           sub("^.*Barcode-", "bc", mle$VariantID)),
                    as.character(pegIndex))
mle <- mle[order(is.na(pegIndex), pegIndex, mle$Index), ] # unindexed last
mle$Index <- factor(mle$Index, levels = unique(mle$Index))

nBars <- nrow(mle)
plotWidth <- max(7, min(48, nBars * 0.055)) # ~0.055 in per bar, clamped
tickSize <- if (nBars > 400) 1.6 else if (nBars > 150) 2.6 else if (nBars > 60) 4 else 6

baseTheme <- theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = tickSize),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5))

subtitle <- sprintf("%s  (%d barcodes with sum1 >= 1000)", basename(opt$mleEffects), nBars)

p1 <- ggplot(mle, aes_string(x = "Index", y = opt$effectColumn)) +
     geom_col(fill = "#4C4C4C", width = 0.8) +
     xlab("pegRNA index") +
     ylab("Effect size on gene expression\n(MLE, vs reference allele)") +
     ggtitle(subtitle) +
     geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
     baseTheme

p2 <- mle %>% mutate(PctEffect = (get(opt$effectColumn) - 1) * 100) %>%
     ggplot(aes_string(x = "Index", y = "PctEffect")) +
     geom_col(fill = "#4C4C4C", width = 0.8) +
     xlab("pegRNA index") +
     ylab("Effect size on gene expression\n(MLE, % change vs reference allele)") +
     ggtitle(subtitle) +
     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
     baseTheme

pdf(file = opt$outfile, width = plotWidth, height = 4.5)
print(p1)
print(p2)
invisible(dev.off())
