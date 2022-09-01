#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("-f", "--file"),
    type = "character", default = NULL,
    help = "dataset file name", metavar = "character"
  ),
  make_option(c("-c", "--control"),
    type = "character", default = NULL,
    help = "control condition", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input file and control condition).n",
    call. = FALSE
  )
}

library(readr)
library(purrr)
library(DescTools)
suppressMessages(library(dplyr))
library(ggplot2)
library(cowplot)
suppressMessages(library(viridis))
library(forcats)
library(grid)
suppressMessages(library(gridExtra))

control <- opt$control

data <- read_csv(
  file = opt$file,
  na = "#DIV/0!",
  show_col_types = FALSE
)
data$condition <- as.factor(data$condition)

# 1-way ANOVA with Dunnett test against control
fit_aov <- function(col) {
  aov(col ~ condition, data = data)
}

anovas <- map(data[, 2:ncol(data)], fit_aov)

res <- data.frame(
  gene = character(),
  anova_pval = numeric(),
  fdr = numeric()
)

for (cond in which(!levels(data$condition) %in% control)) {
  res[[paste(levels(data$condition)[cond], control, sep = "-")]] <- numeric()
}

ra_command <- paste0("add_row(res, gene = i, anova_pval = p")

for (cond in names(res[, 4:ncol(res)])) {
  ra_command <- paste0(
    ra_command,
    ", '",
    cond,
    "' = ",
    "d['",
    cond,
    "', 'pval']"
  )
}
ra_command <- paste0(ra_command, ")")

for (i in names(anovas)) {
  p <- summary(anovas[[i]])[[1]][["Pr(>F)"]][1]
  d <- as.data.frame(DunnettTest(data[[i]], as.factor(data$condition),
    control = control
  )[[1]])
  res <- eval(parse(text = ra_command))
}

res$fdr <- p.adjust(
  p = res$anova_pval,
  method = "fdr"
)

print(res)

write_csv(
  x = res,
  file = paste0(tools::file_path_sans_ext(opt$file), "_anova_dunnett.csv")
)

# Boxplot of each gene
boxplots <- list()
for (gene in colnames(data[, 2:ncol(data)])) {
  boxplots[[gene]] <- local({
    gene <- gene
    p1 <- ggplot(data, aes(
      x = fct_inorder(condition),
      y = data[[gene]],
      fill = condition
    )) +
      geom_boxplot() +
      labs(
        title = gene,
        subtitle = paste0("padj=", signif(res[res$gene == gene, 3], digits = 3))
      ) +
      ylab(NULL) +
      xlab(NULL) +
      scale_y_log10() +
      scale_fill_viridis(discrete = TRUE) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
  })
}

p1 <- suppressWarnings(cowplot::plot_grid(plotlist = boxplots))

title <- ggdraw() +
  draw_label(
    paste0('ANOVA of each gene in "', basename(opt$file), '"'),
    fontface = "bold",
    hjust = 0.5
  )

p1 <- plot_grid(
  title, p1,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

y_title <- textGrob(paste("Expression relative to", control, sep = " "),
  gp = gpar(fontface = "bold", fontsize = 15),
  rot = 90
)

pdf(NULL)
p1 <- suppressWarnings(grid.arrange(arrangeGrob(p1, left = y_title)))
invisible(dev.off())

ggsave(
  plot = p1,
  file = paste0(tools::file_path_sans_ext(opt$file), "_boxplots.pdf"),
  width = 11,
  height = 8.5,
  units = "in",
  dpi = "print"
)
