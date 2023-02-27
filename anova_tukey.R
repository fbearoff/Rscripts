#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("-f", "--file"),
    type = "character", default = NULL,
    help = "dataset file name", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n",
    call. = FALSE
  )
}

library(readr)
library(purrr)
suppressMessages(library(viridis))
library(ggplot2)
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(cowplot))
library(grid)
suppressMessages(library(gridExtra))
library(forcats)

data <- read_csv(
  file = opt$file,
  na = c("", "NA", "N/A", "#DIV/0!"),
  show_col_types = FALSE,
  name_repair = "universal"
)
data$condition <- fct_inorder(data$condition)

# 1-way ANOVA with Tukey HSD
fit_aov <- function(col) {
  aov(col ~ condition, data = data)
}

anovas <- map(data[, 3:ncol(data)], fit_aov)

res <- data.frame(
  gene = character(),
  anova_pval = numeric(),
  fdr = numeric()
)

for (cond in rownames(TukeyHSD(anovas[[1]])$condition)) {
  res[cond] <- numeric()
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
    "', 'p adj']"
  )
}
ra_command <- paste0(ra_command, ")")

for (i in names(anovas)) {
  p <- summary(anovas[[i]])[[1]][["Pr(>F)"]][1]
  d <- as.data.frame(TukeyHSD(anovas[[i]])$condition)
  res <- eval(parse(text = ra_command))
}

res$fdr <- p.adjust(
  p = res$anova_pval,
  method = "fdr"
)

write_csv(
  x = res,
  file = paste0(tools::file_path_sans_ext(opt$file), "_anova_tukey.csv")
)

print(res)

# Boxplot of each gene
boxplots <- list()
suppressWarnings(for (gene in colnames(data[, 3:ncol(data)])) {
  boxplots[[gene]] <- local({
    gene <- gene
    res_aov <- data %>% suppressMessages(anova_test(as.formula(paste(
      gene,
      "~ condition"
    ))))
    stat_test <- data %>% tukey_hsd(as.formula(paste(gene, "~ condition")))
    stat_test <- stat_test %>% add_xy_position(x = "condition")
    stat_test_last <<- stat_test
    p1 <- ggboxplot(data,
      x = paste("condition"),
      y = gene,
      fill = "condition",
      order = levels(data$condition),
      shape = 1,
      add = "jitter"
    ) +
      stat_pvalue_manual(stat_test,
        label = "p.adj.signif",
        tip.length = 0.01,
        hide.ns = TRUE
      ) +
      labs(
        title = gene,
        subtitle = bquote(ANOVA ~ p[adj] == .(format(signif(as.numeric(res[res$gene == gene, 3]), digits = 3), scientific = -2, digits = 2)))
      ) +
      ylab(NULL) +
      xlab(NULL) +
      scale_fill_viridis(discrete = TRUE) +
      theme(
        legend.position = "none",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = rel(0.8)
        )
      )
  })
})

p1 <- suppressWarnings(cowplot::plot_grid(plotlist = boxplots))

title <- ggdraw() +
  draw_label(
    paste0('ANOVA of each gene in "', basename(opt$file), '"'),
    fontface = "bold",
    hjust = 0.5
  )

test_stat_caption <- ggdraw() +
  draw_label(get_pwc_label(stat_test_last),
    fontface = "plain",
    hjust = 0.5
  )

stars <- ggdraw() +
  draw_label("* p<0.05, ** p<0.01, *** p<0.001",
    hjust = 0.5
  )

caption <- plot_grid(
  test_stat_caption, stars,
  ncol = 2
)

p1 <- plot_grid(
  title, p1, caption,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

y_title <- textGrob("Abundance",
  gp = gpar(fontface = "bold", fontsize = 15),
  rot = 90
)

pdf(NULL)
p1 <- suppressWarnings(grid.arrange(arrangeGrob(p1, left = y_title)))
invisible(dev.off())

ggsave(
  plot = p1,
  file = paste0(tools::file_path_sans_ext(opt$file), "_boxplots.pdf"),
  width = 22,
  height = 17,
  units = "in",
  dpi = "print"
)
