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
library(stringr)
suppressMessages(library(dplyr))
library(ggplot2)
library(cowplot)
suppressMessages(library(viridis))
library(forcats)
library(grid)
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))

data <- read_csv(
  file = opt$file,
  na = c("", "NA", "N/A", "#DIV/0!"),
  show_col_types = FALSE,
  name_repair = "universal"
)
data$condition <- as.factor(data$condition)

# t-test per column against condition, assuming equal group variance
res <- data %>% t_test(as.formula(paste(colnames(data[3]), "~ condition")), detailed = TRUE)
for (i in colnames(data[4:ncol(data)])) {
  res <- add_row(res, data %>% t_test(as.formula(paste(i, "~ condition")), detailed = TRUE))
}
res <- res %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
res <- rename(res, term = .y.)

print(res)

write_csv(
  x = res,
  file = paste0(tools::file_path_sans_ext(opt$file), "_multiple_ttests.csv")
)

# Boxplot of each gene
boxplots <- list()
for (gene in colnames(data[, 3:ncol(data)])) {
  boxplots[[gene]] <- local({
    gene <- gene
    stat_test <- data %>% t_test(as.formula(paste(gene, "~ condition")))
    stat_test <- stat_test %>% mutate(y.position = as.numeric(data %>% get_y_position(as.formula(paste(gene, "~ condition"))) %>% select(y.position)))
    stat_test <- stat_test %>% mutate(p.adj = as.character(res %>% filter(term == paste(gene)) %>% select(p.adj)))
    stat_test <- stat_test %>% mutate(p.adj.signif = as.character(res %>% filter(term == paste(gene)) %>% select(p.adj.signif)))
    stat_test <- stat_test %>% rename(condition = .y.)
    stat_test_last <<- stat_test

    p1 <- ggboxplot(data,
      x = paste("condition"),
      y = gene,
      fill = "condition",
      order = levels(data$condition)
    ) +
      stat_pvalue_manual(stat_test,
        label = "p.adj.signif",
        tip.length = 0.01,
      label.size = 6
      ) +
      labs(
        title = str_replace_all(gene, "\\.", " "),
      subtitle = bquote(P[adj]== .(format(signif(as.numeric(res[res$term == gene, 16]), digits = 3), scientific = -2, digits = 2)))
      ) +
      ylab(NULL) +
      xlab(NULL) +
      scale_fill_viridis(discrete = TRUE, option = "viridis", begin = .2, end = .8) +
      theme_pubr() +
      theme(
        legend.position = "none",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = element_text(
          hjust = 0.5
        )
      )
  })
}

p1 <- suppressWarnings(plot_grid(plotlist = boxplots))

title <- ggdraw() +
  draw_label(
    paste0('T-test of each gene in "', basename(opt$file), '"'),
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
