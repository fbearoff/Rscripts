# Rscripts

Scripts I use for semi-automated data analysis in R.

[ANOVA per column with Dunnett test](anova_dunnett.R) 
- Supply CSV with samples as rows, condition as 1st column, and genes as columns.
- Indicate control condition with "-c" flag.

[ANOVA per column with Tukey HSD test](anova_tukey.R) 
- Supply CSV with samples as rows, sample id as first column, condition as second column, and genes/analytes as subsequent columns. Conditions are plotted in order of first appearance.

[Generate tx2gene](tx2gene.R) 
- Stub to generate a tx2gene object as a table
