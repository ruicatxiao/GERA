#ABX RESISTANCE TABLE

library(gt)

abx_table_df <- read.csv("abx_table_r_category.csv")
abx_table <- abx_table_df |>
  gt(rowname_col = "Site")
abx_table <- abx_table |>
  tab_stubhead(label = "Site")
abx_table <- abx_table |>
  tab_spanner(label = "Resistant and Intermediate Isolates, n (%)",
              columns = c(AMC30, AMP10, FOX30, CRO30, MEM10, CN10, S10, AZM15, C30, CIP5, TE30, SXT25))
abx_table <- abx_table |>
  cols_label(
    `Isolates.n` = html("Isolates<br>n"),
    `MDR..n.....` = html("MDR<br>n (%)"),
             Site.Category = html("Site Category"))
abx_table <- abx_table |>
  tab_footnote(footnote = "Defined as non-susceptibility to ≥1 agent in ≥3 antimicrobial categories.",
               locations = cells_column_labels(columns = MDR..n.....))

abx_table <- abx_table |>
  opt_stylize(style = 3, color = "blue") |>
  tab_options(table.font.size = px(11)) |>
  cols_align(align = "center")

abx_table

gtsave(abx_table, filename = "abx_table_category.png", expand = 5, vwidth = 1600, vheight = 1000)


