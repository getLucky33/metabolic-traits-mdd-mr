# Generate Figures 1–4 for the Scientific Reports submission.
# Values and evidence classifications are read from released aggregate tables.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(file_arg[[1]])) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
DATA <- file.path(ROOT, "data", "derived")
cli <- commandArgs(trailingOnly = TRUE)
out_at <- match("--out", cli)
if (!is.na(out_at) && out_at == length(cli)) stop("--out requires a directory")
OUT <- if (is.na(out_at)) file.path(ROOT, "results", "figures") else cli[[out_at + 1L]]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- "sans"

COL <- list(
  navy = "#17324D", blue = "#3B6E8F", sky = "#A9C7D8",
  teal = "#2A7F7F", green = "#5A8F65", orange = "#D9894B",
  red = "#B95C5C", purple = "#7A6E9D", ink = "#25313B",
  grey = "#697680", pale = "#EEF2F4", line = "#D5DDE2", white = "#FFFFFF"
)

theme_sr <- function(base_size = 9) {
  theme_minimal(base_family = FONT, base_size = base_size) +
    theme(
      text = element_text(colour = COL$ink),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "#E7ECEF", linewidth = 0.35),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = COL$ink),
      strip.text = element_text(face = "bold", colour = COL$navy, hjust = 0),
      strip.background = element_rect(fill = "#E8F0F4", colour = NA),
      legend.title = element_blank(),
      legend.key.height = unit(3.2, "mm"),
      plot.margin = margin(7, 9, 7, 7)
    )
}

save_gg <- function(plot, stem, width, height) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(stem, "_600dpi.png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
  svg(file.path(OUT, paste0(stem, ".svg")), width = width, height = height,
      family = FONT, bg = "white")
  print(plot)
  dev.off()
}

manifest <- fread(file.path(DATA, "data_manifest.tsv"))
manifest_paths <- file.path(DATA, manifest$file)
if (any(!file.exists(manifest_paths)) ||
    any(tolower(unname(tools::md5sum(manifest_paths))) != tolower(manifest$md5))) {
  stop("A released input is missing or does not match data_manifest.tsv")
}

sot <- fread(file.path(DATA, "study_counts.tsv"))
sotv <- function(key) {
  value <- sot[id == key, value]
  if (length(value) != 1L) stop("Missing or duplicate SOT id: ", key)
  as.numeric(value)
}
fmt <- function(id) format(sotv(id), big.mark = ",", scientific = FALSE, trim = TRUE)

t2 <- fread(file.path(DATA, "evidence_summary.tsv"))
dictionary <- fread(file.path(DATA, "trait_labels.tsv"))
if (nrow(dictionary) != 15L || uniqueN(dictionary$trait_id) != 15L ||
    uniqueN(dictionary$trait_display) != 15L || any(grepl("_", dictionary$trait_display, fixed = TRUE))) {
  stop("Trait display dictionary failed integrity checks")
}
setkey(dictionary, trait_id)
display_name <- function(x) {
  out <- dictionary[J(x), trait_display]
  if (anyNA(out)) stop("Missing display name for: ", paste(x[is.na(out)], collapse = ", "))
  out
}
figure_full_name <- function(x) {
  out <- display_name(x)
  out <- sub(" \\(GlycA\\)$", "", out)
  out <- gsub("\\bVLDL\\b", "very-low-density lipoprotein", out, perl = TRUE)
  out <- gsub("\\bLDL\\b", "low-density lipoprotein", out, perl = TRUE)
  out <- gsub("\\bIDL\\b", "intermediate-density lipoprotein", out, perl = TRUE)
  out <- gsub("\\bHDL\\b", "high-density lipoprotein", out, perl = TRUE)
  if (any(grepl("\\b(HDL|LDL|IDL|VLDL|GlycA)\\b", out, perl = TRUE))) {
    stop("Figure trait label still contains an abbreviation")
  }
  out
}
wrap_trait <- function(x, width = 46L) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}
t2[, trait_display := display_name(trait)]
t2[, trait_figure_label := wrap_trait(figure_full_name(trait))]
stopifnot(uniqueN(t2$trait_figure_label) == nrow(t2))
t2[, .(display_order = .I, trait_id = trait, trait_display = figure_full_name(trait))] |>
  fwrite(file.path(OUT, "Figure_2_3_trait_labels.tsv"), sep = "\t", eol = "\n")
t2[, group_label := fifelse(reporting_group == "priority_reporting",
                             "Priority reporting group", "Secondary reporting group")]

# Figure 1: study inputs and evidence stages.
draw_figure1 <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.97))

  card <- function(x, y, w, h, title, body, header_fill, body_size = 7.3,
                   border = COL$line, title_size = 9.2) {
    grid.roundrect(x = x, y = y, width = w, height = h, r = unit(2.2, "mm"),
                   gp = gpar(fill = COL$white, col = border, lwd = 1.25))
    header_h <- min(0.045, h * 0.30)
    grid.roundrect(x = x, y = y + (h - header_h) / 2, width = w, height = header_h,
                   r = unit(2.2, "mm"), gp = gpar(fill = header_fill, col = border, lwd = 1.25))
    grid.rect(x = x, y = y + (h - header_h) / 2 - 0.007, width = w - 0.003,
              height = 0.014, gp = gpar(fill = header_fill, col = NA))
    grid.text(title, x = x, y = y + (h - header_h) / 2,
              gp = gpar(fontfamily = FONT, fontsize = title_size,
                        fontface = "bold", col = COL$navy))
    grid.text(body, x = x, y = y - header_h * 0.30, just = "centre",
              gp = gpar(fontfamily = FONT, fontsize = body_size,
                        lineheight = 1.20, col = COL$ink))
  }
  connector <- function(x0, y0, x1, y1, arrow_end = TRUE, lwd = 1.45) {
    grid.lines(x = c(x0, x1), y = c(y0, y1),
               gp = gpar(col = COL$grey, lwd = lwd),
               arrow = if (arrow_end) arrow(type = "closed", length = unit(2.1, "mm")) else NULL)
  }
  connector_path <- function(x, y, arrow_end = TRUE, lwd = 1.45) {
    grid.lines(x = x, y = y, gp = gpar(col = COL$grey, lwd = lwd),
               arrow = if (arrow_end) arrow(type = "closed", length = unit(2.1, "mm")) else NULL)
  }

  # Connectors are drawn before cards so their ends remain visually clean.
  connector(0.27, 0.795, 0.46, 0.728)
  connector(0.73, 0.795, 0.54, 0.728)
  connector(0.50, 0.617, 0.50, 0.595, arrow_end = FALSE)
  connector(0.18, 0.595, 0.82, 0.595, arrow_end = FALSE)
  connector(0.18, 0.595, 0.18, 0.565)
  connector(0.50, 0.595, 0.50, 0.565)
  connector(0.82, 0.595, 0.82, 0.565)
  connector(0.18, 0.395, 0.18, 0.365, arrow_end = FALSE)
  connector(0.50, 0.395, 0.50, 0.365, arrow_end = FALSE)
  connector(0.82, 0.395, 0.82, 0.365, arrow_end = FALSE)
  connector(0.18, 0.365, 0.82, 0.365, arrow_end = FALSE)
  connector(0.50, 0.365, 0.50, 0.339)
  connector(0.50, 0.235, 0.50, 0.195)

  card(
    0.27, 0.865, 0.40, 0.14, "Exposure GWAS",
    paste0(fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits\n",
           "Estonian Biobank + UK Biobank\nN = ", fmt("tambets_meta_eur")),
    "#F5DFA7", body_size = 7.7
  )
  card(
    0.73, 0.865, 0.40, 0.14, "Primary depression outcome GWAS",
    paste0("PGC major depression\n", fmt("pgc_main_cases"), " cases / ",
           fmt("pgc_main_controls"), " controls"),
    "#F5DFA7", body_size = 7.7
  )
  card(
    0.50, 0.670, 0.58, 0.105, "Instrument construction and primary forward MR",
    paste0("Genome-wide significant, LD-clumped instruments\n",
           "Primary outcome: PGC major depression\n",
           fmt("n_screen_bonf"), " Bonferroni-significant traits"),
    "#B9DCEB", body_size = 6.8, border = "#8EB9CC", title_size = 8.8
  )

  card(
    0.18, 0.480, 0.28, 0.17, "Robustness and sensitivity",
    "Harmonization and five estimators\nQ, MR-Egger and MR-PRESSO\nPalindromic-variant sensitivity\nPleiotropic-region and shared-IV exclusions\nPGC excluding UK Biobank",
    "#D8E8EF", body_size = 5.8, title_size = 7.7
  )
  card(
    0.50, 0.480, 0.28, 0.17, "Directionality",
    "Steiger directionality test\nReverse MR with two PGC instrument sets\nAlternative prevalence assumptions",
    "#DDD8E9", body_size = 6.4, title_size = 8.4
  )
  card(
    0.82, 0.480, 0.28, 0.17, "Locus evidence",
    paste0("ABF colocalization and regional checks\nMHC / complex-region assessment\n",
           "Exploratory conditional / SuSiE checks\n", fmt("coloc_records"), " analysis-window records\n",
           fmt("h4_evidence_rows"), " with shared-variant\nposterior support"),
    "#F4DCCB", body_size = 5.6, title_size = 8.4
  )
  card(
    0.50, 0.152, 0.68, 0.085, "Separate alternative-outcome analysis",
    paste0("FinnGen R13 register-based depression: ", fmt("finngen_cases"), " cases / ",
           fmt("finngen_controls"), " controls\n", fmt("fg13_dir_same"), "/",
           fmt("matrix_fg13_rows"), " directionally concordant; reporting groups unchanged\n",
           "Not an independent replication"),
    "#D8E7D7", body_size = 5.9, border = "#AFC7AD", title_size = 8.5
  )
  card(
    0.50, 0.287, 0.78, 0.103, "Evidence integration",
    paste0("Groups defined before FinnGen: ", fmt("groups_priority"), " priority traits | ",
           fmt("groups_secondary"), " secondary traits\n",
           fmt("coloc_records"), " records narrowed to ",
           fmt("h4_evidence_rows"), " with shared-variant posterior support\n",
           fmt("mechanism_eligible_true"), " met the composite mechanism-support criterion"),
    "#C9D8E8", body_size = 7.1, border = "#9CB2C8", title_size = 9.1
  )
  popViewport()
}

mermaid <- c(
  "flowchart TB",
  paste0("  E[\"Exposure GWAS<br/>", fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits<br/>N = ", fmt("tambets_meta_eur"), "; Estonian Biobank + UK Biobank\"]"),
  paste0("  O[\"Primary depression outcome GWAS<br/>PGC major depression: ", fmt("pgc_main_cases"), " cases / ", fmt("pgc_main_controls"), " controls\"]"),
  paste0("  P[\"Instrument construction and primary forward MR<br/>Genome-wide significant, LD-clumped instruments<br/>", fmt("n_screen_bonf"), " Bonferroni-significant traits\"]"),
  "  R[\"Robustness and sensitivity<br/>Harmonization and five estimators<br/>Q, MR-Egger, MR-PRESSO and exclusion analyses<br/>PGC outcome excluding UK Biobank\"]",
  "  D[\"Directionality<br/>Steiger directionality test<br/>Reverse Mendelian randomization with two PGC instrument sets<br/>Alternative prevalence assumptions\"]",
  paste0("  L[\"Locus evidence<br/>ABF colocalization; regional, MHC and complex-region checks<br/>Exploratory conditional / SuSiE checks<br/>", fmt("coloc_records"), " analysis-window records; ", fmt("h4_evidence_rows"), " with shared-variant posterior support\"]"),
  paste0("  F[\"Separate FinnGen R13 alternative-outcome analysis<br/>", fmt("finngen_cases"), " cases / ", fmt("finngen_controls"), " controls<br/>", fmt("fg13_dir_same"), "/", fmt("matrix_fg13_rows"), " directionally concordant; reporting groups unchanged<br/>Not an independent replication\"]"),
  paste0("  G[\"Evidence integration<br/>Groups defined before FinnGen: ", fmt("groups_priority"), " priority / ", fmt("groups_secondary"), " secondary<br/>", fmt("coloc_records"), " records; ", fmt("h4_evidence_rows"), " with shared-variant posterior support; ", fmt("mechanism_eligible_true"), " met the composite criterion\"]"),
  "  E --> P",
  "  O --> P",
  "  P --> R",
  "  P --> D",
  "  P --> L",
  "  R --> G",
  "  D --> G",
  "  L --> G",
  "  G --> F"
)
stopifnot(
  all(c("  P --> R", "  P --> D", "  P --> L",
        "  R --> G", "  D --> G", "  L --> G", "  G --> F") %in% mermaid),
  !any(c("  P --> F", "  R --> F", "  D --> F", "  L --> F", "  F --> G") %in% mermaid)
)
writeBin(
  charToRaw(paste0(paste(mermaid, collapse = "\n"), "\n")),
  file.path(OUT, "Figure_1_source.mermaid")
)
cairo_pdf(file.path(OUT, "Figure_1.pdf"), width = 7.2, height = 9.2, family = FONT, bg = "white")
draw_figure1(); dev.off()
png(file.path(OUT, "Figure_1_600dpi.png"), width = 7.2, height = 9.2,
    units = "in", res = 600, type = "cairo", bg = "white")
draw_figure1(); dev.off()
svg(file.path(OUT, "Figure_1.svg"), width = 7.2, height = 9.2, family = FONT, bg = "white")
draw_figure1(); dev.off()

# Figure 2: three-series forest plot on the interpretable odds-ratio scale.
forest <- fread(file.path(DATA, "forest_estimates.tsv"))
forest <- forest[match(t2$trait, trait)]
stopifnot(identical(forest$trait, t2$trait), !anyNA(forest$b_pgc), !anyNA(forest$b_fg))
forward_screen <- fread(file.path(DATA, "forward_screen_249.tsv"))
noukbb_screen <- fread(file.path(DATA, "forward_noukbb_15.tsv"))
finngen_screen <- fread(file.path(DATA, "finngen_cross_outcome_15.tsv"))
if (nrow(noukbb_screen) != 15L || uniqueN(noukbb_screen$trait) != 15L ||
    !setequal(noukbb_screen$trait, t2$trait)) {
  stop("The UK Biobank-excluded Figure 2 source must contain the 15 unique forward candidates")
}
forward_plot_source <- forward_screen[match(t2$trait, trait)]
noukbb_plot_source <- noukbb_screen[match(t2$trait, trait)]
finngen_plot_source <- finngen_screen[match(t2$trait, trait)]
same_num <- function(a, b) all(is.finite(a)) && all(is.finite(b)) &&
  max(abs(as.numeric(a) - as.numeric(b))) < 1e-12
if (!identical(forward_plot_source$trait, t2$trait) ||
    !identical(noukbb_plot_source$trait, t2$trait) ||
    !identical(finngen_plot_source$trait, t2$trait) ||
    !same_num(forest$b_pgc, forward_plot_source$ivw_b) ||
    !same_num(forest$se_pgc, forward_plot_source$ivw_se) ||
    !same_num(forest$b_fg, finngen_plot_source$fg_b) ||
    !same_num(forest$se_fg, finngen_plot_source$fg_se) ||
    !identical(forest$cross_outcome_support_label,
               finngen_plot_source$cross_outcome_support_label)) {
  stop("Figure 2 values do not match the released primary/FinnGen source tables")
}
f2_source <- rbind(
  forest[, .(trait, reporting_group, dataset = "PGC major depression", n_iv = forward_plot_source$n_iv,
             beta = b_pgc, se = se_pgc, p_value = forward_plot_source$ivw_p)],
  noukbb_plot_source[, .(trait, reporting_group = forest$reporting_group,
                         dataset = "PGC excluding UK Biobank", n_iv,
                         beta = ivw_b, se = ivw_se, p_value = ivw_p)],
  forest[, .(trait, reporting_group, dataset = "FinnGen R13 depression", n_iv = finngen_plot_source$n_iv,
             beta = b_fg, se = se_fg, p_value = finngen_plot_source$fg_p)]
)
if (nrow(f2_source) != 45L || uniqueN(f2_source, by = c("trait", "dataset")) != 45L ||
    any(!is.finite(f2_source$n_iv)) || any(f2_source$n_iv <= 0) ||
    any(!is.finite(f2_source$beta)) || any(!is.finite(f2_source$se)) || any(f2_source$se <= 0) ||
    any(!is.finite(f2_source$p_value)) || any(f2_source$p_value <= 0 | f2_source$p_value > 1) ||
    any(abs(f2_source$p_value - 2 * pnorm(-abs(f2_source$beta / f2_source$se))) > 1e-12) ||
    !identical(sign(noukbb_plot_source$ivw_b) == sign(forest$b_pgc), as.logical(t2$noukbb_ok)) ||
    sum(noukbb_plot_source$ivw_p < 0.05) != 15L ||
    sum(noukbb_plot_source$ivw_p < 0.05 / 249) != 10L) {
  stop("Figure 2 three-dataset source table failed integrity checks")
}
f2 <- f2_source[, .(trait, reporting_group, dataset, n_iv, beta, se, p_value,
                    estimate = exp(beta), lower = exp(beta - 1.96 * se),
                    upper = exp(beta + 1.96 * se))]
f2[, trait_display := wrap_trait(figure_full_name(trait))]
f2[, group_label := fifelse(reporting_group == "priority_reporting",
                             "Priority reporting group", "Secondary reporting group")]
fwrite(
  f2[, .(trait_id = trait, trait_display = figure_full_name(trait), reporting_group, group_label,
         dataset, n_iv, beta, se, p_value, estimate, lower, upper)],
  file.path(OUT, "Figure_2_plot_data.tsv"), sep = "\t", eol = "\n"
)
f2[, trait_display := factor(trait_display, levels = rev(t2$trait_figure_label))]
dataset_levels <- c("PGC major depression", "PGC excluding UK Biobank", "FinnGen R13 depression")
f2[, dataset := factor(dataset, levels = dataset_levels)]
f2[, estimate_text := sprintf("%.3f (%.3f-%.3f)", estimate, lower, upper)]
pd2 <- position_dodge(width = 0.62, reverse = TRUE)
p2_forest <- ggplot(f2, aes(estimate, trait_display, colour = dataset, shape = dataset)) +
  geom_vline(xintercept = 1, colour = "#7F8B93", linewidth = 0.5, linetype = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0,
                position = pd2, linewidth = 0.55) +
  geom_point(position = pd2, size = 2.15, stroke = 0.45) +
  facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(COL$navy, COL$blue, COL$orange), breaks = dataset_levels) +
  scale_shape_manual(values = c(16, 15, 17), breaks = dataset_levels) +
  scale_x_log10(breaks = c(0.88, 0.92, 0.96, 1.00, 1.04, 1.08, 1.12),
                limits = c(0.855, 1.155)) +
  labs(x = "Odds ratio per genetically predicted 1-SD increase (95% CI)", y = NULL) +
  theme_sr(9.3) +
  theme(legend.position = "top", legend.justification = "left",
        axis.text.y = element_text(size = 8.1), panel.spacing.y = unit(3, "mm"))
p2_values <- ggplot(f2, aes(0, trait_display, colour = dataset, group = dataset)) +
  geom_text(aes(label = estimate_text), position = pd2, hjust = 0,
            family = FONT, size = 2.55, show.legend = FALSE) +
  facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(COL$navy, COL$blue, COL$orange), breaks = dataset_levels) +
  scale_x_continuous(limits = c(0, 1), breaks = 0, labels = "OR (95% CI)", position = "top") +
  coord_cartesian(clip = "off") +
  theme_sr(9.3) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.text.x.top = element_text(face = "bold", hjust = 0),
        strip.text = element_text(colour = NA), strip.background = element_rect(fill = NA, colour = NA),
        legend.position = "none", panel.spacing.y = unit(3, "mm"),
        plot.margin = margin(7, 24, 7, 7))
p2 <- p2_forest + p2_values + plot_layout(widths = c(0.69, 0.31), guides = "collect") &
  theme(legend.position = "top")
save_gg(p2, "Figure_2", 13.0, 9.3)

# Figure 3: evidence matrix. Text carries each category; colour repeats it.
rev_main <- fread(file.path(DATA, "reverse_mr_primary.tsv"))
rev_noukbb <- fread(file.path(DATA, "reverse_mr_no_ukbb.tsv"))
setkey(rev_main, trait); setkey(rev_noukbb, trait)
rm <- rev_main[J(t2$trait)]; rn <- rev_noukbb[J(t2$trait)]
if (anyNA(rm$screen_level) || anyNA(rn$screen_level)) stop("Reverse-MR matrix source is incomplete")
reverse_main_full <- fread(file.path(DATA, "reverse_screen_main_249.tsv"))
reverse_noukbb_full <- fread(file.path(DATA, "reverse_screen_noukbb_249.tsv"))
setkey(reverse_main_full, trait); setkey(reverse_noukbb_full, trait)
if (!identical(rm$screen_level, reverse_main_full[J(t2$trait), screen_level]) ||
    !identical(rn$screen_level, reverse_noukbb_full[J(t2$trait), screen_level])) {
  stop("Figure 3 reverse-MR labels do not match the released 249-trait screens")
}
shared_iv <- fread(file.path(DATA, "broad_pleiotropy_sensitivity_15.tsv"))
setkey(shared_iv, trait_id)
if (!identical(t2$pleio_ok, shared_iv[J(t2$trait), after_p] < 0.05)) {
  stop("Figure 3 shared-IV-removal labels do not match the released sensitivity table")
}
reverse_code <- ifelse(rm$screen_level == "bonferroni_hit" & rn$screen_level == "bonferroni_hit", "B/B",
                       ifelse(rm$screen_level == "bonferroni_hit", "B/F",
                              ifelse(rn$screen_level == "bonferroni_hit", "F/B", "F/F")))
main_or <- as.numeric(sub(" .*", "", t2$main_or95))
fg_or <- as.numeric(sub(" .*", "", t2$fg13_or95))
coloc_code <- ifelse(t2$coloc_n_robust > 0, "R", ifelse(t2$coloc_n_prior > 0, "P",
                     ifelse(t2$coloc_n_distinct > 0, "D", "O")))
matrix_wide <- data.table(
  trait = t2$trait,
  `Effect\nIVW dir.` = ifelse(main_or > 1, "+", "-"),
  `Effect\n5 methods` = paste0(t2$five_dir, "/5"),
  `Robustness\nPalindromic exclusion` = ifelse(t2$a3_ok == TRUE, "Yes", "No"),
  `Robustness\nUK Biobank excluded` = ifelse(t2$noukbb_ok == TRUE, "Yes", "No"),
  `Robustness\nMR-PRESSO` = ifelse(t2$presso == TRUE, "Yes", "No"),
  `Sensitivity\nShared-IV removal` = ifelse(t2$pleio_ok == TRUE, "Yes", "No"),
  `Diagnostics\nEgger int.` = ifelse(t2$egger_sig == TRUE, "Sig", "NS"),
  `Diagnostics\nCochran Q` = ifelse(t2$q_sig == TRUE, "Sig", "NS"),
  `Direction\nReverse MR` = reverse_code,
  `Alternative outcome\nFG direction` = ifelse(fg_or > 1, "+", "-"),
  `Alternative outcome\nP level` = fifelse(t2$fg13_p_bh < 0.05, "FDR",
                       fifelse(t2$fg13_p < 0.05, "Nom", "NS")),
  `Locus\nColoc` = coloc_code
)
layer_order <- names(matrix_wide)[-1]
m3 <- melt(matrix_wide, id.vars = "trait", variable.name = "layer", value.name = "value")
m3[, trait_display := wrap_trait(figure_full_name(trait))]
m3[, group_label := t2$group_label[match(trait, t2$trait)]]
m3[, trait_display := factor(trait_display, levels = rev(t2$trait_figure_label))]
m3[, layer := factor(layer, levels = layer_order)]
m3[, state := fcase(
  value %chin% c("Yes", "5/5", "FDR", "R", "B/B"), "support",
  value %chin% c("Nom", "P", "B/F", "F/B"), "qualified",
  value %chin% c("No", "NS", "Sig"), "limitation",
  value %chin% c("+"), "positive",
  value %chin% c("-"), "negative",
  default = "neutral"
)]
fwrite(
  m3[, .(
    trait_id = trait,
    trait_display = figure_full_name(trait),
    reporting_group = t2$reporting_group[match(trait, t2$trait)],
    group_label,
    layer = gsub("\n", " / ", as.character(layer), fixed = TRUE),
    value = as.character(value),
    state,
    definition = fifelse(
      as.character(layer) == "Sensitivity\nShared-IV removal",
      "Yes means P<0.05 after removing variants used as instruments for more than five metabolic traits",
      ""
    )
  )],
  file.path(OUT, "Figure_3_matrix_data.tsv"), sep = "\t", eol = "\n"
)
p3 <- ggplot(m3, aes(layer, trait_display, fill = state, label = value)) +
  geom_tile(colour = "white", linewidth = 0.75) +
  geom_text(family = FONT, size = 2.75, colour = COL$ink) +
  facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(support = "#DCEBDD", qualified = "#F6E8C9",
                               limitation = "#F2D7D5", positive = "#D9E8F2",
                               negative = "#E9E2F2", neutral = "#EDF0F2"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_sr(8.6) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 7.7), panel.spacing.y = unit(3, "mm")) +
  geom_vline(xintercept = c(2.5, 6.5, 8.5, 10.5), colour = "#AAB5BC", linewidth = 0.35)
save_gg(p3, "Figure_3", 12.2, 8.5)

# Figure 4: evidence attrition with exact counts; panel area does not encode counts.
stage1 <- data.table(
  class = factor(c("Robust", "Prior-sensitive", "Distinct signal", "Trait-specific/low power", "Inconclusive"),
                 levels = c("Robust", "Prior-sensitive", "Distinct signal", "Trait-specific/low power", "Inconclusive")),
  n = c(sotv("n_coloc_robust"), sotv("n_coloc_prior"), sotv("n_coloc_distinct"),
        sotv("n_coloc_tsp"), sotv("n_coloc_inconcl"))
)
stage2 <- data.table(
  class = factor(c("Complex-downgraded", "Forward-unassessable", "Prior-forward-consistent", "Forward-opposite"),
                 levels = c("Complex-downgraded", "Forward-unassessable", "Prior-forward-consistent", "Forward-opposite")),
  n = c(sotv("tier_complex"), sotv("tier_unassess"), sotv("tier_prior_consistent"), sotv("tier_opposite"))
)
stopifnot(sum(stage1$n) == sotv("coloc_records"), sum(stage2$n) == sotv("h4_evidence_rows"),
          sotv("mechanism_eligible_true") == 0)
pal1 <- c("#2A7F7F", "#79A8A8", "#D4A35E", "#8DA0AE", "#C5CDD2")
pal2 <- c("#7A6E9D", "#A7A0BF", "#5A8F65", "#B95C5C")
p4a <- ggplot(stage1, aes(x = 1, y = n, fill = class)) +
  geom_col(width = 0.52, colour = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = pal1, labels = paste0(stage1$class, " (", stage1$n, ")")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = scales::comma) +
  labs(title = paste0("a  All colocalization records (n = ", fmt("coloc_records"), ")"), x = NULL, y = NULL) +
  theme_sr(8.6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.4),
        plot.title = element_text(face = "bold", size = 9.5)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p4b <- ggplot(stage2, aes(x = 1, y = n, fill = class)) +
  geom_col(width = 0.52, colour = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = pal2, labels = paste0(stage2$class, " (", stage2$n, ")")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(title = paste0("b  Records with shared-causal-variant posterior support (n = ", fmt("h4_evidence_rows"), ")"), x = NULL, y = NULL) +
  theme_sr(8.6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.4),
        plot.title = element_text(face = "bold", size = 9.5)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p4c <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#F7F4EC", colour = "#DCCDAA", linewidth = 0.8) +
  annotate("text", x = 0.5, y = 0.67, label = paste0(fmt("mechanism_eligible_true"), "/", fmt("h4_evidence_rows")),
           family = FONT, fontface = "bold", colour = COL$navy, size = 11) +
  annotate("text", x = 0.5, y = 0.34, label = "met the composite mechanism-support criterion",
           family = FONT, colour = COL$ink, size = 3.6) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title = "c  Composite mechanism support") +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 9.5, colour = COL$ink),
        plot.margin = margin(7, 9, 20, 7))
p4 <- (p4a / p4b / p4c) + plot_layout(heights = c(1.2, 1.2, 0.85))
save_gg(p4, "Figure_4", 10.8, 8.0)

cat("Scientific Reports Figures 1-4 regenerated from frozen sources\n")
cat(list.files(OUT, pattern = "^Figure_[1-4]", full.names = FALSE), sep = "\n")


expected <- c(
  file.path(OUT, paste0("Figure_", rep(1:4, each = 3),
                       rep(c(".pdf", ".svg", "_600dpi.png"), times = 4))),
  file.path(OUT, "Figure_1_source.mermaid"),
  file.path(OUT, "Figure_2_3_trait_labels.tsv"),
  file.path(OUT, "Figure_2_plot_data.tsv"),
  file.path(OUT, "Figure_3_matrix_data.tsv")
)
if (!all(file.exists(expected)) || any(file.info(expected)$size < 1000)) {
  stop("One or more expected figure outputs are absent or unexpectedly small")
}
checksums <- data.table(
  file = basename(expected),
  bytes = file.info(expected)$size,
  md5 = unname(tools::md5sum(expected))
)
fwrite(checksums, file.path(OUT, "figure_checksums.tsv"), sep = "\t", eol = "\n")
cat("Reproduction checks passed; checksum manifest written to the selected output directory\n")
