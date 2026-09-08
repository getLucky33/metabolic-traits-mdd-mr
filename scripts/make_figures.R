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
  fig <- list(
    navy = "#17324D", ink = "#26333D", slate = "#667782", line = "#C7D2D9",
    pale = "#F5F7F8", white = "#FFFFFF", yellow = "#F3DFA7", yellow_line = "#D6B86E",
    blue = "#B9DCEB", blue_line = "#7FAFC4", teal = "#DCECEF", teal_line = "#93BAC2",
    violet = "#E4DFF0", violet_line = "#AAA0C2", peach = "#F4DFD1", peach_line = "#D7A98A",
    synthesis = "#D5E2EF", synthesis_line = "#8FAAC2", green = "#DDEBD9",
    green_line = "#90B18B", muted_blue = "#E8F0F4"
  )
  support_counts <- table(factor(
    t2$cross_outcome_support_label,
    levels = c("FDR-supported", "nominally-supported", "no-nominal-support")
  ))
  stopifnot(identical(as.integer(support_counts), c(11L, 1L, 3L)))

  rounded_box <- function(x, y, w, h, fill = fig$white, border = fig$line,
                          lwd = 1.2, radius = 2.2, lty = 1) {
    grid.roundrect(
      x = x, y = y, width = w, height = h, r = unit(radius, "mm"),
      gp = gpar(fill = fill, col = border, lwd = lwd, lty = lty)
    )
  }
  text_at <- function(label, x, y, size = 8, face = "plain", colour = fig$ink,
                      just = "centre", lineheight = 1.13) {
    grid.text(
      label, x = x, y = y, just = just,
      gp = gpar(fontfamily = FONT, fontsize = size, fontface = face,
                col = colour, lineheight = lineheight)
    )
  }
  connector <- function(x0, y0, x1, y1, dashed = FALSE, end = TRUE, lwd = 1.35) {
    grid.lines(
      x = c(x0, x1), y = c(y0, y1),
      gp = gpar(col = fig$slate, lwd = lwd, lty = if (dashed) 2 else 1),
      arrow = if (end) arrow(type = "closed", length = unit(1.9, "mm")) else NULL
    )
  }
  stage_heading <- function(label, y) {
    text_at(toupper(label), 0.070, y, size = 7.0, face = "bold", colour = fig$slate,
            just = "left")
  }
  header_card <- function(x, y, w, h, title, body, fill, border, body_size = 7.2) {
    rounded_box(x, y, w, h, fill = fig$white, border = border, lwd = 1.25)
    header_h <- 0.041
    grid.roundrect(
      x = x, y = y + (h - header_h) / 2, width = w, height = header_h,
      r = unit(2.2, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.25)
    )
    grid.rect(x = x, y = y + (h - header_h) / 2 - 0.0065,
              width = w - 0.003, height = 0.013, gp = gpar(fill = fill, col = NA))
    text_at(title, x, y + (h - header_h) / 2, size = 8.8, face = "bold", colour = fig$navy)
    text_at(body, x, y - header_h * 0.32, size = body_size, lineheight = 1.17)
  }
  evaluation_card <- function(x, title, lines, fill, border, footer = NULL) {
    y <- 0.415
    w <- 0.285
    h <- 0.200
    rounded_box(x, y, w, h, fill = fig$white, border = border, lwd = 1.2)
    grid.roundrect(x = x, y = y + 0.0805, width = w, height = 0.039,
                   r = unit(2.2, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.2))
    grid.rect(x = x, y = y + 0.074, width = w - 0.003, height = 0.013,
              gp = gpar(fill = fill, col = NA))
    text_at(title, x, y + 0.0805, size = 8.2, face = "bold", colour = fig$navy)
    body_y <- if (is.null(footer)) y - 0.002 else y + 0.017
    text_at(paste(lines, collapse = "\n"), x, body_y,
            size = if (is.null(footer)) 6.7 else 6.25, lineheight = 1.16)
    if (!is.null(footer)) {
      rounded_box(x, y - 0.062, w - 0.026, 0.052, fill = fig$pale, border = fig$line,
                  lwd = 0.85, radius = 1.5)
      text_at(footer, x, y - 0.062, size = 5.05, lineheight = 1.04)
    }
  }

  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.94, height = 0.965))

  # Fill-only stage bands preserve reading order without numbered badges or crossing borders.
  grid.roundrect(x = 0.5175, y = 0.883, width = 0.925, height = 0.145,
                 r = unit(2.0, "mm"), gp = gpar(fill = "#FBFCFC", col = NA))
  grid.roundrect(x = 0.5175, y = 0.700, width = 0.925, height = 0.200,
                 r = unit(2.0, "mm"), gp = gpar(fill = "#F6FAFC", col = NA))
  grid.roundrect(x = 0.5175, y = 0.445, width = 0.925, height = 0.250,
                 r = unit(2.0, "mm"), gp = gpar(fill = "#FCFCFD", col = NA))
  grid.roundrect(x = 0.5175, y = 0.220, width = 0.925, height = 0.150,
                 r = unit(2.0, "mm"), gp = gpar(fill = "#F7F9FB", col = NA))
  grid.roundrect(x = 0.5175, y = 0.060, width = 0.925, height = 0.105,
                 r = unit(2.0, "mm"), gp = gpar(fill = "#F8FBF7", col = NA))

  stage_heading("GWAS inputs", 0.946)
  stage_heading("Bidirectional MR screening", 0.792)
  stage_heading("Complementary evidence assessment", 0.558)
  stage_heading("Evidence synthesis", 0.285)
  stage_heading("Alternative outcome", 0.105)

  # Connectors are drawn before cards so their ends remain visually clean.
  connector(0.29, 0.825, 0.40, 0.770)
  connector(0.76, 0.825, 0.66, 0.770)
  connector(0.43, 0.657, 0.43, 0.642)
  connector(0.53, 0.605, 0.53, 0.530, end = FALSE)
  connector(0.205, 0.530, 0.855, 0.530, end = FALSE)
  connector(0.205, 0.530, 0.205, 0.515)
  connector(0.53, 0.530, 0.53, 0.515)
  connector(0.855, 0.530, 0.855, 0.515)
  connector(0.205, 0.315, 0.205, 0.300, end = FALSE)
  connector(0.53, 0.315, 0.53, 0.300, end = FALSE)
  connector(0.855, 0.315, 0.855, 0.300, end = FALSE)
  connector(0.205, 0.300, 0.855, 0.300, end = FALSE)
  connector(0.53, 0.300, 0.53, 0.263)
  connector(0.53, 0.153, 0.53, 0.086, dashed = TRUE)

  header_card(
    0.285, 0.880, 0.405, 0.105, "Exposure GWAS",
    paste0(fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits\n",
           "Estonian Biobank and UK Biobank   N = ", fmt("tambets_meta_eur")),
    fig$yellow, fig$yellow_line
  )
  header_card(
    0.755, 0.880, 0.405, 0.105, "Primary depression outcome GWAS",
    paste0("PGC major depression\n", fmt("pgc_main_cases"), " cases   ",
           fmt("pgc_main_controls"), " controls"),
    fig$yellow, fig$yellow_line
  )

  rounded_box(0.53, 0.700, 0.77, 0.140, fill = fig$white, border = fig$blue_line, lwd = 1.3)
  grid.roundrect(x = 0.53, y = 0.751, width = 0.77, height = 0.034,
                 r = unit(2.0, "mm"), gp = gpar(fill = fig$blue, col = fig$blue_line, lwd = 1.3))
  grid.rect(x = 0.53, y = 0.745, width = 0.767, height = 0.011,
            gp = gpar(fill = fig$blue, col = NA))
  text_at("Bidirectional Mendelian randomization screening",
          0.53, 0.751, size = 8.7, face = "bold", colour = fig$navy)

  rounded_box(0.35, 0.687, 0.34, 0.059, fill = fig$muted_blue,
              border = fig$blue_line, lwd = 0.9, radius = 1.5)
  text_at("Primary forward MR screen", 0.35, 0.701, size = 6.9,
          face = "bold", colour = fig$navy)
  text_at(paste0("All ", fmt("n_traits_metabolites"), " traits as exposures; primary PGC outcome\n",
                 "MAF-specific thresholds, F at least 10, LD clumping and random-effects IVW"),
          0.35, 0.678, size = 4.65, lineheight = 1.04)

  rounded_box(0.72, 0.687, 0.34, 0.059, fill = "#F0EDF6",
              border = fig$violet_line, lwd = 0.9, radius = 1.5)
  text_at("Separate reverse-MR screens", 0.72, 0.701, size = 6.9,
          face = "bold", colour = fig$navy)
  text_at(paste0("PGC liability as exposure; all ", fmt("n_traits_metabolites"), " traits as outcomes\n",
                 "Primary and UK Biobank-excluded PGC instrument sets\n",
                 "Corresponding candidate subset used for directionality"),
          0.72, 0.676, size = 4.55, lineheight = 1.03)

  rounded_box(0.53, 0.623, 0.35, 0.036, fill = fig$blue,
              border = fig$blue_line, lwd = 1.1, radius = 1.5)
  text_at(paste0(fmt("n_screen_bonf"), " Bonferroni-significant forward candidates"),
          0.53, 0.623, size = 6.9, face = "bold", colour = fig$navy)

  evaluation_card(
    0.205, "Robustness and sensitivity",
    c("Same 15 forward candidates", "Five-estimator comparison",
      "Cochran's Q, MR-Egger and MR-PRESSO", "Alternative harmonization",
      "Regional and shared-instrument exclusions"),
    fig$teal, fig$teal_line,
    paste0("UK Biobank-excluded PGC sensitivity\n", fmt("pgc_nb_cases"), " cases   ",
           fmt("pgc_nb_controls"), " controls\nSame 15 summarized; full 249-trait screen reported")
  )
  evaluation_card(
    0.53, "Directionality",
    c("Same 15 forward candidates", "Steiger assessment at three prevalences",
      "Candidate-linked reverse-MR comparison", "Corresponding results from the two",
      "complete 249-trait reverse screens"),
    fig$violet, fig$violet_line
  )
  evaluation_card(
    0.855, "Locus evidence",
    c("Same 15 forward candidates", "Analysis-window construction",
      "ABF colocalization in both PGC outcomes", "Regional direction checks",
      "MHC and complex-region assessment"),
    fig$peach, fig$peach_line,
    paste0(fmt("coloc_records"), " trait-specific analysis-window records")
  )

  rounded_box(0.53, 0.205, 0.77, 0.105, fill = fig$white,
              border = fig$synthesis_line, lwd = 1.3)
  grid.roundrect(x = 0.53, y = 0.244, width = 0.77, height = 0.032,
                 r = unit(2.0, "mm"), gp = gpar(fill = fig$synthesis,
                                                col = fig$synthesis_line, lwd = 1.3))
  grid.rect(x = 0.53, y = 0.239, width = 0.767, height = 0.010,
            gp = gpar(fill = fig$synthesis, col = NA))
  text_at("Pre-FinnGen evidence synthesis", 0.53, 0.244, size = 8.5,
          face = "bold", colour = fig$navy)
  grid.lines(x = c(0.53, 0.53), y = c(0.163, 0.230), gp = gpar(col = fig$line, lwd = 1.0))
  text_at("Trait-level reporting", 0.345, 0.222, size = 7.2, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("groups_priority"), " priority-reporting traits\n",
                 fmt("groups_secondary"), " secondary-reporting traits"),
          0.345, 0.191, size = 6.8, face = "bold", colour = fig$navy, lineheight = 1.10)
  text_at("Groups defined before FinnGen", 0.345, 0.161, size = 5.8, colour = fig$slate)
  text_at("Locus-level assessment", 0.715, 0.222, size = 7.2, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("h4_evidence_rows"), " with shared-variant posterior support"),
          0.715, 0.198, size = 6.4, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("mechanism_eligible_true"), "/", fmt("h4_evidence_rows"),
                 " records met the composite mechanism-support criterion"),
          0.715, 0.169, size = 5.55, colour = fig$slate)

  rounded_box(0.53, 0.050, 0.70, 0.065, fill = fig$white,
              border = fig$green_line, lwd = 1.25, lty = 2)
  grid.roundrect(x = 0.53, y = 0.074, width = 0.70, height = 0.026,
                 r = unit(1.8, "mm"), gp = gpar(fill = fig$green,
                                                col = fig$green_line, lwd = 1.2, lty = 2))
  grid.rect(x = 0.53, y = 0.070, width = 0.697, height = 0.008,
            gp = gpar(fill = fig$green, col = NA))
  text_at("Separate FinnGen R13 alternative-outcome comparison",
          0.53, 0.074, size = 7.4, face = "bold", colour = fig$navy)
  text_at(
    paste0(fmt("finngen_cases"), " cases   ", fmt("finngen_controls"), " controls     ",
           fmt("fg13_dir_same"), "/", fmt("matrix_fg13_rows"), " directionally concordant\n",
           support_counts[["FDR-supported"]], " FDR-supported   ",
           support_counts[["nominally-supported"]], " nominally supported   ",
           support_counts[["no-nominal-support"]], " without nominal support\n",
           "Reporting groups unchanged; not an independent replication"),
    0.53, 0.035, size = 5.05, lineheight = 1.06
  )

  popViewport()
}

mermaid <- c(
  "flowchart TB",
  paste0("  E[\"Exposure GWAS<br/>", fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits<br/>N = ", fmt("tambets_meta_eur"), "; Estonian Biobank and UK Biobank\"]"),
  paste0("  O[\"Primary depression outcome GWAS<br/>PGC major depression: ", fmt("pgc_main_cases"), " cases / ", fmt("pgc_main_controls"), " controls\"]"),
  paste0("  N[\"UK Biobank-excluded PGC sensitivity outcome<br/>", fmt("pgc_nb_cases"), " cases / ", fmt("pgc_nb_controls"), " controls\"]"),
  "  FW[\"Primary forward MR screen<br/>All 249 metabolic traits as exposures\"]",
  "  RV[\"Separate reverse-MR screens<br/>All 249 metabolic traits as outcomes; two PGC instrument sets<br/>Corresponding candidate subset used for directionality\"]",
  paste0("  C[\"", fmt("n_screen_bonf"), " Bonferroni-significant forward candidates\"]"),
  "  R[\"Robustness and sensitivity<br/>Same 15 forward candidates; five estimators; Cochran's Q, MR-Egger and MR-PRESSO<br/>UK Biobank-excluded sensitivity; regional and shared-instrument exclusions\"]",
  "  D[\"Directionality<br/>Steiger assessment; candidate-linked results from both complete 249-trait reverse-MR screens\"]",
  paste0("  L[\"Locus evidence<br/>Same 15 forward candidates; ABF colocalization and regional direction checks<br/>MHC and complex-region assessment; ", fmt("coloc_records"), " trait-specific analysis-window records\"]"),
  paste0("  G[\"Pre-FinnGen evidence synthesis<br/>", fmt("groups_priority"), " priority-reporting / ", fmt("groups_secondary"), " secondary-reporting traits<br/>", fmt("h4_evidence_rows"), " with shared-variant posterior support; ", fmt("mechanism_eligible_true"), "/", fmt("h4_evidence_rows"), " met the composite criterion\"]"),
  paste0("  F[\"Separate FinnGen R13 alternative-outcome comparison<br/>", fmt("fg13_dir_same"), "/", fmt("matrix_fg13_rows"), " directionally concordant; reporting groups unchanged<br/>Not an independent replication\"]"),
  "  E --> FW",
  "  O --> FW",
  "  O --> RV",
  "  N --> RV",
  "  N --> R",
  "  FW --> C",
  "  C --> R",
  "  C --> D",
  "  C --> L",
  "  R --> G",
  "  D --> G",
  "  L --> G",
  "  G -.-> F"
)
stopifnot(
  all(c("  FW --> C", "  C --> R", "  C --> D", "  C --> L",
        "  R --> G", "  D --> G", "  L --> G", "  G -.-> F") %in% mermaid),
  !any(c("  RV --> R", "  RV --> D", "  RV --> L", "  F --> G") %in% mermaid)
)
writeBin(
  charToRaw(paste0(paste(mermaid, collapse = "\n"), "\n")),
  file.path(OUT, "Figure_1_source.mermaid")
)
cairo_pdf(file.path(OUT, "Figure_1.pdf"), width = 7.2, height = 9.4, family = FONT, bg = "white")
draw_figure1(); dev.off()
png(file.path(OUT, "Figure_1_600dpi.png"), width = 7.2, height = 9.4,
    units = "in", res = 600, type = "cairo", bg = "white")
draw_figure1(); dev.off()
svg(file.path(OUT, "Figure_1.svg"), width = 7.2, height = 9.4, family = FONT, bg = "white")
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

# Figure 4: colocalization classifications and the disposition of the 101-record subset.
classification <- fread(file.path(DATA, "coloc_classification.tsv"))
integrated <- fread(file.path(DATA, "integrated_mechanism_evidence.tsv"))
if (nrow(classification) != 6434L || uniqueN(classification$locus_key) != 6434L) {
  stop("The colocalization classification must contain 6,434 unique records")
}
if (nrow(integrated) != 101L || uniqueN(integrated$locus_key) != 101L ||
    !all(integrated$locus_key %chin% classification$locus_key)) {
  stop("The integrated review must contain 101 unique records drawn from the classification table")
}
if (any(as.logical(integrated$mechanism_eligible), na.rm = TRUE) ||
    anyNA(as.logical(integrated$mechanism_eligible))) {
  stop("Mechanism-support status must be FALSE for all 101 records")
}

class_order <- c(
  "robust_coloc", "prior_sensitive_coloc", "distinct_signal",
  "trait_specific_or_low_power", "inconclusive"
)
class_labels <- c(
  robust_coloc = "Robust ABF",
  prior_sensitive_coloc = "Prior-sensitive ABF",
  distinct_signal = "Distinct signal",
  trait_specific_or_low_power = "Trait-specific or low-power",
  inconclusive = "Inconclusive"
)
class_colours <- c(
  robust_coloc = "#146C72",
  prior_sensitive_coloc = "#79B7B0",
  distinct_signal = "#D6A04B",
  trait_specific_or_low_power = "#8FA3B0",
  inconclusive = "#C9D1D6"
)
class_expected <- c(
  robust_coloc = 14L,
  prior_sensitive_coloc = 87L,
  distinct_signal = 822L,
  trait_specific_or_low_power = 3430L,
  inconclusive = 2081L
)
class_counts <- classification[, .N, by = abf_class]
class_values <- setNames(class_counts$N, class_counts$abf_class)[class_order]
if (anyNA(class_values) || !identical(as.integer(class_values), unname(class_expected)) ||
    sum(class_values) != sotv("coloc_records")) {
  stop("The five colocalization-class counts do not match the frozen record")
}

tier_order <- c(
  "abf_label_only_complex_downgraded",
  "abf_label_only_forward_unassessable",
  "abf_prior_sensitive_forward_consistent",
  "abf_label_only_forward_opposite"
)
tier_colours <- c(
  abf_label_only_complex_downgraded = "#7A6E9D",
  abf_label_only_forward_unassessable = "#9DA8B1",
  abf_prior_sensitive_forward_consistent = "#4E79A7",
  abf_label_only_forward_opposite = "#B75A5A"
)
tier_expected <- c(
  abf_label_only_complex_downgraded = 55L,
  abf_label_only_forward_unassessable = 35L,
  abf_prior_sensitive_forward_consistent = 10L,
  abf_label_only_forward_opposite = 1L
)
tier_counts <- integrated[, .N, by = evidence_tier]
tier_values <- setNames(tier_counts$N, tier_counts$evidence_tier)[tier_order]
if (anyNA(tier_values) || !identical(as.integer(tier_values), unname(tier_expected)) ||
    sum(tier_values) != sotv("h4_evidence_rows")) {
  stop("The four integrated-disposition counts do not match the frozen record")
}

flow_expected <- data.table(
  abf_class = c(
    "robust_coloc", "robust_coloc",
    "prior_sensitive_coloc", "prior_sensitive_coloc",
    "prior_sensitive_coloc", "prior_sensitive_coloc"
  ),
  evidence_tier = c(
    "abf_label_only_complex_downgraded",
    "abf_label_only_forward_unassessable",
    "abf_label_only_complex_downgraded",
    "abf_label_only_forward_unassessable",
    "abf_prior_sensitive_forward_consistent",
    "abf_label_only_forward_opposite"
  ),
  n = c(4L, 10L, 51L, 25L, 10L, 1L)
)
flow_observed <- integrated[, .N, by = .(abf_class, evidence_tier)]
setnames(flow_observed, "N", "n")
setorder(flow_expected, abf_class, evidence_tier)
setorder(flow_observed, abf_class, evidence_tier)
if (!identical(flow_observed, flow_expected)) {
  stop("The ABF-class to integrated-disposition cross-tabulation has changed")
}

figure4_plot_data <- rbindlist(list(
  data.table(
    panel = "all_records", source = class_order, destination = NA_character_,
    n = as.integer(class_values), denominator = 6434L,
    percent = 100 * as.integer(class_values) / 6434
  ),
  flow_expected[, .(
    panel = "integrated_flow", source = abf_class, destination = evidence_tier,
    n, denominator = 101L, percent = 100 * n / 101
  )],
  data.table(
    panel = "composite_assessment", source = "mechanism_eligible",
    destination = NA_character_, n = 0L, denominator = 101L, percent = 0
  )
), fill = TRUE)
fwrite(figure4_plot_data, file.path(OUT, "Figure_4_plot_data.tsv"), sep = "\t", eol = "\n")

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)

draw_ribbon <- function(x0, x1, lo0, hi0, lo1, hi1, colour) {
  t <- seq(0, 1, length.out = 80L)
  smooth <- 3 * t^2 - 2 * t^3
  x <- x0 + (x1 - x0) * t
  upper <- hi0 + (hi1 - hi0) * smooth
  lower <- lo0 + (lo1 - lo0) * smooth
  grid.polygon(
    x = unit(c(x, rev(x)), "npc"), y = unit(c(upper, rev(lower)), "npc"),
    gp = gpar(fill = adjustcolor(colour, alpha.f = 0.34),
              col = adjustcolor(colour, alpha.f = 0.62), lwd = 0.6)
  )
}

draw_node <- function(x0, x1, lo, hi, fill) {
  grid.roundrect(
    x = (x0 + x1) / 2, y = (lo + hi) / 2,
    width = x1 - x0, height = hi - lo, r = unit(1.0, "mm"),
    gp = gpar(fill = fill, col = COL$white, lwd = 0.9)
  )
}

draw_figure4 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = COL$white, col = NA))

  grid.text("a", x = 0.035, y = 0.965, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 12, fontface = "bold", col = COL$ink))
  grid.text("ABF classifications across 6,434 trait-specific analysis-window records",
            x = 0.057, y = 0.965, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 11.2, fontface = "bold", col = COL$ink))

  bar_left <- 0.055
  bar_width <- 0.89
  bar_y <- 0.865
  bar_height <- 0.050
  cumulative <- 0
  for (id in class_order) {
    width <- bar_width * class_expected[[id]] / 6434
    grid.rect(x = bar_left + cumulative + width / 2, y = bar_y, width = width,
              height = bar_height,
              gp = gpar(fill = class_colours[[id]], col = COL$white, lwd = 0.45))
    cumulative <- cumulative + width
  }
  grid.rect(x = bar_left + bar_width / 2, y = bar_y, width = bar_width,
            height = bar_height, gp = gpar(fill = NA, col = COL$navy, lwd = 0.7))

  support_n <- class_expected[["robust_coloc"]] + class_expected[["prior_sensitive_coloc"]]
  support_mid <- bar_left + bar_width * support_n / 6434 / 2
  grid.lines(x = unit(c(support_mid, support_mid, 0.235), "npc"),
             y = unit(c(bar_y + bar_height / 2, 0.922, 0.922), "npc"),
             gp = gpar(col = COL$navy, lwd = 0.8))
  grid.text(
    paste0(fmt_n(support_n), " / 6,434 (", fmt_pct(support_n / 6434),
           ") entered integrated review"),
    x = 0.242, y = 0.922, just = "left",
    gp = gpar(fontfamily = FONT, fontsize = 8.2, fontface = "bold", col = COL$navy)
  )

  legend_item <- function(x, y, id) {
    grid.rect(x = x, y = y, width = 0.012, height = 0.012,
              gp = gpar(fill = class_colours[[id]], col = NA))
    label <- paste0(class_labels[[id]], "  ", fmt_n(class_expected[[id]]),
                    " (", fmt_pct(class_expected[[id]] / 6434), ")")
    grid.text(label, x = x + 0.011, y = y, just = "left",
              gp = gpar(fontfamily = FONT, fontsize = 7.6, col = COL$ink))
  }
  legend_item(0.060, 0.800, "robust_coloc")
  legend_item(0.330, 0.800, "prior_sensitive_coloc")
  legend_item(0.650, 0.800, "distinct_signal")
  legend_item(0.060, 0.758, "trait_specific_or_low_power")
  legend_item(0.520, 0.758, "inconclusive")

  grid.text("b", x = 0.035, y = 0.700, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 12, fontface = "bold", col = COL$ink))
  grid.text("Integrated disposition of records with posterior support favoring H4",
            x = 0.057, y = 0.700, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 11.2, fontface = "bold", col = COL$ink))
  grid.text("The 101-record subset is magnified; ribbon widths are proportional within this panel.",
            x = 0.057, y = 0.670, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 7.3, col = COL$grey))
  grid.text("ABF classification", x = 0.262, y = 0.630,
            gp = gpar(fontfamily = FONT, fontsize = 7.0, fontface = "bold", col = COL$grey))
  grid.text("Integrated disposition", x = 0.708, y = 0.630,
            gp = gpar(fontfamily = FONT, fontsize = 7.0, fontface = "bold", col = COL$grey))

  scale_y <- 0.25 / 101
  source_gap <- 0.014
  destination_gap <- 0.010
  source_bounds <- list(
    robust_coloc = c(lo = 0.600 - 14 * scale_y, hi = 0.600),
    prior_sensitive_coloc = c(
      lo = 0.600 - 14 * scale_y - source_gap - 87 * scale_y,
      hi = 0.600 - 14 * scale_y - source_gap
    )
  )
  destination_bounds <- list()
  destination_top <- 0.600
  for (id in tier_order) {
    destination_bounds[[id]] <- c(
      lo = destination_top - tier_expected[[id]] * scale_y,
      hi = destination_top
    )
    destination_top <- destination_bounds[[id]][["lo"]] - destination_gap
  }

  flows <- copy(flow_expected)
  flows[, `:=`(source_lo = NA_real_, source_hi = NA_real_,
               destination_lo = NA_real_, destination_hi = NA_real_)]
  source_cursor <- vapply(source_bounds, function(z) z[["hi"]], numeric(1))
  for (i in seq_len(nrow(flows))) {
    id <- flows$abf_class[[i]]
    flows$source_hi[[i]] <- source_cursor[[id]]
    flows$source_lo[[i]] <- source_cursor[[id]] - flows$n[[i]] * scale_y
    source_cursor[[id]] <- flows$source_lo[[i]]
  }
  destination_cursor <- vapply(destination_bounds, function(z) z[["hi"]], numeric(1))
  for (destination in tier_order) {
    for (source in c("robust_coloc", "prior_sensitive_coloc")) {
      i <- which(flows$evidence_tier == destination & flows$abf_class == source)
      if (!length(i)) next
      flows$destination_hi[[i]] <- destination_cursor[[destination]]
      flows$destination_lo[[i]] <- destination_cursor[[destination]] - flows$n[[i]] * scale_y
      destination_cursor[[destination]] <- flows$destination_lo[[i]]
    }
  }
  for (i in order(flows$n, decreasing = TRUE)) {
    draw_ribbon(
      0.275, 0.695, flows$source_lo[[i]], flows$source_hi[[i]],
      flows$destination_lo[[i]], flows$destination_hi[[i]],
      class_colours[[flows$abf_class[[i]]]]
    )
  }
  for (id in names(source_bounds)) {
    bounds <- source_bounds[[id]]
    draw_node(0.250, 0.275, bounds[["lo"]], bounds[["hi"]], class_colours[[id]])
    grid.text(paste0(class_labels[[id]], "\n", fmt_n(class_expected[[id]])),
              x = 0.235, y = mean(bounds), just = "right",
              gp = gpar(fontfamily = FONT, fontsize = 7.4, fontface = "bold",
                        col = COL$ink, lineheight = 1.05))
  }

  tier_text <- c(
    abf_label_only_complex_downgraded = paste(
      "Complex-region downgraded", "55 total", "4 robust + 51 prior-sensitive", sep = "\n"
    ),
    abf_label_only_forward_unassessable = paste(
      "Regional forward direction", "unassessable (35 total)",
      "10 robust + 25 prior-sensitive", sep = "\n"
    ),
    abf_prior_sensitive_forward_consistent = paste(
      "Prior-sensitive and", "forward-consistent (10)", "all prior-sensitive", sep = "\n"
    ),
    abf_label_only_forward_opposite = paste(
      "Opposite regional direction", "1 prior-sensitive record", sep = "\n"
    )
  )
  tier_label_y <- vapply(destination_bounds, mean, numeric(1))
  tier_label_y[["abf_prior_sensitive_forward_consistent"]] <- 0.350
  tier_label_y[["abf_label_only_forward_opposite"]] <- 0.300
  for (id in tier_order) {
    bounds <- destination_bounds[[id]]
    draw_node(0.695, 0.720, bounds[["lo"]], bounds[["hi"]], tier_colours[[id]])
    grid.lines(x = unit(c(0.720, 0.730, 0.735), "npc"),
               y = unit(c(mean(bounds), mean(bounds), tier_label_y[[id]]), "npc"),
               gp = gpar(col = tier_colours[[id]], lwd = 0.7))
    grid.text(tier_text[[id]], x = 0.740, y = tier_label_y[[id]], just = "left",
              gp = gpar(fontfamily = FONT, fontsize = 6.8, col = COL$ink, lineheight = 1.03))
  }

  grid.text("c", x = 0.035, y = 0.260, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 12, fontface = "bold", col = COL$ink))
  grid.text("Composite mechanism-support assessment", x = 0.057, y = 0.260, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 11.2, fontface = "bold", col = COL$ink))
  grid.roundrect(x = 0.50, y = 0.155, width = 0.89, height = 0.155,
                 r = unit(2.2, "mm"),
                 gp = gpar(fill = "#F7F9FA", col = "#9CB2C8", lwd = 1.1))
  grid.text("0 / 101", x = 0.160, y = 0.166,
            gp = gpar(fontfamily = FONT, fontsize = 28, fontface = "bold", col = COL$navy))
  grid.text("No record met every documented criterion", x = 0.285, y = 0.190, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 10.2, fontface = "bold", col = COL$ink))
  grid.text(
    paste(
      "Required jointly: robust ABF classification; a non-complex region outside predefined",
      "pleiotropic regions; and assessable, directionally concordant forward and reverse",
      "regional IVW sensitivity estimates.", sep = "\n"
    ),
    x = 0.285, y = 0.135, just = "left",
    gp = gpar(fontfamily = FONT, fontsize = 7.6, col = COL$ink, lineheight = 1.12)
  )
  grid.text("ABF, approximate Bayes factor; H4, shared-causal-variant hypothesis.",
            x = 0.055, y = 0.052, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 7.2, col = COL$grey))
  grid.text("Counts are trait-specific analysis-window records, not independent physical loci.",
            x = 0.055, y = 0.032, just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 7.2, col = COL$grey))
  grid.text(
    "ABF classifications are model- and prior-dependent and do not establish mediation or a trait-specific causal mechanism.",
    x = 0.055, y = 0.014, just = "left",
    gp = gpar(fontfamily = FONT, fontsize = 7.2, col = COL$grey)
  )
}

figure4_width <- 7.2
figure4_height <- 6.6
cairo_pdf(file.path(OUT, "Figure_4.pdf"), width = figure4_width,
          height = figure4_height, bg = COL$white)
draw_figure4()
dev.off()
png(file.path(OUT, "Figure_4_600dpi.png"), width = figure4_width,
    height = figure4_height, units = "in", res = 600, bg = COL$white)
draw_figure4()
dev.off()
svg(file.path(OUT, "Figure_4.svg"), width = figure4_width,
    height = figure4_height, family = FONT, bg = COL$white)
draw_figure4()
dev.off()

cat("Scientific Reports Figures 1-4 regenerated from frozen sources\n")
cat(list.files(OUT, pattern = "^Figure_[1-4]", full.names = FALSE), sep = "\n")


expected <- c(
  file.path(OUT, paste0("Figure_", rep(1:4, each = 3),
                       rep(c(".pdf", ".svg", "_600dpi.png"), times = 4))),
  file.path(OUT, "Figure_1_source.mermaid"),
  file.path(OUT, "Figure_2_3_trait_labels.tsv"),
  file.path(OUT, "Figure_2_plot_data.tsv"),
  file.path(OUT, "Figure_3_matrix_data.tsv"),
  file.path(OUT, "Figure_4_plot_data.tsv")
)
expected_sizes <- file.info(expected)$size
small_text_output <- grepl("\\.(tsv|mermaid)$", expected)
if (!all(file.exists(expected)) ||
    any(expected_sizes[small_text_output] < 100) ||
    any(expected_sizes[!small_text_output] < 1000)) {
  stop("One or more expected figure outputs are absent or unexpectedly small")
}
checksums <- data.table(
  file = basename(expected),
  bytes = file.info(expected)$size,
  md5 = unname(tools::md5sum(expected))
)
fwrite(checksums, file.path(OUT, "figure_checksums.tsv"), sep = "\t", eol = "\n")
cat("Reproduction checks passed; checksum manifest written to the selected output directory\n")
