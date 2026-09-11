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
full_dictionary <- fread(file.path(DATA, "trait_display_dictionary.tsv"))
if (nrow(dictionary) != 15L || uniqueN(dictionary$trait_id) != 15L ||
    uniqueN(dictionary$trait_display) != 15L || any(grepl("_", dictionary$trait_display, fixed = TRUE))) {
  stop("Trait display dictionary failed integrity checks")
}
if (nrow(full_dictionary) != 249L || uniqueN(full_dictionary$trait_id) != 249L ||
    uniqueN(full_dictionary$trait_display) != 249L || any(grepl("_", full_dictionary$trait_display, fixed = TRUE))) {
  stop("The full trait display dictionary failed integrity checks")
}
setkey(dictionary, trait_id)
setkey(full_dictionary, trait_id)
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
  header_card <- function(x, y, w, h, title, body, fill, border, body_size = 7.2,
                          body_lineheight = 1.17) {
    rounded_box(x, y, w, h, fill = fig$white, border = border, lwd = 1.25)
    header_h <- 0.041
    grid.roundrect(
      x = x, y = y + (h - header_h) / 2, width = w, height = header_h,
      r = unit(2.2, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.25)
    )
    grid.rect(x = x, y = y + (h - header_h) / 2 - 0.0065,
              width = w - 0.003, height = 0.013, gp = gpar(fill = fill, col = NA))
    text_at(title, x, y + (h - header_h) / 2, size = 8.8, face = "bold", colour = fig$navy)
    text_at(body, x, y - header_h * 0.32, size = body_size,
            lineheight = body_lineheight)
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
      rounded_box(x, y - 0.062, w - 0.026, 0.058, fill = fig$pale, border = fig$line,
                  lwd = 0.85, radius = 1.5)
      text_at(footer, x, y - 0.062, size = 6.4, lineheight = 1.04)
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
  stage_heading("Candidate-level evidence assessment", 0.558)
  stage_heading("Evidence synthesis", 0.285)
  stage_heading("Alternative outcome", 0.105)

  # The two MR directions have separate outputs. Only the primary forward
  # screen selects the candidate set used by all three downstream assessments.
  connector(0.29, 0.825, 0.40, 0.785)
  connector(0.76, 0.825, 0.66, 0.785)
  connector(0.35, 0.647, 0.35, 0.642)
  connector(0.72, 0.647, 0.72, 0.642)
  connector(0.35, 0.604, 0.35, 0.530, end = FALSE)
  connector(0.205, 0.530, 0.760, 0.530, end = FALSE)
  connector(0.205, 0.530, 0.205, 0.515)
  connector(0.53, 0.530, 0.53, 0.515)
  connector(0.760, 0.530, 0.760, 0.515)
  connector(0.205, 0.315, 0.205, 0.300, end = FALSE)
  connector(0.53, 0.315, 0.53, 0.300, end = FALSE)
  connector(0.855, 0.315, 0.855, 0.300, end = FALSE)
  connector(0.205, 0.300, 0.855, 0.300, end = FALSE)
  connector(0.53, 0.300, 0.53, 0.263)
  connector(0.53, 0.153, 0.53, 0.086, dashed = TRUE)

  header_card(
    0.285, 0.880, 0.405, 0.105, "Circulating metabolic-trait GWAS",
    paste0(fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits\n",
           "Estonian Biobank and UK Biobank   N = ", fmt("tambets_meta_eur")),
    fig$yellow, fig$yellow_line
  )
  header_card(
    0.755, 0.880, 0.405, 0.105, "PGC major-depression GWAS",
    paste0("Primary: ", fmt("pgc_main_cases"), " cases   ",
           fmt("pgc_main_controls"), " controls\n",
           "UK Biobank-excluded: ", fmt("pgc_nb_cases"), " cases\n",
           fmt("pgc_nb_controls"), " controls"),
    fig$yellow, fig$yellow_line, body_size = 6.35, body_lineheight = 1.04
  )

  rounded_box(0.53, 0.690, 0.77, 0.180, fill = fig$white,
              border = fig$blue_line, lwd = 1.3)
  grid.roundrect(x = 0.53, y = 0.751, width = 0.77, height = 0.034,
                 r = unit(2.0, "mm"), gp = gpar(fill = fig$blue, col = fig$blue_line, lwd = 1.3))
  grid.rect(x = 0.53, y = 0.745, width = 0.767, height = 0.011,
            gp = gpar(fill = fig$blue, col = NA))
  text_at("Parallel direction-specific MR analyses",
          0.53, 0.751, size = 8.7, face = "bold", colour = fig$navy)
  text_at("Exposure-specific instruments selected separately in each direction",
          0.53, 0.727, size = 5.8, colour = fig$slate)

  rounded_box(0.35, 0.685, 0.34, 0.076, fill = fig$muted_blue,
              border = fig$blue_line, lwd = 0.9, radius = 1.5)
  text_at("Forward MR", 0.35, 0.706, size = 7.2,
          face = "bold", colour = fig$navy)
  text_at(
    paste0(fmt("n_traits_metabolites"), " metabolic traits as exposures\n",
           "PGC major depression as outcome\n",
           "Primary and UK Biobank-excluded screens"),
    0.35, 0.674, size = 6.5, lineheight = 1.05
  )

  rounded_box(0.72, 0.685, 0.34, 0.076, fill = "#F0EDF6",
              border = fig$violet_line, lwd = 0.9, radius = 1.5)
  text_at("Reverse MR", 0.72, 0.706, size = 7.2,
          face = "bold", colour = fig$navy)
  text_at(
    paste0("Genetic liability to MDD as exposure\n",
           fmt("n_traits_metabolites"), " metabolic traits as outcomes\n",
           "Primary and UK Biobank-excluded screens"),
    0.72, 0.674, size = 6.5, lineheight = 1.05
  )

  rounded_box(0.35, 0.623, 0.30, 0.040, fill = fig$blue,
              border = fig$blue_line, lwd = 1.1, radius = 1.5)
  text_at(paste0(fmt("n_screen_bonf"), " forward-selected traits"),
          0.35, 0.630, size = 6.25, face = "bold", colour = fig$navy)
  text_at("Primary forward screen only",
          0.35, 0.614, size = 5.5, colour = fig$slate)

  rounded_box(0.72, 0.623, 0.30, 0.040, fill = fig$violet,
              border = fig$violet_line, lwd = 1.1, radius = 1.5)
  text_at("Reverse MR for 15 forward-selected traits",
          0.72, 0.630, size = 5.55, face = "bold", colour = fig$navy)
  text_at("Both complete 249-trait screens",
          0.72, 0.614, size = 5.5, colour = fig$slate)

  # Reverse MR informs only the directional assessment. Its path is routed to
  # the right of the forward-candidate branch so the connectors do not cross.
  reverse_path_x <- c(0.72, 0.72, 0.90, 0.90)
  reverse_path_y <- c(0.603, 0.548, 0.548, 0.515)
  grid.lines(
    x = reverse_path_x, y = reverse_path_y,
    gp = gpar(col = fig$white, lwd = 3.6)
  )
  grid.lines(
    x = reverse_path_x, y = reverse_path_y,
    gp = gpar(col = fig$violet_line, lwd = 1.25),
    arrow = arrow(type = "closed", length = unit(1.7, "mm"))
  )

  evaluation_card(
    0.205, "Robustness and sensitivity",
    c("Same 15 forward candidates", "Five-estimator comparison",
      "Cochran's Q, MR-Egger and MR-PRESSO", "Alternative harmonization",
      "Regional and shared-instrument exclusions"),
    fig$teal, fig$teal_line,
    paste0("UK Biobank-excluded PGC sensitivity\n", fmt("pgc_nb_cases"), " cases   ",
           fmt("pgc_nb_controls"), " controls\nSame 15 summarized; full screen reported")
  )
  evaluation_card(
    0.53, "Locus evidence",
    c("Same 15 forward candidates", "Analysis-window construction",
      "ABF colocalization in both PGC outcomes", "Regional direction checks",
      "MHC and complex-region assessment"),
    fig$peach, fig$peach_line,
    paste0(fmt("coloc_records"), " trait-specific\nanalysis-window records")
  )
  evaluation_card(
    0.855, "Directional evidence assessment",
    c("15 forward-selected traits",
      "Steiger directionality test",
      "under three MDD prevalence assumptions",
      "Matched reverse-MR estimates",
      "from both complete 249-trait screens"),
    fig$violet, fig$violet_line
  )

  rounded_box(0.53, 0.205, 0.77, 0.105, fill = fig$white,
              border = fig$synthesis_line, lwd = 1.3)
  grid.roundrect(x = 0.53, y = 0.244, width = 0.77, height = 0.032,
                 r = unit(2.0, "mm"), gp = gpar(fill = fig$synthesis,
                                                col = fig$synthesis_line, lwd = 1.3))
  grid.rect(x = 0.53, y = 0.239, width = 0.767, height = 0.010,
            gp = gpar(fill = fig$synthesis, col = NA))
  text_at("Integrated evidence assessment", 0.53, 0.244, size = 8.5,
          face = "bold", colour = fig$navy)
  grid.lines(x = c(0.53, 0.53), y = c(0.163, 0.230), gp = gpar(col = fig$line, lwd = 1.0))
  text_at("Trait-level reporting", 0.345, 0.222, size = 7.2, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("groups_priority"), " priority-reporting traits\n",
                 fmt("groups_secondary"), " secondary-reporting traits"),
          0.345, 0.191, size = 6.8, face = "bold", colour = fig$navy, lineheight = 1.10)
  text_at("Locus-level assessment", 0.715, 0.222, size = 7.2, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("h4_evidence_rows"), " with shared-variant posterior support"),
          0.715, 0.198, size = 6.4, face = "bold", colour = fig$navy)
  text_at(paste0(fmt("mechanism_eligible_true"), "/", fmt("h4_evidence_rows"),
                 " records met the composite\nmechanism-support criterion"),
          0.715, 0.169, size = 6.3, colour = fig$slate, lineheight = 1.06)

  rounded_box(0.53, 0.050, 0.70, 0.075, fill = fig$white,
              border = fig$green_line, lwd = 1.25, lty = 2)
  grid.roundrect(x = 0.53, y = 0.078, width = 0.70, height = 0.026,
                 r = unit(1.8, "mm"), gp = gpar(fill = fig$green,
                                                col = fig$green_line, lwd = 1.2, lty = 2))
  grid.rect(x = 0.53, y = 0.074, width = 0.697, height = 0.008,
            gp = gpar(fill = fig$green, col = NA))
  text_at("Separate FinnGen R13 alternative-outcome comparison",
          0.53, 0.078, size = 7.4, face = "bold", colour = fig$navy)
  text_at(
    paste0(fmt("finngen_cases"), " cases   ", fmt("finngen_controls"), " controls\n",
           fmt("fg13_dir_same"), "/", fmt("matrix_fg13_rows"), " directionally concordant\n",
           support_counts[["FDR-supported"]], " FDR-supported   ",
           support_counts[["nominally-supported"]], " nominally supported   ",
           support_counts[["no-nominal-support"]], " without nominal support\n",
           "Reporting groups unchanged; not an independent replication"),
    0.53, 0.041, size = 6.2, lineheight = 1.03
  )

  popViewport()
}

mermaid <- c(
  "flowchart TB",
  paste0("  E[\"Circulating metabolic-trait GWAS<br/>", fmt("n_traits_metabolites"), " NMR-derived circulating metabolic traits<br/>N = ", fmt("tambets_meta_eur"), "; Estonian Biobank and UK Biobank\"]"),
  paste0("  O[\"PGC major-depression GWAS<br/>Primary: ", fmt("pgc_main_cases"), " cases / ", fmt("pgc_main_controls"), " controls<br/>UK Biobank-excluded: ", fmt("pgc_nb_cases"), " cases / ", fmt("pgc_nb_controls"), " controls\"]"),
  "  FW[\"Forward MR<br/>249 metabolic traits as exposures; PGC major depression as outcome<br/>Primary and UK Biobank-excluded screens\"]",
  "  RV[\"Reverse MR<br/>Genetic liability to MDD as exposure; 249 metabolic traits as outcomes<br/>Primary and UK Biobank-excluded screens\"]",
  paste0("  C[\"", fmt("n_screen_bonf"), " forward-selected traits<br/>Primary forward screen only\"]"),
  "  RM[\"Reverse MR for 15 forward-selected traits<br/>Both complete 249-trait screens\"]",
  "  R[\"Robustness and sensitivity<br/>Same 15 forward candidates; five estimators; Cochran's Q, MR-Egger and MR-PRESSO<br/>UK Biobank-excluded sensitivity; regional and shared-instrument exclusions\"]",
  paste0("  L[\"Locus evidence<br/>Same 15 forward candidates; ABF colocalization and regional direction checks<br/>MHC and complex-region assessment; ", fmt("coloc_records"), " trait-specific analysis-window records\"]"),
  "  D[\"Directional evidence assessment<br/>15 forward-selected traits; Steiger directionality test at three MDD prevalences<br/>Matched reverse-MR estimates from both complete 249-trait screens\"]",
  paste0("  G[\"Integrated evidence assessment<br/>", fmt("groups_priority"), " priority-reporting / ", fmt("groups_secondary"), " secondary-reporting traits<br/>", fmt("h4_evidence_rows"), " with shared-variant posterior support; ", fmt("mechanism_eligible_true"), "/", fmt("h4_evidence_rows"), " met the composite criterion\"]"),
  paste0("  F[\"Separate FinnGen R13 alternative-outcome comparison<br/>", fmt("fg13_dir_same"), "/", fmt("matrix_fg13_rows"), " directionally concordant; reporting groups unchanged<br/>Not an independent replication\"]"),
  "  E --> FW",
  "  O --> FW",
  "  E --> RV",
  "  O --> RV",
  "  FW --> C",
  "  RV --> RM",
  "  C --> R",
  "  C --> L",
  "  C --> D",
  "  RM --> D",
  "  R --> G",
  "  L --> G",
  "  D --> G",
  "  G -.-> F"
)
stopifnot(
  all(c("  E --> FW", "  O --> FW", "  E --> RV", "  O --> RV",
        "  FW --> C", "  RV --> RM", "  C --> R", "  C --> L",
        "  C --> D", "  RM --> D", "  R --> G", "  L --> G",
        "  D --> G", "  G -.-> F") %in% mermaid),
  !any(c("  RV --> C", "  RM --> R", "  RM --> L", "  F --> G") %in% mermaid)
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

# Figure 2: full-screen overview followed by three-series candidate estimates.
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

screen_level_order <- c("bonferroni_hit", "fdr_only", "nominal", "null")
screen_level_labels <- c(
  bonferroni_hit = "Bonferroni-significant (15)",
  fdr_only = "FDR-only (53)",
  nominal = "Nominal only (35)",
  null = "No nominal association (146)"
)
screen_level_colours <- c(
  bonferroni_hit = COL$navy,
  fdr_only = COL$teal,
  nominal = COL$orange,
  null = "#B8C1C7"
)
screen_level_shapes <- c(bonferroni_hit = 16, fdr_only = 17, nominal = 15, null = 1)
screen_plot <- copy(forward_screen)
screen_plot[, trait_display := full_dictionary[J(trait), trait_display]]
if (nrow(screen_plot) != 249L || uniqueN(screen_plot$trait) != 249L ||
    anyNA(screen_plot$trait_display) ||
    !identical(screen_plot[, .N, by = screen_level][match(screen_level_order, screen_level), N],
               c(15L, 53L, 35L, 146L))) {
  stop("Figure 2 full-screen source failed the 249-trait and significance-class checks")
}
screen_plot[, screen_level := factor(screen_level, levels = screen_level_order)]
fwrite(
  screen_plot[, .(
    trait_id = trait, trait_display, n_iv, beta = ivw_b, se = ivw_se, p_value = ivw_p,
    estimate = ivw_or, lower = ivw_or_lci, upper = ivw_or_uci,
    p_fdr_bh, p_bonf, screen_level = as.character(screen_level)
  )],
  file.path(OUT, "Figure_2_screen_data.tsv"), sep = "\t", eol = "\n"
)

p2_overview <- ggplot(screen_plot, aes(ivw_b, -log10(ivw_p), colour = screen_level,
                                       shape = screen_level)) +
  geom_hline(yintercept = -log10(0.05), colour = "#AAB5BC", linewidth = 0.35,
             linetype = 3) +
  geom_hline(yintercept = -log10(0.05 / 249), colour = COL$navy, linewidth = 0.45,
             linetype = 2) +
  geom_vline(xintercept = 0, colour = "#AAB5BC", linewidth = 0.35) +
  geom_point(size = 2.1, stroke = 0.55, alpha = 0.92) +
  scale_colour_manual(values = screen_level_colours, labels = screen_level_labels,
                      breaks = screen_level_order) +
  scale_shape_manual(values = screen_level_shapes, labels = screen_level_labels,
                     breaks = screen_level_order) +
  scale_x_continuous(breaks = seq(-0.08, 0.08, 0.04)) +
  labs(
    tag = "a",
    x = "Primary PGC log odds ratio per genetically predicted 1-SD increase",
    y = expression(-log[10](italic(P))),
    subtitle = "All 249 circulating metabolic traits; dashed line indicates P = 0.05/249"
  ) +
  theme_sr(9.5) +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.box = "horizontal", plot.subtitle = element_text(colour = COL$grey, size = 8.3)
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
  labs(tag = "b", x = "Odds ratio per genetically predicted 1-SD increase (95% CI)", y = NULL) +
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
p2_candidate <- p2_forest + p2_values + plot_layout(widths = c(0.69, 0.31), guides = "collect") &
  theme(legend.position = "top")
p2 <- p2_overview / p2_candidate +
  plot_layout(heights = c(0.32, 0.68)) &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = 12))
save_gg(p2, "Figure_2", 13.0, 12.0)

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
m3[, domain := fcase(
  layer %chin% c("Effect\nIVW dir.", "Effect\n5 methods"), "Association direction",
  layer %chin% c("Robustness\nPalindromic exclusion", "Robustness\nUK Biobank excluded",
                 "Robustness\nMR-PRESSO", "Sensitivity\nShared-IV removal"),
    "Robustness and sensitivity",
  layer %chin% c("Diagnostics\nEgger int.", "Diagnostics\nCochran Q"), "Diagnostics",
  layer == "Direction\nReverse MR", "Directionality",
  layer %chin% c("Alternative outcome\nFG direction", "Alternative outcome\nP level"),
    "Alternative outcome",
  layer == "Locus\nColoc", "Locus evidence"
)]
domain_order <- c(
  "Association direction", "Robustness and sensitivity", "Diagnostics",
  "Directionality", "Alternative outcome", "Locus evidence"
)
layer_labels <- c(
  "Effect\nIVW dir." = "IVW\ndirection",
  "Effect\n5 methods" = "Five-method\nagreement",
  "Robustness\nPalindromic exclusion" = "Palindromic\nexclusion",
  "Robustness\nUK Biobank excluded" = "UK Biobank\nexcluded",
  "Robustness\nMR-PRESSO" = "MR-PRESSO",
  "Sensitivity\nShared-IV removal" = "Shared-instrument\nremoval",
  "Diagnostics\nEgger int." = "MR-Egger\nintercept",
  "Diagnostics\nCochran Q" = "Cochran's Q",
  "Direction\nReverse MR" = "Reverse MR",
  "Alternative outcome\nFG direction" = "FinnGen\ndirection",
  "Alternative outcome\nP level" = "FinnGen\nsignificance",
  "Locus\nColoc" = "Colocalization"
)
m3[, trait_display := wrap_trait(figure_full_name(trait))]
m3[, group_label := fifelse(
  t2$reporting_group[match(trait, t2$trait)] == "priority_reporting",
  "Priority group", "Secondary group"
)]
m3[, trait_display := factor(trait_display, levels = rev(t2$trait_figure_label))]
m3[, layer := factor(layer, levels = layer_order)]
m3[, domain := factor(domain, levels = domain_order)]
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
    domain = as.character(domain),
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
  facet_grid(group_label ~ domain, scales = "free", space = "free") +
  scale_fill_manual(values = c(support = "#DCEBDD", qualified = "#F6E8C9",
                               limitation = "#F2D7D5", positive = "#D9E8F2",
                               negative = "#E9E2F2", neutral = "#EDF0F2"), guide = "none") +
  scale_x_discrete(labels = layer_labels) +
  labs(x = NULL, y = NULL) +
  theme_sr(8.6) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.7),
    panel.spacing = unit(1.5, "mm"),
    strip.text.x = element_text(hjust = 0.5, size = 8.0),
    strip.text.y = element_text(angle = 0, size = 7.8)
  )
save_gg(p3, "Figure_3", 13.0, 8.8)

# Figure 4: aggregate, trait-level and integrated colocalization summaries.
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

support_classes <- c("robust_coloc", "prior_sensitive_coloc")
trait_support <- data.table(
  trait = rep(t2$trait, each = length(support_classes)),
  abf_class = rep(support_classes, times = nrow(t2)),
  n = 0L
)
trait_support[
  classification[abf_class %chin% support_classes, .N, by = .(trait, abf_class)],
  on = .(trait, abf_class), n := i.N
]
robust_by_trait <- trait_support[abf_class == "robust_coloc", n]
prior_by_trait <- trait_support[abf_class == "prior_sensitive_coloc", n]
if (!identical(robust_by_trait, as.integer(t2$coloc_n_robust)) ||
    !identical(prior_by_trait, as.integer(t2$coloc_n_prior)) ||
    sum(trait_support$n) != 101L) {
  stop("Trait-level shared-variant-support counts do not match the released evidence summary")
}
trait_support[, `:=`(
  trait_display = t2$trait_figure_label[match(trait, t2$trait)],
  reporting_group = t2$reporting_group[match(trait, t2$trait)],
  group_label = fifelse(
    t2$reporting_group[match(trait, t2$trait)] == "priority_reporting",
    "Priority group", "Secondary group"
  )
)]

figure4_plot_data <- rbindlist(list(
  data.table(
    panel = "all_records", source = class_order, destination = NA_character_,
    n = as.integer(class_values), denominator = 6434L,
    percent = 100 * as.integer(class_values) / 6434
  ),
  trait_support[, .(
    panel = "trait_support", source = trait, destination = abf_class,
    n, denominator = 101L, percent = 100 * n / 101,
    trait_display = figure_full_name(trait), reporting_group
  )],
  data.table(
    panel = "integrated_disposition", source = tier_order, destination = NA_character_,
    n = as.integer(tier_values), denominator = 101L,
    percent = 100 * as.integer(tier_values) / 101
  ),
  data.table(
    panel = "composite_assessment", source = "mechanism_eligible",
    destination = NA_character_, n = 0L, denominator = 101L, percent = 0
  )
), fill = TRUE)
fwrite(figure4_plot_data, file.path(OUT, "Figure_4_plot_data.tsv"), sep = "\t", eol = "\n")

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

all_records_plot <- data.table(
  abf_class = factor(class_order, levels = class_order),
  n = as.integer(class_values)
)
all_records_legend <- setNames(
  paste0(class_labels[class_order], " (", fmt_n(class_expected[class_order]), ")"),
  class_order
)
p4a <- ggplot(all_records_plot, aes(n, "All records", fill = abf_class)) +
  geom_col(width = 0.58, colour = COL$white, linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  geom_text(
    data = all_records_plot[n >= 500],
    aes(label = fmt_n(n)), position = position_stack(vjust = 0.5, reverse = TRUE),
    family = FONT, size = 2.8, colour = COL$white, fontface = "bold"
  ) +
  scale_fill_manual(values = class_colours, labels = all_records_legend, breaks = class_order) +
  scale_x_continuous(labels = function(x) fmt_n(x), expand = expansion(mult = c(0, 0.015))) +
  labs(
    tag = "a",
    title = "ABF classifications across 6,434 trait-specific analysis-window records",
    subtitle = "The 101 robust or prior-sensitive records (1.57%) entered integrated review",
    x = "Analysis-window records", y = NULL
  ) +
  theme_sr(9.2) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom", legend.justification = "left",
    legend.text = element_text(size = 7.5),
    plot.subtitle = element_text(colour = COL$grey, size = 8.1),
    axis.text.y = element_text(face = "bold")
  )

trait_support[, abf_class := factor(abf_class, levels = support_classes)]
trait_support[, trait_display := factor(trait_display, levels = rev(t2$trait_figure_label))]
support_colours <- c(robust_coloc = class_colours[["robust_coloc"]],
                     prior_sensitive_coloc = class_colours[["prior_sensitive_coloc"]])
support_labels <- c(robust_coloc = "Robust ABF", prior_sensitive_coloc = "Prior-sensitive ABF")
p4b <- ggplot(trait_support, aes(n, trait_display, fill = abf_class)) +
  geom_col(width = 0.68, colour = COL$white, linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  geom_text(
    data = trait_support[n > 0], aes(label = n, colour = abf_class),
    position = position_stack(vjust = 0.5, reverse = TRUE), family = FONT, size = 2.65,
    fontface = "bold", show.legend = FALSE
  ) +
  facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = support_colours, labels = support_labels, breaks = support_classes) +
  scale_colour_manual(values = c(robust_coloc = COL$white, prior_sensitive_coloc = COL$ink)) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, 0.01))) +
  labs(
    tag = "b",
    title = "Distribution of the 101 shared-variant-support records across the 15 traits",
    x = "Records per trait", y = NULL
  ) +
  theme_sr(9.0) +
  theme(
    legend.position = "top", legend.justification = "left",
    axis.text.y = element_text(size = 7.6), panel.spacing.y = unit(2.5, "mm")
  )

tier_labels <- c(
  abf_label_only_complex_downgraded = "Complex-region downgraded",
  abf_label_only_forward_unassessable = "Forward direction unassessable",
  abf_prior_sensitive_forward_consistent = "Prior-sensitive and forward-consistent",
  abf_label_only_forward_opposite = "Opposite forward direction"
)
disposition_plot <- data.table(
  evidence_tier = factor(tier_order, levels = rev(tier_order)),
  n = as.integer(tier_values)
)
p4c <- ggplot(disposition_plot, aes(n, evidence_tier, fill = evidence_tier)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), hjust = -0.25, family = FONT, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = tier_colours, guide = "none") +
  scale_y_discrete(labels = tier_labels) +
  scale_x_continuous(limits = c(0, 62), expand = expansion(mult = c(0, 0))) +
  labs(tag = "c", title = "Integrated disposition of the 101-record subset",
       x = "Records", y = NULL) +
  theme_sr(8.8) +
  theme(axis.text.y = element_text(size = 7.5), panel.grid.major.y = element_blank())

p4d <- ggplot() +
  annotate("text", x = 0.5, y = 0.70, label = "0 / 101", family = FONT,
           size = 11, fontface = "bold", colour = COL$navy) +
  annotate("text", x = 0.5, y = 0.45,
           label = "met every documented\nmechanism-support criterion",
           family = FONT, size = 3.6, fontface = "bold", colour = COL$ink) +
  annotate("text", x = 0.5, y = 0.18,
           label = "The result is specific to the stated criteria\nand available summary data.",
           family = FONT, size = 2.8, colour = COL$grey) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT) +
  theme(plot.background = element_rect(fill = "#F7F9FA", colour = "#9CB2C8", linewidth = 0.8),
        plot.margin = margin(9, 9, 9, 9))

p4 <- p4a / p4b / (p4c | p4d) +
  plot_layout(heights = c(0.24, 0.52, 0.24), widths = c(0.68, 0.32)) +
  plot_annotation(
    caption = paste(
      "Counts are trait-specific analysis-window records, not independent physical loci.",
      "ABF classifications are prior- and model-dependent and do not establish mediation or a trait-specific causal mechanism."
    )
  ) &
  theme(
    plot.tag = element_text(family = FONT, face = "bold", size = 12),
    plot.title = element_text(family = FONT, face = "bold", size = 10.5),
    plot.caption = element_text(family = FONT, size = 7.2, colour = COL$grey, hjust = 0)
  )
save_gg(p4, "Figure_4", 13.0, 10.8)

cat("Scientific Reports Figures 1-4 regenerated from frozen sources\n")
cat(list.files(OUT, pattern = "^Figure_[1-4]", full.names = FALSE), sep = "\n")


expected <- c(
  file.path(OUT, paste0("Figure_", rep(1:4, each = 3),
                       rep(c(".pdf", ".svg", "_600dpi.png"), times = 4))),
  file.path(OUT, "Figure_1_source.mermaid"),
  file.path(OUT, "Figure_2_3_trait_labels.tsv"),
  file.path(OUT, "Figure_2_plot_data.tsv"),
  file.path(OUT, "Figure_2_screen_data.tsv"),
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
