## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 4.8,
  dpi = 96,
  fig.retina = 2
)

# Number format for p-values and relative errors: three decimals when
# x >= 0.0005, otherwise one decimal times a power of ten (e.g. 3.3 x 10^-4)
fmt_p <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (v == 0) return("0")
    if (abs(v) >= 5e-4) return(sprintf("%.3f", v))
    e <- floor(log10(abs(v)))
    m <- round(abs(v) / 10^e, 1)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("%s%.1f &times; 10<sup>&minus;%d</sup>", if (v < 0) "&minus;" else "", m, -e)
  }, character(1))
}

## ----data---------------------------------------------------------------------
library(MultinomFourier)

votes <- matrix(
  c(
    # Candidate A
     37,  34,  31,  34,  37,  37,  35,  32,  40,  37,  37,  38,  36, 22, 25,
    # Candidate B
     39,  37,  36,  36,  26,  29,  28,  39,  39,  34,  28,  32,  35, 15, 27,
    # Candidate C
     99, 107, 110, 108, 127, 107, 108, 125,  89, 100, 109, 103,  98, 74, 70,
    # Candidate D
     83,  84,   8,  87,  70, 102,  88,  74,  82,  86,  69,  86,  86, 49, 44,
    # Null and blank votes
     46,  37,  42,  41,  34,  37,  45,  26,  62,  45,  48,  49,  41, 21, 21,
    # Abstention
     96, 101, 173,  94, 106,  88,  96, 104,  88,  98, 109,  92, 104, 57, 50
  ),
  nrow = 6,
  byrow = TRUE
)

rownames(votes) <- c(
  "Candidate A",
  "Candidate B",
  "Candidate C",
  "Candidate D",
  "Null/blank",
  "Abstention"
)

colnames(votes) <- as.character(276:290)

polling_station_sizes <- colSums(votes)
polling_station_sizes

## ----heatmap, fig.cap = "Vote proportions by polling station at Escuela Centenario."----
prop <- sweep(votes, 2, polling_station_sizes, "/")

heat_data <- as.data.frame(as.table(prop))
names(heat_data) <- c("Option", "Station", "Proportion")

# Cell labels: percentage with one decimal, in white on dark cells
heat_data$Label <- sprintf("%.1f", 100 * heat_data$Proportion)
dark_cell <- heat_data$Proportion > min(prop) + 0.55 * diff(range(prop))
heat_data$TextColor <- ifelse(dark_cell, "white", "black")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(
    heat_data,
    ggplot2::aes(x = Station, y = Option, fill = Proportion)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label, color = TextColor),
      size = 2.6
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_distiller(
      palette = "YlOrRd",
      direction = 1,
      labels = function(x) paste0(100 * x, "%")
    ) +
    ggplot2::scale_y_discrete(limits = rev(rownames(votes))) +
    ggplot2::labs(
      x = "Polling station",
      y = NULL,
      fill = "Vote proportion"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )
} else {
  # Rows reversed so that Candidate A appears at the top
  z <- t(prop[rev(rownames(prop)), , drop = FALSE])
  op <- par(mar = c(5, 8, 2, 2))
  image(
    seq_len(nrow(z)), seq_len(ncol(z)), z,
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    axes = FALSE,
    xlab = "Polling station",
    ylab = ""
  )
  axis(1, at = seq_len(nrow(z)), labels = rownames(z), las = 1)
  axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 1)
  lab <- matrix(sprintf("%.1f", 100 * z), nrow = nrow(z))
  dark <- z > min(z) + 0.55 * diff(range(z))
  text(
    rep(seq_len(nrow(z)), ncol(z)), rep(seq_len(ncol(z)), each = nrow(z)),
    lab, cex = 0.6, col = ifelse(dark, "white", "black")
  )
  par(op)
}

## ----loo-probabilities--------------------------------------------------------
loo_p <- sapply(
  seq_len(ncol(votes)),
  function(m) {
    other_counts <- rowSums(votes[, -m, drop = FALSE])
    other_counts / sum(other_counts)
  }
)

rownames(loo_p) <- rownames(votes)
colnames(loo_p) <- colnames(votes)

# Reference distribution for station 278
round(loo_p[, "278"], 4)

## ----expected-278-------------------------------------------------------------
comparison_278 <- data.frame(
  Option = rownames(votes),
  Observed = votes[, "278"],
  Expected = round(polling_station_sizes["278"] * loo_p[, "278"], 1),
  Difference = round(votes[, "278"] - polling_station_sizes["278"] * loo_p[, "278"], 1),
  row.names = NULL
)

knitr::kable(
  comparison_278,
  format = "html",
  escape = FALSE,
  align = "l",
  caption = "Observed vs. expected counts for Polling Station 278."
)

## ----fourier-278--------------------------------------------------------------
fit_278 <- pval_fourier(
  x = votes[, "278"],
  p = loo_p[, "278"],
  rel_eps = 1e-6
)

fit_278

## ----exact-278----------------------------------------------------------------
exact_fit_278 <- pval_exact(
  x = votes[, "278"],
  p = loo_p[, "278"]
)

exact_fit_278$pval

## ----flexible-278-------------------------------------------------------------
pval_flexible(
  x = votes[, "278"],
  p = loo_p[, "278"],
  rel_eps = 1e-6
)$method

## ----screening----------------------------------------------------------------
p_fourier_all <- sapply(colnames(votes), function(m) {
  pval_fourier(votes[, m], loo_p[, m], rel_eps = 1e-6)$pval
})

p_exact_all <- sapply(colnames(votes), function(m) {
  # Station 278 was already computed above
  if (m == "278") return(exact_fit_278$pval)
  pval_exact(votes[, m], loo_p[, m])$pval
})

rel_err_all <- abs(p_fourier_all / p_exact_all - 1)

screening <- data.frame(
  Station = colnames(votes),
  Voters = polling_station_sizes,
  p_Fourier = fmt_p(p_fourier_all),
  p_Exact = fmt_p(p_exact_all),
  Rel_Error = fmt_p(rel_err_all),
  row.names = NULL
)

knitr::kable(
  screening,
  format = "html",
  escape = FALSE,
  align = "l",
  caption = "Leave-one-out p-values for all stations (LLR statistic): Fourier inversion vs. exact computation."
)

## ----screening-plot, fig.cap = "Leave-one-out -log10(p) across polling stations (LLR): Fourier inversion vs. exact computation."----
screening_plot_data <- data.frame(
  Station = rep(as.integer(colnames(votes)), 2),
  Method = rep(c("Fourier (pval_fourier)", "Exact (pval_exact)"), each = ncol(votes)),
  neg_log10_p = -log10(c(p_fourier_all, p_exact_all))
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(
    screening_plot_data,
    ggplot2::aes(x = Station, y = neg_log10_p, color = Method, shape = Method, linetype = Method)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = 276:290) +
    ggplot2::scale_color_manual(values = c("Fourier (pval_fourier)" = "#1f77b4", "Exact (pval_exact)" = "#ff7f0e")) +
    ggplot2::labs(
      title = "LLR p-values by polling station at Escuela Centenario",
      x = "Polling station",
      y = expression(-log[10](p)),
      color = "Method",
      shape = "Method",
      linetype = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
} else {
  matplot(
    276:290,
    -log10(cbind(p_fourier_all, p_exact_all)),
    type = "b",
    pch = c(19, 17),
    col = c("blue", "orange"),
    lty = c(1, 2),
    xlab = "Polling station",
    ylab = expression(-log[10](p)),
    xaxt = "n",
    main = "LLR p-values by polling station at Escuela Centenario"
  )
  axis(1, at = 276:290, labels = 276:290)
  legend("topright", legend = c("Fourier (pval_fourier)", "Exact (pval_exact)"),
         col = c("blue", "orange"), pch = c(19, 17), lty = c(1, 2), bty = "n")
}

## ----llr-ranking, include = FALSE---------------------------------------------
p_sorted <- sort(p_exact_all)

## ----compare-tests------------------------------------------------------------
test_specs <- list(
  "PMF"                   = list(stat = "pmf"),
  "Pearson (Chi2)"        = list(stat = "chi2"),
  "Log-Likelihood Ratio"  = list(stat = "llr"),
  "Power Divergence (2/3)"= list(stat = "pd", lambda = 2/3)
)

# Run a p-value function on station m with the statistic described by spec
run_test <- function(fun, m, spec, ...) {
  args <- list(x = votes[, m], p = loo_p[, m], stat = spec$stat, ...)
  if (!is.null(spec$lambda)) args$lambda <- spec$lambda
  do.call(fun, args)
}

fits_by_test <- lapply(
  test_specs,
  function(spec) lapply(colnames(votes), function(m) run_test(pval_fourier, m, spec, rel_eps = 1e-6))
)

p_by_test <- sapply(fits_by_test, function(fits) sapply(fits, `[[`, "pval"))
status_by_test <- sapply(fits_by_test, function(fits) sapply(fits, `[[`, "status"))
rownames(p_by_test) <- rownames(status_by_test) <- colnames(votes)

# A dagger marks results where the series did not converge
p_display <- matrix(fmt_p(p_by_test), nrow = nrow(p_by_test), dimnames = dimnames(p_by_test))
p_display[status_by_test != 0] <- paste0(p_display[status_by_test != 0], " †")

comparison <- data.frame(
  Station = colnames(votes),
  p_display,
  check.names = FALSE,
  row.names = NULL
)

knitr::kable(
  comparison,
  format = "html",
  escape = FALSE,
  align = "l",
  caption = "Leave-one-out p-values across different test statistics (Fourier inversion). † marks results where the series did not converge (status != 0)."
)

## ----compare-tests-278--------------------------------------------------------
exact_278 <- sapply(
  test_specs,
  function(spec) run_test(pval_exact, "278", spec)$pval
)

fourier_278 <- lapply(fits_by_test, function(fits) fits[[which(colnames(votes) == "278")]])

comparison_278 <- data.frame(
  Statistic = names(test_specs),
  p_Fourier = fmt_p(sapply(fourier_278, `[[`, "pval")),
  Terms = sapply(fourier_278, `[[`, "n_terms"),
  Status = sapply(fourier_278, `[[`, "status"),
  p_Exact = fmt_p(exact_278),
  row.names = NULL
)

knitr::kable(
  comparison_278,
  format = "html",
  escape = FALSE,
  align = "l",
  caption = "Station 278: Fourier inversion vs. exact p-values for each statistic."
)

