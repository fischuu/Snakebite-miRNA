# Shared helpers for the per-module PDF reports (workflow/scripts/R/*.Rmd).
# Sourced from every module report; keeps benchmark-aggregation and common
# file parsers in one place instead of duplicated across seven Rmd files.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
})

# ---- Runtime / resource statistics (from Snakemake `benchmark:` files) ----

#' Every benchmark TSV for one module lives directly under
#' "<project_folder>/benchmark/<module>/", named "<rule_slug>.<wildcards...>.tsv"
#' (or just "<rule_slug>.tsv" for wildcard-free rules). Combine them into one
#' data frame with a "rule" column derived from the filename.
discover_benchmarks <- function(project_folder, module) {
  root <- file.path(project_folder, "benchmark", module)
  if (!dir.exists(root)) {
    return(empty_benchmark_table())
  }

  paths <- list.files(root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(paths) == 0) {
    return(empty_benchmark_table())
  }

  rows <- map(paths, function(p) {
    df <- tryCatch(read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df$rule <- benchmark_label(p)
    df
  })
  rows <- compact(rows)
  if (length(rows) == 0) {
    return(empty_benchmark_table())
  }
  bind_rows(rows)
}

empty_benchmark_table <- function() {
  tibble(
    s = numeric(0), `h:m:s` = character(0), max_rss = numeric(0),
    max_vms = numeric(0), max_uss = numeric(0), max_pss = numeric(0),
    io_in = numeric(0), io_out = numeric(0), mean_load = numeric(0),
    cpu_time = numeric(0), rule = character(0)
  )
}

#' Derive the rule label from a benchmark file's name: everything before the
#' first "." (wildcards -- sample_id, library_id -- follow as further
#' dot-separated segments and never appear in the rule slug itself).
#'   "bowtie_mature.Sample1_test_S12.tsv" -> "bowtie_mature"
#'   "join_loci.tsv"                     -> "join_loci"
benchmark_label <- function(path) {
  base <- sub("\\.tsv$", "", basename(path))
  sub("\\..*$", "", base)
}

#' One row per rule: number of jobs, total/mean/max runtime, peak memory.
summarise_runtime <- function(bench_df) {
  if (nrow(bench_df) == 0) {
    return(tibble(
      rule = character(0), n_jobs = integer(0), total_runtime = character(0),
      mean_runtime = character(0), max_runtime = character(0), peak_mem_gb = numeric(0)
    ))
  }
  bench_df %>%
    group_by(rule) %>%
    summarise(
      n_jobs = n(),
      total_s = sum(s, na.rm = TRUE),
      mean_s = mean(s, na.rm = TRUE),
      max_s = max(s, na.rm = TRUE),
      peak_mem_gb = round(max(max_rss, na.rm = TRUE) / 1024, 2),
      .groups = "drop"
    ) %>%
    mutate(
      total_runtime = format_hms(total_s),
      mean_runtime = format_hms(mean_s),
      max_runtime = format_hms(max_s)
    ) %>%
    select(rule, n_jobs, total_runtime, mean_runtime, max_runtime, peak_mem_gb) %>%
    arrange(desc(total_runtime))
}

format_hms <- function(seconds) {
  seconds[is.na(seconds)] <- 0
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- round(seconds %% 60)
  sprintf("%02d:%02d:%02d", h, m, s)
}

#' kableExtra escapes table *cell* content for LaTeX automatically, but not
#' the `caption=` string itself -- a literal "_" there (e.g. a module name
#' like "novel_mirna") crashes LaTeX with "Missing $ inserted". Every caption
#' in this file is routed through here rather than escaped by hand.
escape_latex <- function(text) {
  gsub("([_%&#{}])", "\\\\\\1", text)
}

#' Nicely formatted runtime/resource table for the report.
runtime_table <- function(bench_df, caption = "Runtime and peak memory per rule") {
  tbl <- summarise_runtime(bench_df)
  if (nrow(tbl) == 0) {
    return(cat("_No benchmark data recorded for this module yet._\n"))
  }
  colnames(tbl) <- c("Rule", "Jobs", "Total runtime", "Mean runtime", "Max runtime", "Peak RAM (GB)")
  kable(tbl, caption = escape_latex(caption), booktabs = TRUE) %>%
    kable_styling(latex_options = c("striped", "hold_position"), font_size = 9)
}

#' Horizontal bar chart of total runtime per rule, longest first.
runtime_plot <- function(bench_df, title = "Total runtime per rule") {
  tbl <- bench_df %>%
    group_by(rule) %>%
    summarise(total_s = sum(s, na.rm = TRUE), .groups = "drop")
  if (nrow(tbl) == 0) {
    return(invisible(NULL))
  }
  tbl <- tbl %>% mutate(rule = reorder(rule, total_s))
  ggplot(tbl, aes(x = rule, y = total_s / 3600)) +
    geom_col(fill = "#00B5E2") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Total runtime (hours)") +
    theme_minimal(base_size = 11)
}

# ---- Allocated vs. actual resource usage (for calibrating escalation.yaml) ----

#' Load the pipeline's resource tiers (config.yaml's resource_sets) and the
#' per-rule escalation order (escalation.yaml). Both live at a fixed,
#' project-relative path regardless of which config files were passed to
#' this particular report (same convention the Snakefile itself uses).
load_resource_config <- function(project_folder) {
  list(
    resource_sets = yaml::read_yaml(file.path(project_folder, "config", "config.yaml"))$resource_sets,
    escalation = yaml::read_yaml(file.path(project_folder, "config", "escalation.yaml"))
  )
}

#' Best-effort match from a `discover_benchmarks()` rule label (e.g.
#' "bowtie_mature", "samtools_star_mature_flagstat") back to its real
#' escalation.yaml key (e.g. "align__bowtie__mature",
#' "align__samtools__star_mature_flagstat"), so its configured resource tier
#' can be looked up. Every word in the label must appear in the candidate
#' key. Ambiguous or unmatched rules return NA and the caller falls back to
#' the pipeline's default 3-tier escalation, rather than guessing wrong.
match_escalation_key <- function(rule_label, module, escalation) {
  words <- tolower(strsplit(gsub("[:_]", " ", rule_label), "\\s+")[[1]])
  words <- words[nchar(words) > 0]
  candidates <- names(escalation)[startsWith(names(escalation), paste0(module, "__"))]
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  hits <- vapply(candidates, function(k) all(vapply(words, grepl, logical(1), x = tolower(k), fixed = TRUE)), logical(1))
  matched <- candidates[hits]
  if (length(matched) == 1) matched else NA_character_
}

#' For each rule with benchmark data, compare its actual peak memory/runtime
#' against what its first-attempt resource tier allocates. This is the
#' escalation.yaml/config.yaml "small/medium/large" tier a rule starts on
#' before any retry-escalation -- a rule sitting near 100% utilization here
#' is a candidate for a bigger starting tier; one sitting well under is a
#' candidate for a smaller (cheaper, faster-scheduled) one.
resource_allocation_table <- function(bench_df, module, project_folder,
                                      caption = "Allocated vs. actual resource usage") {
  if (nrow(bench_df) == 0) {
    return(cat("_No benchmark data recorded for this module yet._\n"))
  }
  cfg <- load_resource_config(project_folder)
  default_tiers <- c("small", "medium", "large")

  actual <- bench_df %>%
    group_by(rule) %>%
    summarise(actual_mem_gb = max(max_rss, na.rm = TRUE) / 1024,
              actual_runtime_s = max(s, na.rm = TRUE), .groups = "drop")

  rows <- pmap_dfr(actual, function(rule, actual_mem_gb, actual_runtime_s) {
    key <- match_escalation_key(rule, module, cfg$escalation)
    tiers <- if (!is.na(key)) cfg$escalation[[key]] else default_tiers
    tier1 <- trimws(tiers[[1]])
    alloc <- cfg$resource_sets[[tier1]]
    if (is.null(alloc)) {
      return(NULL)
    }
    alloc_mem_gb <- alloc$mem_mb / 1024
    alloc_runtime_s <- alloc$runtime * 60
    mem_pct <- round(100 * actual_mem_gb / alloc_mem_gb, 1)
    runtime_pct <- round(100 * actual_runtime_s / alloc_runtime_s, 1)
    tibble(
      rule = rule,
      tier = tier1,
      allocated_mem_gb = round(alloc_mem_gb, 1),
      actual_mem_gb = round(actual_mem_gb, 1),
      mem_used_pct = mem_pct,
      allocated_runtime = format_hms(alloc_runtime_s),
      actual_runtime = format_hms(actual_runtime_s),
      runtime_used_pct = runtime_pct,
      flag = case_when(
        pmax(mem_pct, runtime_pct) >= 90 ~ "near limit -- consider a bigger starting tier",
        pmax(mem_pct, runtime_pct) <= 15 ~ "over-provisioned -- consider a smaller tier",
        TRUE ~ ""
      )
    )
  })

  if (is.null(rows) || nrow(rows) == 0) {
    return(cat("_Could not match any rule in this module to a configured resource tier._\n"))
  }

  colnames(rows) <- c("Rule", "Starting tier", "Allocated RAM (GB)", "Actual RAM (GB)", "RAM used (%)",
                       "Allocated runtime", "Actual runtime", "Runtime used (%)", "Note")
  kable(rows, caption = escape_latex(caption), booktabs = TRUE) %>%
    kable_styling(latex_options = c("striped", "hold_position", "scale_down"), font_size = 7)
}

# ---- Generic table helpers ----

#' Read a TSV if it exists, else return an empty (zero-row) tibble so
#' downstream code can render "no data yet" instead of erroring.
safe_read_tsv <- function(path, ...) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble())
  }
  tryCatch(read_tsv(path, show_col_types = FALSE, ...), error = function(e) tibble())
}

#' Preview the first n rows of a (possibly large) table.
preview_table <- function(df, n = 10, caption = NULL) {
  if (nrow(df) == 0) {
    return(cat("_No data available yet._\n"))
  }
  if (!is.null(caption)) caption <- escape_latex(caption)
  kable(head(df, n), caption = caption, booktabs = TRUE) %>%
    kable_styling(latex_options = c("striped", "hold_position", "scale_down"), font_size = 8)
}

# ---- Tool-specific parsers ----

#' FastQC's fastqc_data.txt lives inside each *_fastqc.zip; extract the
#' ">>Basic Statistics" block as a one-row data frame.
parse_fastqc_zip <- function(zip_path) {
  if (!file.exists(zip_path)) {
    return(tibble())
  }
  inner <- grep("fastqc_data.txt$", unzip(zip_path, list = TRUE)$Name, value = TRUE)
  if (length(inner) == 0) {
    return(tibble())
  }
  con <- unz(zip_path, inner[[1]])
  lines <- readLines(con, warn = FALSE)
  close(con)
  start <- grep("^>>Basic Statistics", lines)
  end <- grep("^>>END_MODULE", lines)
  end <- end[end > start][1]
  block <- lines[(start + 2):(end - 1)]
  kv <- strsplit(block, "\t")
  df <- as_tibble(setNames(
    as.list(map_chr(kv, 2)),
    map_chr(kv, 1)
  ))
  df$sample <- basename(zip_path)
  df
}

#' samtools flagstat's plain-text output: total, mapped, and mapped % reads.
parse_flagstat <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  lines <- readLines(path, warn = FALSE)
  total_line <- grep("in total", lines, value = TRUE)
  mapped_line <- grep("mapped \\(", lines, value = TRUE)
  total <- if (length(total_line) > 0) as.numeric(sub("^([0-9]+).*", "\\1", total_line[1])) else NA
  mapped <- if (length(mapped_line) > 0) as.numeric(sub("^([0-9]+).*", "\\1", mapped_line[1])) else NA
  pct <- if (length(mapped_line) > 0) {
    m <- regmatches(mapped_line[1], regexpr("[0-9.]+(?=%)", mapped_line[1], perl = TRUE))
    if (length(m) > 0) as.numeric(m) else NA
  } else {
    NA
  }
  tibble(total_reads = total, mapped_reads = mapped, mapped_pct = pct)
}

#' Collect one flagstat per sample into a single sample x metric table.
collect_flagstats <- function(project_folder, path_template, samples) {
  rows <- map(samples, function(s) {
    path <- file.path(project_folder, sprintf(path_template, s))
    df <- parse_flagstat(path)
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$sample_id <- s
    df
  })
  rows <- compact(rows)
  if (length(rows) == 0) {
    return(tibble())
  }
  bind_rows(rows) %>% select(sample_id, everything())
}

#' seqkit bam -C's stderr-derived count report: one line per reference.
parse_seqkit_bam_count <- function(path) {
  df <- safe_read_tsv(path, col_names = FALSE)
  if (nrow(df) == 0) {
    return(df)
  }
  df
}

#' featureCounts summary table (rows: Assigned/Unassigned_* ; one count column).
parse_featurecounts_summary <- function(path) {
  summary_path <- paste0(path, ".summary")
  df <- safe_read_tsv(summary_path)
  if (nrow(df) == 0) {
    return(df)
  }
  colnames(df)[1:2] <- c("status", "count")
  df
}

#' mirdb.stats holds 4 raw counts (mature, mature_species, hairpin, hairpin_species),
#' one integer per line, in that fixed order (see reference/mirdb.smk).
parse_mirdb_stats <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  counts <- suppressWarnings(as.numeric(readLines(path, warn = FALSE)))
  tibble(
    reference = c("mature.fa", "mature_species.fa", "hairpin.fa", "hairpin_species.fa")[seq_along(counts)],
    n_sequences = counts
  )
}
