# 1 Set up
# 1.1 Installation
packages <- c("CBPS", "cobalt", "data.table", "dplyr", "DT", "DoubleML", "ebal", "estimatr", "ggplot2", "gridExtra", "grf", "hbal", "highr", "highs", 
              "kableExtra", "MatchIt", "Matching", "mlr3", "mlr3learners", "OutcomeWeights", "optmatch", "optweight", "quickmatch", 
              "readr", "rgenoud", "tidyr", "tidyverse", "ppcor", "patchwork", "WeightIt"
)

#### install_all()
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}

install_all(packages)

# load packages
library(CBPS)
library(cobalt)
library(data.table)
library(dplyr)
library(DT)
library(DoubleML)
library(ebal)
library(ggplot2)
library(gridExtra)
library(grf)
library(hbal)
library(highr)
library(highs)
library(kableExtra)
library(MatchIt)
library(Matching)
library(mlr3)
library(mlr3learners)
library(optmatch)
library(optweight)
library(quickmatch)
library(readr)
library(rgenoud)
library(tidyr)
library(tidyverse)
library(ppcor)
library(patchwork)
library(WeightIt)
library(OutcomeWeights)

# 1.2 Data inspection
#### inspect_datasets()
inspect_data <- function(data, treat = "treat") {
  if (is.data.frame(data)) {
    data <- list(dataset = data)
  }
  data.frame(
    dataset = names(data),
    num_obs = sapply(data, nrow),
    num_treated = sapply(data, function(df) sum(df[[treat]] == 1, na.rm = TRUE)),
    num_controls = sapply(data, function(df) sum(df[[treat]] == 0, na.rm = TRUE)),
    num_vars = sapply(data, ncol),
    name_vars = sapply(data, function(df) paste(names(df), collapse = ", ")),
    row.names = NULL
  )
}

#### save_overlap_panels()
save_overlap_panels <- function(
    data_list,
    treat,
    covar,
    plot_titles = NULL,
    prefix = "overlap_panels",
    folder = "graphs/lalonde",
    plots_per_page = 1
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n <- length(data_list)
  pages <- ceiling(n / plots_per_page)
  for (p in seq_len(pages)) {
    start <- (p - 1) * plots_per_page + 1
    end   <- min(p * plots_per_page, n)
    file_name <- file.path(folder, paste0(prefix, "_", p, ".pdf"))
    pdf(file_name, width = 8, height = 2.75 * (end - start + 1))
    par(mfrow = c(end - start + 1, 1), mar = c(4, 5, 3, 2))
    for (i in start:end) {
      assess_overlap(data = data_list[[i]], treat = treat, cov = covar)
      if (!is.null(plot_titles) && length(plot_titles) >= i) {
        title(main = plot_titles[i])
      } else {
        title(main = paste("Panel", i))
      }
    }
    dev.off()
  }
}

# 2. Improving covariate balance and overlap
## 2.3 Trimming
### 2.3.1 Propensity score threshold trimming (similar to tutorial of Imbens & Xu (2024))
# ps_trim()*** 
ps_trim <- function(data, ps = "ps_assoverlap", threshold = 0.9) { 
  sub <- data[which(data[, ps] < threshold), ]
  return(sub)
}

### 2.3.2 Common range trimming
#### common_range_trim()
common_range_trim <- function(data, ps = "ps_assoverlap", treat = "treat") {
  lower_cut <- min(data[[ps]][data[[treat]] == 1], na.rm = TRUE)
  upper_cut <- max(data[[ps]][data[[treat]] == 0], na.rm = TRUE)
  sub <- data[data[[ps]] >= lower_cut & data[[ps]] <= upper_cut, ]
  return(sub)
}

### 2.3.3 Crump trimming
#### crump_trim()
crump_trim <- function(data, ps = "ps_assoverlap", lower = 0.1, upper = 0.9) {
  sub <- data[data[[ps]] >= lower & data[[ps]] <= upper, ]
  return(sub)
}

### 2.3.4 Stuermer trimming
#### stuermer_trim()
stuermer_trim <- function(data, treat = "treat", ps = "ps_assoverlap", 
                          lower_percentile = 0.01, upper_percentile = 0.99) {
  treated_ps   <- data[[ps]][data[[treat]] == 1]
  untreated_ps <- data[[ps]][data[[treat]] == 0]
  lower_cutoff <- quantile(treated_ps, probs = lower_percentile, na.rm = TRUE)
  upper_cutoff <- quantile(untreated_ps, probs = upper_percentile, na.rm = TRUE)
  sub <- data[data[[ps]] >= lower_cutoff & data[[ps]] <= upper_cutoff, ]
  return(sub)
}

### 2.3.5 Walker trimming
#### walker_trim()
walker_trim <- function(data, treat = "treat", ps = "ps_assoverlap", 
                        lower_cutoff = 0.1, upper_cutoff = 0.9) {
  treat_prevalence  <- mean(data[[treat]], na.rm = TRUE)
  logit_ps          <- log(data[[ps]] / (1 - data[[ps]]))
  logit_prevalence  <- log(treat_prevalence / (1 - treat_prevalence))
  preference_score  <- 1 / (1 + exp(-(logit_ps - logit_prevalence)))
  sub <- data[preference_score >= lower_cutoff & preference_score <= upper_cutoff, ]
  return(sub)
}

#### ps_estimate()***
ps_estimate <- function(data, treat, cov, num.trees = NULL, seed = 42) {
  if(is.null(num.trees)){
    p.forest1 <- probability_forest(
      X = data[, cov],
      Y = as.factor(data[,treat]), seed = seed)
  } else {
    p.forest1 <- probability_forest(
      X = data[, cov],
      Y = as.factor(data[,treat]), seed = seed, num.trees = num.trees)
  }
  data$ps_assoverlap <- p.forest1$predictions[,2]
  data$ps_assoverlap[which(abs(data$ps_assoverlap) <= 1e-7)] <- 1e-7
  return(data)
}

## 2.5 Integrated methods
### 2.5.1 Trimming and matching
#### attach_matchit()
attach_matchit <- function(model, data_list, ..., verbose = FALSE) {
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("attach_matchit(): 'data_list' must be a non-empty list", call. = FALSE)
  }
  matchit_results <- vector("list", length(data_list))
  names(matchit_results) <- names(data_list)
  dims_report <- character(length(data_list))
  for (i in seq_along(data_list)) {
    data_i <- data_list[[i]]
    label_i <- if (length(names(matchit_results))) names(matchit_results)[i] else i
    res <- tryCatch(
      {
        obj <- matchit(model, data = data_i, ...)
        attr(obj, "match_source") <- data_i
        obj
      },
      error = function(e) {
        warning(sprintf("attach_matchit(): '%s' failed: %s", label_i, e$message), call. = FALSE)
        NULL
      }
    )
    matchit_results[[i]] <- res
    source_rows <- tryCatch(nrow(data_i), error = function(e) NA_integer_)
    matched_rows <- if (inherits(res, "matchit")) {
      tryCatch(nrow(match.data(res, data = data_i)), error = function(e) NA_integer_)
    } else {
      NA_integer_
    }
    dims_report[[i]] <- sprintf("%s: source=%s matched=%s", label_i, source_rows, matched_rows)
  }
  if (isTRUE(verbose)) {
    message("attach_matchit(): ", paste(dims_report, collapse = "; "))
  }
  attr(matchit_results, "match_dims") <- dims_report
  # report summary
  null_count <- sum(sapply(matchit_results, is.null))
  success_count <- length(matchit_results) - null_count
  if (null_count > 0) {
    failed_names <- names(matchit_results)[sapply(matchit_results, is.null)]
    cat(sprintf("MATCHING SUMMARY: %d succeeded, %d failed. Failed: %s\n", 
                success_count, null_count, paste(failed_names, collapse = ", ")))
  }
  return(matchit_results)
}

# 3. Assessing methods
## 3.1 Trimming 
### 3.1.1 SMD
#### compute_abs_smd_trim()
compute_abs_smd_trim <- function(data_list, treat = "treat", covar) {
  smd_list <- lapply(names(data_list), function(name) {
    data <- data_list[[name]]
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      un = TRUE,
      stats = "mean.diffs",
      s.d.denom = "treated"
    )
    smds <- bal_obj$Balance$Diff.Un   
    smd_vals <- abs(smds)
    mean_smd <- mean(smd_vals, na.rm = TRUE)
    max_smd  <- max(smd_vals, na.rm = TRUE)
    return(data.frame(
      Method   = name,
      Mean_Abs_SMD = mean_smd,
      Max_Abs_SMD  = max_smd
    ))
  })
  smd_summary <- do.call(rbind, smd_list)
  rownames(smd_summary) <- NULL
  return(smd_summary)
}

### 3.1.2 OVL
#### compute_ovl_trim()
compute_ovl_trim <- function(data_list, ps = "ps_assoverlap", treat = "treat", n_points = 512) {
  ovl_res <- lapply(names(data_list), function(meth) {
    dat <- data_list[[meth]]
    ps <- as.numeric(dat[[ps]])
    treat <- dat[[treat]]
    valid_idx <- !is.na(ps) & !is.na(treat) & (ps >= 0 & ps <= 1)
    ps <- ps[valid_idx]
    treat <- treat[valid_idx]
    treated_scores <- ps[treat == 1]
    control_scores <- ps[treat == 0]
    dens_treated <- density(treated_scores, from = 0, to = 1, n = n_points)
    dens_control <- density(control_scores, from = 0, to = 1, n = n_points)
    x_grid <- dens_treated$x
    min_density <- pmin(dens_treated$y, dens_control$y)
    bin_width <- x_grid[2] - x_grid[1]
    ovl <- sum(min_density) * bin_width
    data.frame(Method = meth, OVL = ovl)
  })
  res <- do.call(rbind, ovl_res)
  rownames(res) <- NULL
  return(res)
}

## 3.2 Trimming and matching
### 3.2.1 SMD
#### compute_abs_smd_matchit()
compute_abs_smd_matchit <- function(match_list, data_list) {
  if (!is.list(match_list) || length(match_list) == 0) {
    stop("compute_abs_smd_matchit(): 'match_list' must be a non-empty list", call. = FALSE)
  }
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("compute_abs_smd_matchit(): 'data_list' must be a non-empty list", call. = FALSE)
  }
  data_len <- length(data_list)
  data_labels <- if (length(names(data_list))) names(data_list) else seq_len(data_len)
  smd_results <- list()
  for (method_name in names(match_list)) {
    method_sublist <- match_list[[method_name]]
    if (!is.list(method_sublist)) {
      warning(sprintf("compute_abs_smd_matchit(): '%s' is not a list of matchit objects; skipping", method_name), call. = FALSE)
      next
    }
    method_len <- length(method_sublist)
    if (method_len != data_len) {
      warning(sprintf(
        "compute_abs_smd_matchit(): '%s' match list length (%d) differs from data_list length (%d); only the first %d elements are considered.",
        method_name,
        method_len,
        data_len,
        min(method_len, data_len)
      ), call. = FALSE)
    }
    max_idx <- min(method_len, data_len)
    for (i in seq_len(max_idx)) {
      match_obj <- method_sublist[[i]]
      data_i <- data_list[[i]]
      label_i <- if (i <= length(data_labels)) data_labels[[i]] else i
      if (is.null(match_obj)) {
        warning(sprintf(
          "compute_abs_smd_matchit(): skipping '%s' match object %d (%s) because it is NULL; check attach_matchit() failures.",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      if (is.null(data_i)) {
        warning(sprintf(
          "compute_abs_smd_matchit(): skipping '%s' source data %d (%s) because it is NULL; verify trimming outputs.",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      bal <- bal.tab(
        match_obj,
        data = data_i,
        stats = "mean.diffs",
        un = TRUE,
        s.d.denom = "treated"
      )
      smds <- bal$Balance$Diff.Adj
      smd_vals <- abs(smds)
      mean_smd <- mean(smd_vals, na.rm = TRUE)
      max_smd <- max(smd_vals, na.rm = TRUE)
      smd_results[[paste(method_name, label_i, sep = "_")]] <- data.frame(
        Method = method_name,
        MatchIndex = label_i,
        Mean_Abs_SMD = mean_smd,
        Max_Abs_SMD = max_smd
      )
    }
  }
  if (length(smd_results) == 0) {
    return(data.frame(
      Method = character(0),
      MatchIndex = integer(0),
      Mean_Abs_SMD = numeric(0),
      Max_Abs_SMD = numeric(0)
    ))
  }
  smd_summary <- do.call(rbind, smd_results)
  rownames(smd_summary) <- NULL
  return(smd_summary)
}

### 3.2.2 OVL
#### compute_ovl_matchit()
compute_ovl_matchit <- function(match_list, data_list, ps = "ps_assoverlap", treat = "treat",
                                covar = NULL, num.trees = NULL, seed = 42, n_points = 512) {
  if (!is.list(match_list) || length(match_list) == 0) {
    stop("compute_ovl_matchit(): 'match_list' must be a non-empty list", call. = FALSE)
  }
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("compute_ovl_matchit(): 'data_list' must be a non-empty list", call. = FALSE)
  }
  data_len <- length(data_list)
  data_labels <- if (length(names(data_list))) names(data_list) else seq_len(data_len)
  ovl_results <- list()
  for (method_name in names(match_list)) {
    method_sublist <- match_list[[method_name]]
    if (!is.list(method_sublist)) {
      warning(sprintf("compute_ovl_matchit(): '%s' is not a list of matchit objects; skipping", method_name), call. = FALSE)
      next
    }
    method_len <- length(method_sublist)
    if (method_len != data_len) {
      warning(sprintf(
        "compute_ovl_matchit(): '%s' match list length (%d) differs from data_list length (%d); only the first %d elements are considered.",
        method_name,
        method_len,
        data_len,
        min(method_len, data_len)
      ), call. = FALSE)
    }
    max_idx <- min(method_len, data_len)
    for (i in seq_len(max_idx)) {
      match_obj <- method_sublist[[i]]
      data_i <- data_list[[i]]
      label_i <- if (i <= length(data_labels)) data_labels[[i]] else i
      if (is.null(match_obj)) {
        warning(sprintf(
          "compute_ovl_matchit(): skipping '%s' match object %d (%s) because it is NULL; check attach_matchit() failures.",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      if (!inherits(match_obj, "matchit")) {
        warning(sprintf(
          "compute_ovl_matchit(): skipping '%s' match object %d (%s) because it is class '%s' not 'matchit'",
          method_name,
          i,
          label_i,
          paste(class(match_obj), collapse = "/")
        ), call. = FALSE)
        next
      }
      if (is.null(data_i)) {
        data_i <- attr(match_obj, "match_source")
        if (is.null(data_i)) {
          warning(sprintf(
            "compute_ovl_matchit(): source data missing for method '%s' index %d (%s); cannot recover matched data.",
            method_name,
            i,
            label_i
          ), call. = FALSE)
          next
        }
      }
      if (!is.data.frame(data_i)) {
        warning(sprintf(
          "compute_ovl_matchit(): data_list[[%d]] for method '%s' is type '%s' not data.frame",
          i,
          method_name,
          paste(class(data_i), collapse = "/")
        ), call. = FALSE)
        next
      }
      if (!nrow(data_i)) {
        warning(sprintf(
          "compute_ovl_matchit(): data_list[[%d]] for method '%s' has zero rows",
          i,
          method_name
        ), call. = FALSE)
        next
      }
      matched_data <- match.data(match_obj, data = data_i)
      if (!is.data.frame(matched_data) || !nrow(matched_data)) {
        warning(sprintf(
          "compute_ovl_matchit(): matched data empty for method '%s' index %d (%s)",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      if (ps %in% names(matched_data)) {
        ps_vals <- as.numeric(matched_data[[ps]])
      } else {
        if (is.null(covar)) {
          stop(sprintf(
            "compute_ovl_matchit(): propensity scores missing for '%s' index %d (%s) and 'covar' not provided",
            method_name,
            i,
            label_i
          ), call. = FALSE)
        }
        X <- matched_data[, covar, drop = FALSE]
        Y <- as.factor(matched_data[[treat]])
        if (is.null(num.trees)) {
          p.forest1 <- probability_forest(X = X, Y = Y, seed = seed)
        } else {
          p.forest1 <- probability_forest(X = X, Y = Y, seed = seed, num.trees = num.trees)
        }
        ps_vals <- p.forest1$predictions[, 2]
        ps_vals[which(abs(ps_vals) <= 1e-7)] <- 1e-7
      }
      treat_vals <- matched_data[[treat]]
      valid_idx <- !is.na(ps_vals) & !is.na(treat_vals) & (ps_vals >= 0 & ps_vals <= 1)
      ps_vals <- ps_vals[valid_idx]
      treat_vals <- treat_vals[valid_idx]
      if (!length(ps_vals) || length(unique(treat_vals)) < 2) {
        warning(sprintf(
          "compute_ovl_matchit(): insufficient treated/control overlap after cleaning for method '%s' index %d (%s)",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      treated_scores <- ps_vals[treat_vals == 1]
      control_scores <- ps_vals[treat_vals == 0]
      if (length(treated_scores) < 2 || length(control_scores) < 2) {
        warning(sprintf(
          "compute_ovl_matchit(): Not enough treated/control units for density estimation in method '%s' index %d (%s)",
          method_name,
          i,
          label_i
        ), call. = FALSE)
        next
      }
      dens_treated <- density(treated_scores, from = 0, to = 1, n = n_points)
      dens_control <- density(control_scores, from = 0, to = 1, n = n_points)
      x_grid <- dens_treated$x
      min_density <- pmin(dens_treated$y, dens_control$y)
      bin_width <- x_grid[2] - x_grid[1]
      ovl <- sum(min_density) * bin_width
      ovl_results[[paste(method_name, label_i, sep = "_")]] <- data.frame(
        Method = paste(method_name, label_i, sep = "_"),
        OVL = ovl
      )
    }
  }
  if (length(ovl_results) == 0) {
    return(data.frame(Method = character(0), OVL = numeric(0)))
  }
  final_df <- do.call(rbind, ovl_results)
  rownames(final_df) <- NULL
  return(final_df)
}

## 4. Identifying best methods
### 4.1 Ranking
#### combine_results()
combine_results <- function(dataset_name) {
  dataset_lower <- tolower(dataset_name)
  # retrieve individual method results
  smd_trimming <- get(paste0("smd_trim.", dataset_lower))
  ovl_trimming <- get(paste0("ovl_trim.", dataset_lower))
  smd_trim_match_combined <- get(paste0("smd_trim_match_comb.", dataset_lower))
  ovl_trim_match_combined <- get(paste0("ovl_trim_match_comb.", dataset_lower))
  format_smd_results <- function(df) {
    if (is.null(df) || !nrow(df)) {
      return(data.frame(Method = character(0), Mean_Abs_SMD = numeric(0), Max_Abs_SMD = numeric(0)))
    }
    out <- df
    if ("MatchIndex" %in% names(out)) {
      out$Method <- paste(out$Method, out$MatchIndex, sep = "_")
      out$MatchIndex <- NULL
    }
    out[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD"), drop = FALSE]
  }
  format_ovl_results <- function(df) {
    if (is.null(df) || !nrow(df)) {
      return(data.frame(Method = character(0), OVL = numeric(0)))
    }
    df[, c("Method", "OVL"), drop = FALSE]
  }
  smd_all <- dplyr::bind_rows(
    format_smd_results(smd_trimming),
    format_smd_results(smd_trim_match_combined)
  )
  ovl_all <- dplyr::bind_rows(
    format_ovl_results(ovl_trimming),
    format_ovl_results(ovl_trim_match_combined)
  )
  # merge absolute SMD and OVL results by method
  final_df <- dplyr::full_join(smd_all, ovl_all, by = "Method")
  if (!nrow(final_df)) {
    return(final_df)
  }
  final_df$Method <- as.character(final_df$Method)
  # remove dataset suffixes for clean labels
  final_df$Method <- gsub("\\.psid", "", final_df$Method, ignore.case = TRUE)
  final_df$Method <- gsub("\\.cps", "", final_df$Method, ignore.case = TRUE)
  # reset row names
  rownames(final_df) <- NULL
  return(final_df)
}

#### save_csv()
save_csv <- function(data, filename) {
  folder <- "tables"
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  file_csv <- file.path(folder, paste0(filename, ".csv"))
  write.csv(data, file = file_csv, row.names = FALSE)
}

#### assess_methods()
assess_methods <- function(data) {
  data %>%
    dplyr::select(Method, OVL, Mean_Abs_SMD) %>%
    arrange(desc(OVL))
}

#### get_top_methods()
get_top_methods <- function(summary_df, top_n = 5, score_col = NULL) {
  if (!"Method" %in% names(summary_df)) {
    stop("Data frame must contain column 'Method'")
  }
  if (is.null(score_col)) {
    if ("Score" %in% names(summary_df)) {
      score_col <- "Score"
    } else if ("OVL" %in% names(summary_df)) {
      score_col <- "OVL"
    } else {
      stop("Data frame must contain a scoring column ('Score' or 'OVL'), or specify 'score_col'")
    }
  }
  if (!score_col %in% names(summary_df)) {
    stop(sprintf("Column '%s' not found in summary_df", score_col))
  }
  summary_df %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    arrange(dplyr::desc(.data[[score_col]])) %>%
    head(top_n) %>%
    dplyr::pull(Method)
}

#### rerank_methods()
rerank_methods <- function(top_methods, ranked_table) {
  filtered <- ranked_table[ranked_table$Method %in% top_methods, ]
  filtered <- filtered[order(-filtered$Mean_Abs_SMD), ]
  filtered$Method
}

### 4.2 Dataset construction
#### wrap_match_entries()
wrap_match_entries <- function(match_list, source_list, prefix) {
  if (!is.list(match_list)) {
    return(list())
  }
  max_idx <- min(length(match_list), length(source_list))
  source_labels <- if (length(names(source_list))) names(source_list) else seq_len(max_idx)
  entries <- vector("list", max_idx)
  names(entries) <- paste0(prefix, "_", source_labels[seq_len(max_idx)])
  for (i in seq_len(max_idx)) {
    match_obj <- match_list[[i]]
    if (is.null(match_obj)) {
      next
    }
    entries[[i]] <- list(matchit = match_obj, data = source_list[[i]])
  }
  Filter(Negate(is.null), entries)
}

#### create_top5_datasets()
create_top5_datasets <- function(dataset_list, top_method_names) {
  if (!is.list(dataset_list) || length(dataset_list) == 0) {
    stop("create_top5_datasets(): 'dataset_list' must be a non-empty list", call. = FALSE)
  }
  if (length(top_method_names) == 0) {
    stop("create_top5_datasets(): 'top_method_names' must contain at least one method name", call. = FALSE)
  }
  missing_methods <- setdiff(top_method_names, names(dataset_list))
  if (length(missing_methods)) {
    stop(sprintf(
      "create_top5_datasets(): the following methods are missing from 'dataset_list': %s",
      paste(missing_methods, collapse = ", ")
    ), call. = FALSE)
  }
  out <- lapply(top_method_names, function(method_name) {
    ds <- dataset_list[[method_name]]
    if (is.null(ds)) {
      stop(sprintf("create_top5_datasets(): dataset for method '%s' is NULL", method_name), call. = FALSE)
    }
    if (inherits(ds, "matchit")) {
      source_data <- attr(ds, "match_source")
      if (is.null(source_data)) {
        stop(sprintf(
          "create_top5_datasets(): matchit object '%s' is missing source data. Regenerate it with attach_matchit() so 'match_source' attribute is available or provide a data.frame entry instead.",
          method_name
        ), call. = FALSE)
      }
      return(as.data.frame(match.data(ds, data = source_data)))
    }
    if (is.list(ds) && inherits(ds$matchit, "matchit") && is.data.frame(ds$data)) {
      return(as.data.frame(match.data(ds$matchit, data = ds$data)))
    }
    if (is.data.frame(ds)) {
      return(ds)
    }
    stop(sprintf("create_top5_datasets(): unsupported entry type for method '%s'", method_name), call. = FALSE)
  })
  names(out) <- top_method_names
  out
}  

#### save_top5_datasets()
save_top5_datasets <- function(combined_methods_list, top_method_names, prefix) {
  dir.create("tutorial/data", showWarnings = FALSE, recursive = TRUE)
  for (i in seq_along(top_method_names)) {
    method_name <- top_method_names[i]
    if (!method_name %in% names(combined_methods_list)) {
      warning(paste0("Method '", method_name, "' not found in combined methods list"))
      next
    }
    dataset_to_save <- combined_methods_list[[method_name]]
    file_name <- sprintf("tutorial/data/top%d_%s_method_%s.RData", i, prefix, method_name)
    save(dataset_to_save, file = file_name)
  }
}

# 5. Estimating
## 5.1 ATT
#### quiet()***
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#### diff()***
diff <- function(data, Y, treat) {
  fml <- as.formula(paste(Y, "~", treat))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  return(out) # extract coef, se, ci.lower, ci.upper
}

#### reg()***
reg <- function(data, Y, treat, covar) {
  fml <- as.formula(paste(Y, "~", treat, "+", paste(covar, collapse = " + ")))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  # extract coef, se, ci.lower, ci.upper
  return(out)
}

# library(Matching)
#### matching()***
matching <- function(data, Y, treat, covar) {
  tryCatch({
    m.out <- Match(Y = data[, Y], Tr = data[, treat], X = data[, covar], Z = data[, covar],
                   estimand = "ATT", M = 5, replace = TRUE, ties = TRUE, BiasAdjust = TRUE)
    out <- c(m.out$est[1], m.out$se[1], m.out$est[1] - 1.96 * m.out$se[1],
             m.out$est[1] + 1.96 * m.out$se[1])
    return(out)
  }, error = function(e) {
    cat("Warning in matching method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

#### psm()***
psm <- function(data, Y, treat, covar) {
  ps <- probability_forest(X = data[, covar],
                           Y = as.factor(data[,treat]), seed = 42, num.trees = 4000)$predictions[,2]
  m.out <- Match(Y = data[, Y], Tr = data[, treat], X = matrix(ps, nrow(data), 1),
                 estimand = "ATT", M = 5, replace = TRUE, ties = FALSE, BiasAdjust = FALSE)
  if (is.null(m.out$se)==FALSE) {
    se <- m.out$se[1]
  } else {
    se <- m.out$se.standard[1]
  }
  out <- c(m.out$est[1], se, m.out$est[1] - 1.96 * se,
           m.out$est[1] + 1.96 * se)
  return(out)
}

#### om.reg()***
om.reg <- function(data, Y, treat, covar) {
  tr <- which(data[, treat] == 1)
  co <- which(data[, treat] == 0)
  fml <- as.formula(paste(Y, "~", paste(covar, collapse = " + ")))
  out.co <- lm(fml, data = data[co, ])
  Y.tr.hat <- predict(out.co, newdata = data[tr, covar, drop = FALSE])
  newdata <- cbind.data.frame(Y = c(data[tr, Y], Y.tr.hat), treat = rep(c(1, 0), each = length(tr)))
  out <- summary(lm_robust(Y ~ treat, data = newdata, se_type = "stata"))$coefficients["treat", c(1, 2, 5, 6)]
  return(out)
}

# library(grf)
#### om.grf()***
om.grf <- function(data, Y, treat, covar) {
  tr <- which(data[, treat] == 1)
  co <- which(data[, treat] == 0)
  out.co <- regression_forest(X = data[co, covar, drop = FALSE], Y = as.vector(data[co, Y]) )
  Y.tr.hat <- as.vector(unlist(predict(out.co, newdata = data[tr, covar, drop = FALSE])))
  newdata <- cbind.data.frame(Y = c(data[tr, Y], Y.tr.hat), treat = rep(c(1, 0), each = length(tr)))
  out <- summary(lm_robust(Y ~ treat, data = newdata, se_type = "stata"))$coefficients["treat", c(1, 2, 5, 6)]
  return(out)
}

#### ipw()***
ipw <- function(data, Y, treat, covar) {
  ps <- probability_forest(X = data[, covar, drop = FALSE], Y = as.factor(data[, treat]), seed = 42)$predictions[,2]
  fml <- as.formula(paste(Y, "~", treat))
  weights <- rep(1, nrow(data))
  co <- which(data[, treat] == 0)
  weights[co] <- ps[co]/(1-ps[co])
  out <- summary(lm_robust(fml, data = data, weights = weights, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  # extract coef, se, ci.lower, ci.upper
  return(out)
}

# library("CBPS")
#### cbps()***
cbps <- function(data, Y, treat, covar) {
  tryCatch({
    fml <- as.formula(paste(treat, "~", paste(covar, collapse = " + ")))
    ps <- quiet(CBPS(fml, data = data, standardize = TRUE)$fitted.values)
    fml <- as.formula(paste(Y, "~", treat))
    weights <- rep(1, nrow(data))
    co <- which(data[, treat] == 0)
    weights[co] <- ps[co]/(1-ps[co])
    out <- summary(lm_robust(fml, data = data, weights = weights, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
    return(out)
  }, error = function(e) {
    cat("Warning in cbps method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

# library(hbal)
#### ebal()***
ebal <- function(data, Y, treat, covar) {
  tryCatch({
    ebal.out <- hbal::hbal(Y = Y, Treat = treat, X = covar,  data = data, expand.degree = 1)
    out <- hbal::att(ebal.out, dr = FALSE)[1, c(1, 2, 5, 6)]
    return(out)
  }, error = function(e) {
    cat("Warning in ebal method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

#### hbal()***
# hbal <- function(data, Y, treat, covar) {
#   hbal.out <- hbal::hbal(Y = Y, Treat = treat, X = covar,  data = data, expand.degree = 2, # cv = TRUE)
#   out <- hbal::att(hbal.out, dr = FALSE)[1, c(1, 2, 5, 6)]
#   return(out)
# }

#### aipw_grf()***
aipw <- function(data, Y, treat, covar) {
  tryCatch({
    #library("grf")
    for (var in c(Y, treat, covar)) {
      data[, var] <- as.vector(data[, var])
    }
    c.forest <- causal_forest(X = data[, covar, drop = FALSE], Y = data[, Y],
                              W = data[, treat], seed = 42)
    att <- average_treatment_effect(c.forest, target.sample = "treated", method = "AIPW")
    att <- c(att, att[1] - 1.96 * att[2], att[1] + 1.96 * att[2])
    return(att)
  }, error = function(e) {
    cat("Error in aipw_grf method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

#### aipw.match()***
aipw.match <- function(data, Y, treat, covar) { # match on ps
  ps <- probability_forest(X = data[, covar], Y = as.factor(data[, treat]), seed = 42)$predictions[,2]
  m.out <- Match(Y = data[, Y], Tr = data[, treat], X = ps,
                 estimand = "ATT", M = 1, replace = FALSE, ties = FALSE, BiasAdjust = FALSE)
  mb <- quiet(MatchBalance(treat ~ ps, data = data, match.out = m.out, nboots= 0))
  ks <- mb$AfterMatching[[1]]$ks$ks$statistic
  s <- data[c(m.out$index.treated, m.out$index.control), ]
  out <- aipw(s, Y, treat, covar)
  # return(out)
  return(c(out, ks))
}

#### aipw_ow()
aipw_ow <- function(data, Y, treat, covar) {
  tryCatch({
    for (var in c(Y, treat, covar)) {
      data[, var] <- as.vector(data[, var])
    }
    X <- data[, covar, drop = FALSE]
    Y <- data[, Y]
    W <- data[, treat]
    # run dml_with_smoother with AIPW_ATT
    dml_fit <- dml_with_smoother(
      Y = Y, D = W, X = X,
      estimators = c("AIPW_ATT"),
      smoother = "honest_forest",
      n_cf_folds = 5,
      n_reps = 1)
    # extract estimate and SE from summary
    summ <- summary(dml_fit, quiet = TRUE)
    est <- summ["AIPW-ATT", "Estimate"]
    se <- summ["AIPW-ATT", "SE"]
    ci_lower <- est - 1.96 * se
    ci_upper <- est + 1.96 * se
    return(c(est, se, ci_lower, ci_upper))
  }, error = function(e) {
    cat("Error in aipw_ow method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

### this script checks for robustness by estimating original model
### using double/debiased machine learning using DoubleML package
#### dml()***
dml <-function(data, Y = NULL, treat = NULL, covar = NULL, clust_var = NULL, ml_l = lrn("regr.cv_glmnet"), ml_m = lrn("regr.cv_glmnet")){
  tryCatch({
    if(is.null(covar)){
      stop("No controls in specification.")
    }
    #require(DoubleML)
    #require(mlr3learners) 
    #require(fixest)
    #require(ggplot2)
    if(is.null(clust_var) == TRUE){
      dat = data[,c(Y,treat,covar)]
      dat = na.omit(dat)
      dml_dat = DoubleMLData$new(
        dat,
        y_col = Y,
        d_cols = treat,
        use_other_treat_as_covariate = FALSE,
        x_cols = covar)
    }else{
      dat = data[,c(Y, treat, covar, clust_var)]
      dat[,clust_var] = as.numeric(factor(dat[,clust_var]))
      dat = dat[is.na(dat[,Y]) == FALSE,]
      dat = dat[is.na(dat[,D]) == FALSE,]
      features = data.frame(model.matrix(formula(paste(c('~ 1',treat,covar), collapse="+")), dat))
      dat = cbind(dat[,c(Y,clust_var)],features)
      dml_dat = DoubleMLClusterData$new(
        dat,
        y_col = Y,
        d_cols = treat,
        cluster_cols = clust_var,
        use_other_treat_as_covariate = FALSE,
        x_cols = covar)
    }
    # set active treatment treatment
    dml_dat$set_data_model(treat)
    # estimate with DML
    set.seed(pi)
    dml_mod = DoubleMLPLR$new(dml_dat, ml_l=ml_l, ml_m=ml_m)
    quiet(dml_mod$fit())
    out = c(dml_mod$coef[treat], dml_mod$se[treat], dml_mod$confint()[treat,])
    return(out)
  }, error = function(e) {
    cat("Error in dml method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

#### estimate_all***
estimate_all <- function(data, Y, treat, covar, 
                         methods = c("diff", "reg", "om.reg", "om.grf",
                                     "matching", "psm", "ipw", "cbps", "ebal", 
                                     "dml", "aipw_grf", "aipw_ow")) {
  
  results <- as.data.frame(matrix(NA, length(methods), 4))
  rownames(results) <- methods
  colnames(results) <- c("Estimate", "SE", "CI_lower", "CI_upper")
  m <- 1
  if ("diff" %in% methods) {
    results[m, ] <- diff(data, Y, treat) 
    m <- m + 1
  }
  if ("reg" %in% methods) {
    results[m, ] <- reg(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("om.reg" %in% methods) {
    results[m, ] <- om.reg(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("om.grf" %in% methods) {
    results[m, ] <- om.grf(data, Y, treat, covar) 
    m <- m + 1
  } 
  if ("matching" %in% methods) {
    results[m, ] <- matching(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("psm" %in% methods) {
    results[m, ] <- psm(data, Y, treat, covar) 
    m <- m + 1
  }  
  if ("ipw" %in% methods) {
    results[m, ] <- ipw(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("cbps" %in% methods) {
    results[m, ] <- cbps(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("ebal" %in% methods) {
    results[m, ] <- quiet(ebal(data, Y, treat, covar))
    m <- m + 1
  }
  # if ("hbal" %in% methods) {
  #   results[m, ] <- quiet(hbal(data, Y, treat, covar))
  #   m <- m + 1
  # }
  if ("dml" %in% methods) {
    results[m, ] <-dml(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("aipw_grf" %in% methods) {
    results[m, ] <- aipw(data, Y, treat, covar) 
    m <- m + 1
  }
  if ("aipw_ow" %in% methods) {
    results[m, ] <- aipw_ow(data, Y, treat, covar)
    m <- m + 1
  }
  return(results)
}

#### plot_coef()***
plot_coef <- function(out, 
                      methods = c("diff", "reg", "om.reg", "om.grf", 
                                  "matching", "psm", "ipw", "cbps", "ebal", 
                                  "dml", "aipw_grf", "aipw_ow"),
                      labels = c("Diff-in-Means", "Reg", "OM: Reg", "OM: GRF",
                                 "NN\nMatching", "PS\nMatching",
                                 "IPW", "CBPS", "Ebal", "DML\nElasnet", "AIPW-GRF", "AIPW-OW"),
                      main = NULL,
                      ylab = "Estimate",
                      band = NULL,
                      line = NULL,
                      grid = TRUE,
                      main.pos = 1,
                      main.line = -2,
                      ylim = NULL,
                      textsize = 0.8
) {
  if (is.null(methods) == TRUE) {
    methods <- rownames(out)
  }
  if (is.null(labels) == TRUE) {
    labels <- methods
  }
  # # check
  # if (is.null(out)==FALSE) {
  #   if (inherits(out, "ivDiag") == FALSE) {stop("\"out\" needs to be a \"ltz\" object.")}
  # }
  # 
  # # title
  # if (is.null(main)==TRUE) {
  #   main <- "Estimates with 95% CIs"
  # }
  # data for the plot
  data <- out
  rg <- range(data[,c(3,4)], na.rm = TRUE)
  adj <- rg[2] - rg[1]
  if (is.null(ylim) == TRUE) {
    ylim  <- c(min(0, rg[1] - 0.3*adj), max(0, rg[2] + 0.35*adj))
  }
  adj2 <- ylim[2] - ylim[1] 
  # Set up the plot
  ncoefs <- length(methods)
  par(mar = c(2.5, 4, 1, 2))
  plot(1: ncoefs, data[, 1], xlim = c(0.5, ncoefs + 0.5), ylim = ylim,
       ylab = "", xlab = "", main = "", 
       axes = FALSE, xaxt = "n", yaxt = "n", type = "n")
  axis(1, at = 1: ncoefs, labels =  labels, las = 1, cex.axis = 0.8)
  axis(2, cex.axis = 0.7)
  mtext(main, main.pos, line = main.line, cex = textsize)
  mtext(ylab, 2, line = 2.5)
  if (is.null(band) == FALSE) {
    rect(-0.5, band[1], ncoefs + 1, band[2], col = "#ff000030", border = "white") # label at bottom
  }
  if (is.null(line) == FALSE) {
    abline(h = line, col = "red", lty = 2)
  }
  if (grid == TRUE) {
    abline(h = axTicks(2), lty = "dotted", col = "gray50")
    abline(v = c(0.5, c(1: ncoefs) + 0.5), lty = "dotted", col = "gray50") # horizontal grid
  }
  abline(h = 0, col = "red", lwd = 2, lty = "solid")
  segments(y0 = data[, 3], x0 = c(1: ncoefs), y1 = data[, 4], x1 = c(1: ncoefs), lwd = 2) #CI
  points(1: ncoefs, data[, 1], pch = 16, col = 1, cex = 1.2) #point coefs
  box()
}

#### save_att_panels()
save_att_panels <- function(
    out_list, plot_titles, band_list, est_list,
    prefix = "ldw_model_a_plus",
    plots_per_page = 4,
    ylim = c(-15500, 5500),
    folder = "graphs/lalonde") {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n <- length(out_list)
  pages <- ceiling(n / plots_per_page)
  
  for (p in seq_len(pages)) {
    start <- (p - 1) * plots_per_page + 1
    end <- min(p * plots_per_page, n)
    file_name <- file.path(folder, paste0(prefix, "_", p, ".pdf"))
    pdf(file_name, width = 8, height = 11)
    par(mfrow = c(plots_per_page, 1), mar = c(3, 4, 3, 2))
    for (i in start:end) {
      plot_coef(out_list[[i]], band = band_list[[i]], line = est_list[[i]],
                ylim = ylim, main = plot_titles[i])
    }
    dev.off()
  }
}

#### create_matrix_results()
create_matrix_results <- function(all_outs, all_out_mat, sample_names) {
  n_samples <- length(sample_names)
  n_estimators <- nrow(all_outs[[1]])
  result_mat <- matrix("", nrow = n_estimators + 1, ncol = n_samples * 2)
  cnames <- character(n_samples * 2)
  for (j in seq_along(sample_names)) {
    cnames[(j - 1) * 2 + 1] <- sample_names[j]
    cnames[(j - 1) * 2 + 2] <- ""
  }
  colnames(result_mat) <- cnames
  estimator_names <- rownames(all_outs[[1]])
  rownames(result_mat) <- c("Experimental Benchmark", estimator_names)
  # all_out_mat contains 13 entries
  bench.exp      <- if(!is.null(all_out_mat[[1]])) all_out_mat[[1]][1, 1:2] else c(NA, NA)
  bench.cps      <- if(length(all_out_mat) >= 6) lapply(all_out_mat[2:6], function(x) if(!is.null(x)) x[1, 1:2] else c(NA, NA)) else NULL
  bench.psid     <- if(length(all_out_mat) >= 11) lapply(all_out_mat[7:11], function(x) if(!is.null(x)) x[1, 1:2] else c(NA, NA)) else NULL
  bench.cps_plus <- if(length(all_out_mat) >= 12 && !is.null(all_out_mat[[12]])) all_out_mat[[12]][1, 1:2] else c(NA, NA)
  bench.psid_plus<- if(length(all_out_mat) >= 13 && !is.null(all_out_mat[[13]])) all_out_mat[[13]][1, 1:2] else c(NA, NA)
  # all_outs contains 15 entries
  for (j in seq_along(all_outs)) {
    out <- all_outs[[j]]
    # assign benchmark according to position
    # use bench.exp for the first 3 entries in all_outs
    if (j <= 3 && !any(is.na(bench.exp))) {
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.exp[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.exp[2]), ")")
    }
    # use bench.cps for the next 5 entries in all_outs
    # use all entries in bench.cps one after another
    if (j >= 4 && j <= 8 && !is.null(bench.cps) && length(bench.cps) >= (j - 3)) {
      b <- bench.cps[[j - 3]]
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", b[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", b[2]), ")")
    }
    # use bench.psid for the next 5 entries in all_outs
    # use all entries in bench.psid one after another
    if (j >= 9 && j <= 13 && !is.null(bench.psid) && length(bench.psid) >= (j - 8)) {
      b <- bench.psid[[j - 8]]
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", b[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", b[2]), ")")
    }
    # use bench.cps_plus for the second last entry in all_outs
    if (j == n_samples - 1 && !any(is.na(bench.cps_plus))) {
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.cps_plus[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.cps_plus[2]), ")")
    }
    # use bench.psid_plus for the last entry in all_outs
    if (j == n_samples && !any(is.na(bench.psid_plus))) {
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.psid_plus[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.psid_plus[2]), ")")
    }
    # fill in estimates + SEs
    for (i in 2:(n_estimators + 1)) {
      if ((i-1) <= nrow(out)) {
        result_mat[i, (j - 1) * 2 + 1] <- sprintf("%.2f", out[i - 1, 1])
        result_mat[i, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", out[i - 1, 2]), ")")
      } else {
        result_mat[i, (j - 1) * 2 + 1] <- ""
        result_mat[i, (j - 1) * 2 + 2] <- ""
      }
    }
  }
  return(result_mat)
}

#### eval_att()
eval_att <- function(result) {
  data.frame(
    Mean_SE = mean(result[, "SE"], na.rm = TRUE), # mean standard error
    Min_Estimate = min(result[, "Estimate"], na.rm = TRUE), # maximum estimate 
    Max_Estimate = max(result[, "Estimate"], na.rm = TRUE), # minimum estimate 
    Diff_Estimate = max(result[, "Estimate"], na.rm = TRUE) - min(result[, "Estimate"], na.rm = TRUE)
  )
}

## 5.2 CATT
#### plot_catt_panels()***
plot_catt_panels <- function(exp_catt_list, catt_list, plot_titles, plots_per_page = 4, range = c(-8000, 8000)) {
  n <- length(catt_list)
  num_pages <- ceiling(n / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx   <- min(page * plots_per_page, n)
    par(mfrow = c(2, 2), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
      catt_ldw <- exp_catt_list[[i]]$catt
      att_ldw  <- exp_catt_list[[i]]$att[1]
      id_ldw   <- if (!is.null(exp_catt_list[[i]]$id)) exp_catt_list[[i]]$id else seq_along(catt_ldw)
      catt2 <- catt_list[[i]]$catt
      att2  <- catt_list[[i]]$att[1]
      id2   <- if (!is.null(catt_list[[i]]$id)) catt_list[[i]]$id else seq_along(catt2)
      common_ids <- intersect(id_ldw, id2)
      idx_ldw    <- match(common_ids, id_ldw)
      idx_other  <- match(common_ids, id2)
      catt1_plot <- catt_ldw[idx_ldw]
      catt2_plot <- catt2[idx_other]
      plot_catt(
        catt1 = catt1_plot,
        catt2 = catt2_plot,
        att1  = att_ldw,
        att2  = att2,
        xlab  = "CATT (Experimental)",
        ylab  = plot_titles[i],
        main  = plot_titles[i],
        axes.range = range
      )
    }
  }
}

#### save_catt_panels()***
save_catt_panels <- function(
    exp_catt_list,   
    catt_list, 
    plot_titles, 
    prefix = "catt_top5", 
    plots_per_page = 1, 
    range = c(-8000, 8000), 
    folder = "graphs/lalonde") {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n <- length(catt_list)
  num_pages <- ceiling(n / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx   <- min(page * plots_per_page, n)
    file_name <- file.path(folder, paste0(prefix, "_", page, ".pdf"))
    pdf(file = file_name, width = 10, height = 12)
    par(mfrow = c(plots_per_page, 1), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
      # Use matching elements from both lists
      catt_ldw <- exp_catt_list[[i]]$catt
      att_ldw  <- exp_catt_list[[i]]$att[1]
      id_ldw   <- if (!is.null(exp_catt_list[[i]]$id)) exp_catt_list[[i]]$id else seq_along(catt_ldw)
      catt2 <- catt_list[[i]]$catt
      att2  <- catt_list[[i]]$att[1]
      id2   <- if (!is.null(catt_list[[i]]$id)) catt_list[[i]]$id else seq_along(catt2)
      common_ids <- intersect(id_ldw, id2)
      idx_ldw    <- match(common_ids, id_ldw)
      idx_other  <- match(common_ids, id2)
      catt1_plot <- catt_ldw[idx_ldw]
      catt2_plot <- catt2[idx_other]
      plot_catt(
        catt1 = catt1_plot,
        catt2 = catt2_plot,
        att1  = att_ldw,
        att2  = att2,
        xlab  = "CATT (Experimental)",
        ylab  = plot_titles[i],
        main  = plot_titles[i],
        axes.range = range
      )
    }
    # fill empty panels if last page is not full
    if ((end_idx - start_idx + 1) < plots_per_page) {
      for (k in seq_len(plots_per_page - (end_idx - start_idx + 1))) plot.new()
    }
    dev.off()
  }
}

#### save_main_catt_panels()***
save_main_catt_panels <- function(
    catt_refs,
    catt_comps, 
    ylabels,
    prefix = "catt_main", 
    plots_per_page = 1, 
    main_titles = NULL, 
    range = c(-8000, 8000), 
    folder = "graphs/lalonde"
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n_panels <- length(catt_comps)
  num_pages <- ceiling(n_panels / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_panels)
    file_name <- file.path(folder, paste0(prefix, "_", page, ".pdf"))
    pdf(file = file_name, width = 10, height = 3 * plots_per_page)
    par(mfrow = c(plots_per_page, 1), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
      ref_idx <- min(i, length(catt_refs))
      catt1 <- catt_refs[[ref_idx]]$catt
      catt2 <- catt_comps[[i]]$catt
      att1  <- catt_refs[[ref_idx]]$att[1]
      att2  <- catt_comps[[i]]$att[1]
      plot_catt(
        catt1 = catt1,
        catt2 = catt2,
        att1  = att1,
        att2  = att2,
        xlab  = "CATT (Reference)",
        ylab  = ylabels[i],
        main  = if (is.null(main_titles)) paste("Comparison", i) else main_titles[i],
        axes.range = range
      )
    }
    # fill empty panels if last page not full
    if ((end_idx - start_idx + 1) < plots_per_page) {
      for (k in seq_len(plots_per_page - (end_idx - start_idx + 1))) plot.new()
    }
    dev.off()
  }
}

#### save_plus_catt_panels()***
save_plus_catt_panels <- function(
    catt1_list, 
    catt2_list, 
    ylabels, 
    prefix = "catt_plus", 
    plots_per_page = 1, 
    main_titles = NULL, 
    range = c(-8000, 8000), 
    folder = "graphs/lalonde"
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n_panels <- min(length(catt1_list), length(catt2_list))
  num_pages <- ceiling(n_panels / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_panels)
    file_name <- file.path(folder, paste0(prefix, "_", page, ".pdf"))
    pdf(file = file_name, width = 10, height = 3 * plots_per_page)
    par(mfrow = c(plots_per_page, 1), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
      catt1 <- catt1_list[[i]]$catt
      catt2 <- catt2_list[[i]]$catt
      att1  <- catt1_list[[i]]$att[1]
      att2  <- catt2_list[[i]]$att[1]
      plot_catt(
        catt1 = catt1,
        catt2 = catt2,
        att1  = att1,
        att2  = att2,
        xlab  = "CATT (Reference/Trimmed)",
        ylab  = ylabels[i],
        main  = if (is.null(main_titles)) paste("Trimmed Panel", i) else main_titles[i],
        axes.range = range
      )
    }
    # fill empty panels if last page not full
    if ((end_idx - start_idx + 1) < plots_per_page) {
      for (k in seq_len(plots_per_page - (end_idx - start_idx + 1))) plot.new()
    }
    dev.off()
  }
}

#### eval_catt()
eval_catt <- function(all_catt, plot_titles) {
  do.call(rbind, lapply(seq_along(all_catt), function(i) {
    catt_vec <- all_catt[[i]]$catt
    data.frame(
      Method = plot_titles[i],
      Min_Catt = min(catt_vec, na.rm = TRUE),
      Max_Catt = max(catt_vec, na.rm = TRUE),
      Mean_Catt = mean(catt_vec, na.rm = TRUE),
      Diff_Catt = max(catt_vec, na.rm = TRUE) - min(catt_vec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

## 5.3 QTET
#### est_qte_safe()***
est_qte_safe <- function(Y, treat, covar, data, cores = 1) {
  tryCatch({
    est_qte(Y, treat, covar, data = data, cores = cores)
  }, error = function(e) {
    cat("Warning in est_qte:", e$message, "\n")
    return(NULL)
  })
}

#### plot_qte_top()***
plot_qte_top <- function(qtet_top, qtet_top0, bm_list, plot_titles, main_start = 1, 
                         ylim = c(-25000, 15000), col = NULL) {
  n <- length(qtet_top)
  for (i in seq_len(n)) {
    main_title <- plot_titles[main_start + i - 1]
    mod <- qtet_top[[i]]
    mod2 <- qtet_top0[[i]]
    bm  <- bm_list[[i]]
    plot_qte(mod, mod2, bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
  }
}

#### save_qtet()***
save_qtet <- function(plots, plot_titles = NULL, main_start = 1, 
                      ylim = c(-25000, 15000), col = NULL, prefix = "ldw") {
  dir.create("graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  n <- length(plots)
  for (i in seq_len(n)) {
    p <- plots[[i]]
    main_title <- if (is.null(plot_titles)) p$main else plot_titles[main_start + i - 1]
    clean_title <- gsub("[^a-zA-Z0-9]", "_", main_title)
    file_name <- sprintf("graphs/lalonde/%s_qtet_estimates_%s.pdf", prefix, clean_title)
    pdf(file = file_name, width = 7, height = 5)
    plot_qte(p$mod, p$mod0, p$bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
    dev.off()
  }
}

#### save_qte_top()***
save_qte_top <- function(qtet_top, qtet_top0, bm_list, plot_titles, main_start = 1,
                         ylim = c(-25000, 15000), col = NULL, prefix = "ldw_top") {
  n <- length(qtet_top)
  dir.create("graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(n)) {
    mod <- qtet_top[[i]]
    mod2 <- qtet_top0[[i]]
    bm <- bm_list[[i]]                 # Use the i-th entry in bm_list
    main_title <- plot_titles[main_start + i - 1]
    clean_title <- gsub("[^a-zA-Z0-9]", "_", main_title)
    file_name <- sprintf("graphs/lalonde/%s_qte_estimates_%s.pdf", prefix, clean_title)
    pdf(file = file_name, width = 7, height = 5)
    plot_qte(mod, mod2, bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
    dev.off()
  }
}

## 5.4 Assessing outcome weights (OW)
#### get_res_att()
get_res_att <- function(dataset_list, Y, treat, covar,
                        estimator = "AIPW_ATT", 
                        smoother = "honest_forest", 
                        n_cf_folds = 5,
                        n_reps = 1) {
      results <- lapply(seq_along(dataset_list),
          function(i) {
            data <- dataset_list[[i]]
            n_obs <- nrow(data)
            # adjust n_cf_folds for small samples
            folds <- if (n_obs < 100) {
              max(2, min(3, n_cf_folds))
            } else if (n_obs < 300) {
              max(3, min(4, n_cf_folds))
            } else {
              n_cf_folds
            }
            tryCatch({
              dml_with_smoother(
                Y = data[[Y]],
                D = data[[treat]],
                X = data[, covar, drop = FALSE],
                estimator = estimator,
                smoother = smoother,
                n_cf_folds = folds,
                n_reps = n_reps
              )
            }, error = function(e) {
              cat(sprintf("WARNING: get_res_att failed for dataset %d (n=%d): %s\n", i, n_obs, e$message))
              return(NULL)
            })
          }
      )
    # summary report
    null_count <- sum(sapply(results, is.null))
    success_count <- length(results) - null_count
    cat(sprintf("OUTCOME WEIGHTS SUMMARY: %d succeeded, %d failed out of %d datasets\n", 
                success_count, null_count, length(results)))
    return(results)
}

#### derive_ow()
derive_ow <- function(results_list) {
lapply(results_list, function(res) get_outcome_weights(res))
}

#### plot_ow()
plot_ow <- function(outcome_weights, plot_titles = NULL, breaks = 50, 
                    col = "#ff000080", xlab = "Outcome Weight", 
                    estimator = "AIPW-ATT") {
  N <- length(outcome_weights)
  for (i in seq_len(N)) {
    weights <- outcome_weights[[i]]$omega[estimator, ]
    main_title <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    hist(weights, breaks = breaks, main = main_title, xlab = xlab, col = col)
    mtext(paste("N =", length(weights)), side = 3, line = -1.5, cex = 0.8)
  }
  par(mfrow = c(1, 1))
}

#### eval_ow()
eval_ow <- function(outcome_weights, dataset_list, plot_titles = NULL, treat = "treat", estimator = "AIPW-ATT") {
  results <- lapply(seq_along(outcome_weights), function(i) {
    ow <- outcome_weights[[i]]$omega[estimator, ]
    treat <- dataset_list[[i]][[treat]]
    method <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    sum_treated <- sum(ow[treat == 1])
    sum_untreated <- sum(ow[treat == 0])
    data.frame(
      Method = method,
      Sum_Treated = sum_treated,
      Sum_Untreated = sum_untreated,
      Total_Sum = sum_treated + sum_untreated,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#### save_ow()
save_ow <- function(outcome_weights, plot_titles = NULL,
                    breaks = 50, col = "#ff000080", xlab = "Outcome Weight",
                    prefix = "model_a", estimand = "AIPW-ATT") {
  dir.create("graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  N <- length(outcome_weights)
  for (i in seq_len(N)) {
    file_name <- sprintf("graphs/lalonde/%s_outcomewt_%d.pdf", prefix, i)
    pdf(file = file_name, width = 8, height = 6)
    weights <- outcome_weights[[i]]$omega[estimand, ]
    main_title <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    hist(weights, breaks = breaks, main = main_title, xlab = xlab, col = col)
    mtext(paste("N =", length(weights)), side = 3, line = -1.5, cex = 0.8)
    dev.off()
  }
}

# 6. Sensitivity Analysis
#### check_filter_datasets()
check_filter_datasets <- function(datasets, Y, treat, covar, bm) {
  valid_datasets <- list()
  for (i in seq_along(datasets)) {
    data <- datasets[[i]]
    name <- names(datasets)[i]
    if (is.null(name) || name == "") name <- paste0("dataset_", i)  # provide a name if missing
    vars_needed <- c(Y, treat, covar, bm)
    if (!all(vars_needed %in% names(data))) { # check all variables exist
      message("Removed ", name, ": missing required variables") 
      next
    }
    sub <- data[, vars_needed, drop = FALSE]
    if (any(is.na(sub))) { # check no missing values
      message("Removed ", name, ": contains missing values") 
      next
    }
    tvals <- unique(sub[[treat]])  # check treatment binary
    if (!all(tvals %in% c(0, 1)) && !all(tvals %in% c(TRUE, FALSE))) {
      message("Removed ", name, ": treatment variable not binary")
      next
    }
    if (any(sapply(sub[bm], function(x) length(unique(x)) <= 1))) { # check baseline variables have variation
      message("Removed ", name, ": baseline variable(s) lack variation") 
      next
    }
    if (any(sapply(sub[covar], function(x) length(unique(x)) <= 1))) { # check covariates have variation
      message("Removed ", name, ": covariate variable(s) lack variation")
      next
    }
    valid_datasets[[name]] <- data # if all passed dataset is kept
  }
  return(valid_datasets)
}

#### save_sensitivity_plots()
save_sensitivity_plots <- function(filtered_datasets, Y, treat, covar, bm, plot_titles, prefix) {
  folder <- "graphs/lalonde"
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  success_count <- 0
  fail_count <- 0
  for (i in seq_along(filtered_datasets)) {
    idx <- i
    file_name <- file.path(folder, paste0(prefix, "_sensitivity_", idx, ".pdf"))
    pdf(file_name, width = 8, height = 8)
    tryCatch({
      sens_ana(filtered_datasets[[i]], Y, treat, covar, bm, kd = 1:3)
      if (!is.null(plot_titles) && length(plot_titles) >= idx) {
        title(main = plot_titles[idx])
      }
      success_count <- success_count + 1
    }, error = function(e) {
      plot.new()
      title_text <- if (!is.null(plot_titles) && length(plot_titles) >= idx) plot_titles[idx] else paste("Dataset", idx)
      text(0.5, 0.5, paste("ERROR:", title_text, "\n", e$message), cex = 0.8)
      cat(sprintf("WARNING: Sensitivity plot %d (%s) failed: %s\n", idx, title_text, e$message))
      fail_count <- fail_count + 1
    })
    dev.off()
  }
  cat(sprintf("SUMMARY: Sensitivity plots - %d succeeded, %d failed out of %d total\n", 
              success_count, fail_count, length(filtered_datasets)))
}

# 7. Balance
#### save_balance()
save_balance <- function(
    datasets,
    method_names,
    balance_var = "re75",
    prefix = "balance_panels",
    plots_per_page = 5,
    folder = "graphs/lalonde"
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n <- length(datasets)
  pages <- ceiling(n / plots_per_page)
  for (p in seq_len(pages)) {
    start <- (p - 1) * plots_per_page + 1
    end <- min(p * plots_per_page, n)
    file_name <- file.path(folder, paste0(prefix, "_", balance_var, "_", p, ".pdf"))
    plot_list <- list()
    for (i in start:end) {
      plot_title <- paste0("(", LETTERS[i], ") ", method_names[i], " - ", balance_var)
      pl <- tryCatch({
        form <- as.formula(paste0("treat ~ ", balance_var))
        bal.plot(form, data = datasets[[i]], which = "both") +
          ggtitle(plot_title)
      }, error = function(e) {
        message(sprintf("Could not plot %s: %s", method_names[i], e$message))
        ggplot() + ggtitle(sprintf("Error: %s", method_names[i]))
      })
      plot_list[[length(plot_list)+1]] <- pl
    }
    pdf(file_name, width = 8, height = 11)
    print(wrap_plots(plotlist = plot_list, ncol = 1))
    dev.off()
  }
}
