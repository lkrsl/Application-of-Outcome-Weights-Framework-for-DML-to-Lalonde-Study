# 1 Set up
# 1.1 Installation
packages <- c("data.table", "dplyr", "ggplot2", "gridExtra", "highr", "MatchIt", "optmatch", "quickmatch", 
  "readr", "rgenoud", "tidyr", "tidyverse", "WeightIt"
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
library(cobalt)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(highr)
library(MatchIt)
library(optmatch)
library(quickmatch)
library(readr)
library(tidyr)
library(tidyverse)
library(WeightIt)

#### source files for aipw_ow extension
skripte <- list.files(
  "../package/OutcomeWeights/R",  
  pattern = "\\.R$", 
  full.names = TRUE)
for (f in skripte) source(f, local = knitr::knit_global())

# 1.2 Data inspection
#### inspect_datasets()
inspect_data <- function(data_list) {
  if (is.data.frame(data_list)) {
    data_list <- list(dataset = data_list)
  }
  data.frame(
    dataset = names(data_list),
    num_obs = sapply(data_list, nrow),
    num_vars = sapply(data_list, ncol),
    name_vars = sapply(data_list, function(df) paste(names(df), collapse = ", ")),
    row.names = NULL
  )
}

# 2. Improving covariate balance and overlap
## 2.1 Matching
### 2.1.1 Profile matching
#### matchit_profile()
matchit_profile <- function(data, treat, covar) {
  overall_means <- colMeans(data[, covar], na.rm = TRUE) 
  cov_matrix <- as.matrix(data[, covar])   
  dist_from_target <- apply(cov_matrix, 1, function(x) sqrt(sum((x - overall_means)^2)))
  data$dist_from_target <- dist_from_target 
  formula <- as.formula(paste(treat, "~", paste(c(covar, "dist_from_target"), collapse = "+")))
  match_out <- matchit(formula, data = data, method = "nearest", distance = "glm")
  return(match_out)
}

## 2.2 Weighting 
### 2.2.1 Standardized mortality ratio (SMR) treated weights
#### create_smr_weights()
create_smr_weights <- function(data, formula, estimand = "ATT") {
  ps_model <- glm(formula, data = data, family = binomial())
  ps <- predict(ps_model, type = "response")
  if (estimand == "ATT") {
    weights <- ifelse(data$treat == 1, 1/ps, 1/(1-ps)) 
  } else if (estimand == "ATE") {
    weights <- ifelse(data$treat == 1, 1/ps, 1/(1-ps)) 
  }
  return(weights)
}

#### create_overlap_weights()
create_overlap_weights <- function(data, formula) {
  ps_model <- glm(formula, data = data, family = binomial())
  ps <- predict(ps_model, type = "response")
  return(ifelse(data$treat == 1, 1 - ps, ps))
}

## 2.3 Truncation
### 2.3.1 Fixed maximum value truncation
#### truncate_weights_fixed()
truncate_weights_fixed <- function(data, weight_col, max_weight = 10) {
  data[[weight_col]] <- pmin(data[[weight_col]], max_weight)
  return(data)
}

### 2.3.2 At percentile truncation
#### truncate_weights_percentile()
truncate_weights_percentile <- function(data, weight_col, percentile = 0.99) {
  weight_cutoff <- quantile(data[[weight_col]], probs = percentile, na.rm = TRUE)
  data[[weight_col]] <- pmin(data[[weight_col]], weight_cutoff)
  return(data)
}  

### 2.3.3 Adaptive weight truncation
#### check_weights()
check_weights <- function(data, weight_col = "weight") {
  w <- data[[weight_col]]
  cat("Variance of", weight_col, ":", var(w, na.rm = TRUE), "\n")
}  

#### truncate_weights_adaptive()
truncate_weights_adaptive <- function(data, weight_col, c = 3) {
  w <- data[[weight_col]]
  # Only apply if variance is positive
  if (var(w, na.rm = TRUE) > 0) {
    cutoff <- mean(w, na.rm = TRUE) + c * sd(w, na.rm = TRUE)
    data[[weight_col]] <- pmin(w, cutoff)
  }
  return(data)
}

## 2.4 Trimming
### 2.4.1 Propensity score threshold trimming (similar to tutorial of Imbens & Xu (2024))
ps_trim <- function(data, ps = "ps_assoverlap", threshold = 0.9) { 
  sub <- data[which(data[, ps] < threshold), ]
  return(sub)
}

### 2.4.2 Common range trimming
#### common_range_trim ()
common_range_trim <- function(data, ps_col = "ps_assoverlap", treat_col = "treat") {
  lower_cut <- max(
    min(data[[ps_col]][data[[treat_col]] == 1], na.rm = TRUE),
    min(data[[ps_col]][data[[treat_col]] == 0], na.rm = TRUE)
  )
  upper_cut <- min(
    max(data[[ps_col]][data[[treat_col]] == 1], na.rm = TRUE),
    max(data[[ps_col]][data[[treat_col]] == 0], na.rm = TRUE)
  )
  sub <- data[data[[ps_col]] >= lower_cut & data[[ps_col]] <= upper_cut, ]
  return(sub)
}

### 2.4.3 Crump trimming
#### crump_trim ()
crump_trim <- function(data, ps_col = "ps_assoverlap", lower = 0.1, upper = 0.9) {
  sub <- data[data[[ps_col]] >= lower & data[[ps_col]] <= upper, ]
  return(sub)
}

### 2.4.4 Stuermer trimming
#### stuermer_trim ()
stuermer_trim <- function(data, treat_col = "treat", ps_col = "ps_assoverlap", 
                         lower_percentile = 0.05, upper_percentile = 0.95) {
  treated_ps   <- data[[ps_col]][data[[treat_col]] == 1]
  untreated_ps <- data[[ps_col]][data[[treat_col]] == 0]
  lower_cutoff <- quantile(treated_ps, probs = lower_percentile, na.rm = TRUE)
  upper_cutoff <- quantile(untreated_ps, probs = upper_percentile, na.rm = TRUE)
  sub <- data[data[[ps_col]] >= lower_cutoff & data[[ps_col]] <= upper_cutoff, ]
  return(sub)
}

### 2.4.5 Walker trimming
#### walker_trim ()
walker_trim <- function(data, treat_col = "treat", ps_col = "ps_assoverlap", 
                        lower_cutoff = 0.3, upper_cutoff = 0.7) {
  treat_prevalence  <- mean(data[[treat_col]], na.rm = TRUE)
  logit_ps          <- log(data[[ps_col]] / (1 - data[[ps_col]]))
  logit_prevalence  <- log(treat_prevalence / (1 - treat_prevalence))
  preference_score  <- 1 / (1 + exp(-(logit_ps - logit_prevalence)))
  sub <- data[preference_score >= lower_cutoff & preference_score <= upper_cutoff, ]
  return(sub)
}

## 2.5 Combination of methods
#### trim_attach_weights()
trim_attach_weights <- function(trimmed_data, original_data, weight_col){
trimmed_data$orig_row   <- as.integer(rownames(trimmed_data))
original_data$orig_row  <- as.integer(rownames(original_data))
trimmed_data <- merge(
  trimmed_data,
  original_data[, c("orig_row", weight_col)],
  by = "orig_row",
  all.x = TRUE
)
colnames(trimmed_data)[colnames(trimmed_data) == weight_col] <- "weight"
trimmed_data$orig_row <- NULL
original_data$orig_row <- NULL
return(trimmed_data)
}

# 3. Reassessing methods
## 3.1 Matching
#### get_smd_stats()
get_smd_stats <- function(match_object) {
  bal <- bal.tab(match_object, stats = "mean.diffs", un = TRUE, s.d.denom = "treated")
  smds <- bal$Balance$Diff.Adj
  smds <- smds[!(rownames(bal$Balance) %in% c("distance"))]
  mean_smd <- mean(abs(smds), na.rm = TRUE)
  max_smd <- max(abs(smds), na.rm = TRUE)
  return(c(Mean_Abs_SMD = mean_smd, Max_Abs_SMD = max_smd))
}

### 3.1.1 SMD
#### compute_abs_smd_matchit()
compute_abs_smd_matchit <- function(match_list) {
  smd_list <- lapply(match_list, get_smd_stats)       
  smd_mat <- do.call(rbind, smd_list)                 
  smd_df <- data.frame(
    Method = names(match_list),
    Mean_Abs_SMD = smd_mat[, "Mean_Abs_SMD"],
    Max_Abs_SMD  = smd_mat[, "Max_Abs_SMD"],
    row.names = NULL
  )
  return(smd_df)
}  

### 3.1.2 ESS
#### compute_ess_matchit()
compute_ess_matchit <- function(bal_tab_object) {
  samples <- bal_tab_object$Observations
  df <- as.data.frame(samples)
  df$Method <- rownames(samples)
  df <- df[, c("Method","Control", "Treated")]
  rownames(df) <- NULL
  df
}

### 3.1.3 Visuals
#### plot_matchit()
plot_matchit <- function(match_list, dataset_name) {
  for (method_name in names(match_list)) {
    match_obj <- match_list[[method_name]]
    if (method_name == "subcl") {
      plot(summary(match_obj, subclass = TRUE), 
           main = paste0(dataset_name, " - ", method_name," matching"))
    } else {
      plot(summary(match_obj), 
           main = paste0(dataset_name, " - ", method_name," matching"))
    }
  }
} 

## 3.2 Weighting
### 3.2.1 SMD
#### compute_abs_smd_weight()
compute_abs_smd_weight <- function(data, treat, covar, weight_cols) {
  smd_list <- lapply(weight_cols, function(wcol) {
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = data[[wcol]],
      un = TRUE,
      s.d.denom = "treated"
    )
    smd_df <- as.data.frame(bal_obj$Balance)
    smd_vals <- abs(smd_df$Diff.Adj)
    mean_smd <- mean(smd_vals, na.rm = TRUE)
    max_smd  <- max(smd_vals, na.rm = TRUE)
    return(data.frame(
      Method = wcol,
      Mean_Abs_SMD = mean_smd,
      Max_Abs_SMD  = max_smd
    ))
  })
  do.call(rbind, smd_list)
}

### 3.2.2 ESS
#### compute_ess_weight()
compute_ess_weight <- function(data, treat, covar, weight_cols) {
  ess_list <- lapply(weight_cols, function(wcol) {
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = data[[wcol]],
      un = FALSE
    )
    samples <- bal_obj$Observations
    if ("Adjusted" %in% rownames(samples)) {
      df <- samples["Adjusted", c("Control", "Treated"), drop = FALSE]
    } else {
      df <- samples[1, c("Control", "Treated"), drop = FALSE]
    }
    df <- cbind(Method = wcol, df)
    rownames(df) <- NULL
    return(df)
  })
  do.call(rbind, ess_list)
}

### 3.2.3 Visuals
#### plot_weighting_methods()
plot_weighting_methods <- function(data, treat, covar, weight_list, dataset_name = NULL) {
  for (wcol in names(weight_list)) {
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = weight_list[[wcol]],
      un = TRUE,
      s.d.denom = "treated"
    )
    title_text <- if (!is.null(dataset_name)) paste(dataset_name, "-", wcol, "weighting") else wcol
    lp <- cobalt::love.plot(
      bal_obj,
      stats = "mean.diffs",
      abs = TRUE,
      var.order = "unadjusted",
      thresholds = c(m = 0.1),
      title = title_text
    )
    print(lp)
  }
}

## 3.3 Truncation
### 3.3.1 SMD
#### compute_abs_smd_trunc()
compute_abs_smd_trunc <- function(trunc_list, treat, covar, weight_cols) {
  all_smd <- list()
  for(trunc_name in names(trunc_list)) {
    dataset <- trunc_list[[trunc_name]]
    smd_list <- lapply(weight_cols, function(wcol) {
      if (wcol %in% names(dataset)) {
        bal_obj <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = dataset,
          weights = dataset[[wcol]],
          un = TRUE,
          s.d.denom = "treated"
        )
        smd_df <- as.data.frame(bal_obj$Balance)
        smd_vals <- abs(smd_df$Diff.Adj)
        data.frame(
          Method = paste(trunc_name, wcol, sep = "_"),
          Mean_Abs_SMD = mean(smd_vals, na.rm = TRUE),
          Max_Abs_SMD = max(smd_vals, na.rm = TRUE)
        )
      } else {
        NULL
      }
    })
    all_smd[[trunc_name]] <- do.call(rbind, smd_list)
  }
  smd_summary <- do.call(rbind, all_smd)
  rownames(smd_summary) <- NULL  
  return(smd_summary)
}

### 3.3.2 ESS
#### compute_ess_trunc()
compute_ess_trunc <- function(trunc_list, treat, covar, weight_cols) {
  all_ess <- list()
  for(trunc_name in names(trunc_list)) {
    dataset <- trunc_list[[trunc_name]]
    ess_list <- lapply(weight_cols, function(wcol) {
      if (wcol %in% names(dataset)) {
        bal_obj <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = dataset,
          weights = dataset[[wcol]],
          un = FALSE,
          s.d.denom = "treated"
        )
        samples <- bal_obj$Observations
        if ("Adjusted" %in% rownames(samples)) {
          df <- data.frame(
            Method = paste(trunc_name, wcol, sep = "_"),
            Control = samples["Adjusted", "Control"],
            Treated = samples["Adjusted", "Treated"]
          )
        } else {
          df <- data.frame(
            Method = paste(trunc_name, wcol, sep = "_"),
            Control = samples[1, "Control"],
            Treated = samples[1, "Treated"]
          )
        }
        df
      } else {
        NULL
      }
    })
    all_ess[[trunc_name]] <- do.call(rbind, ess_list)
  }
  ess_summary <- do.call(rbind, all_ess)
  rownames(ess_summary) <- NULL   
  return(ess_summary)
}

### 3.3.3 Visuals
#### plot_trunc_methods()
plot_trunc_methods <- function(trunc_list, treat, covar, weight_cols, dataset_name = NULL) {
  for(trunc_name in names(trunc_list)) {
    dataset <- trunc_list[[trunc_name]]
    for(wcol in weight_cols) {
      if(wcol %in% names(dataset)) {
        bal_obj <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = dataset,
          weights = dataset[[wcol]],
          un = TRUE,
          s.d.denom = "treated"
        )
        title_text <- if (!is.null(dataset_name)) {
          paste(dataset_name, "-", trunc_name, "-", wcol, "truncation")
        } else {
          paste(trunc_name, wcol, "truncation")
        }
        lp <- cobalt::love.plot(
          bal_obj,
          stats = "mean.diffs",
          abs = TRUE,
          var.order = "unadjusted",
          thresholds = c(m = 0.1),
          stars = "none",
          title = title_text
        )
        print(lp)
      }
    }
  }
}

## 3.4 Trimming
### 3.4.1 SMD
#### compute_abs_smd_trim()
compute_abs_smd_trim <- function(trimming_list, treat, covar) {
  smd_list <- lapply(names(trimming_list), function(name) {
    data <- trimming_list[[name]]
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      un = TRUE, 
      s.d.denom = "treated"
    )
    smd_df <- as.data.frame(bal_obj$Balance)
    smd_vals <- abs(smd_df$Diff.Un)
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

### 3.4.2 ESS
#### compute_ess_trim()
compute_ess_trim <- function(trimming_list, treat, covar) {
  ess_list <- lapply(trimming_list, function(data) {
    bal_obj <- cobalt::bal.tab(as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
                               data = data, un = TRUE)
    samples <- bal_obj$Observations
    df <- as.data.frame(samples)[c("Control", "Treated")]
    return(df)
  })
  ess_df <- do.call(rbind, ess_list)
  ess_df$Method <- names(trimming_list)
  rownames(ess_df) <- NULL
  ess_df <- ess_df[, c("Method", "Control", "Treated")]
  return(ess_df)
}

### 3.4.3 Visuals
#### plot_trim()
plot_trim <- function(data_list, treat, covar) {
  par(mfrow = c(2,3))
  for (data_obj in data_list) {
    assess_overlap(data_obj, treat = treat, cov = covar)
  }
}

## 3.5 Combined methods
### 3.5.1 SMD
#### compute_smd_all_datasets()
compute_smd_all_datasets <- function(combined_list, treat, covar) {
  smd_list <- lapply(names(combined_list), function(weight_method) {
    method_list <- combined_list[[weight_method]]
    res <- lapply(names(method_list), function(trim_method) {
      data <- method_list[[trim_method]]
      # Use equal weights if weight column missing
      if (!"weight" %in% colnames(data)) {
        data$weight <- rep(1, nrow(data))
      }
      bal_obj <- cobalt::bal.tab(
        as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
        data = data,
        weights = data$weight,
        un = TRUE,
        s.d.denom = "treated"
      )
      smd_df <- as.data.frame(bal_obj$Balance)
      smd_vals <- abs(smd_df$Diff.Adj)
      method_name <- paste(weight_method, trim_method, sep = "_")
      data.frame(
        Method = method_name,
        Mean_Abs_SMD = mean(smd_vals, na.rm = TRUE),
        Max_Abs_SMD = max(smd_vals, na.rm = TRUE)
      )
    })
    do.call(rbind, res)
  })
  smd_summary <- do.call(rbind, smd_list)
  rownames(smd_summary) <- NULL
  return(smd_summary)
}

### 3.5.2 ESS
#### compute_ess_all_datasets()
compute_ess_all_datasets <- function(combined_list, treat, covar) {
  ess_list <- lapply(names(combined_list), function(weight_method) {
    method_list <- combined_list[[weight_method]]
    res <- lapply(names(method_list), function(trim_method) {
      data <- method_list[[trim_method]]
      # Use equal weights if weight column missing
      if (!"weight" %in% colnames(data)) {
        data$weight <- rep(1, nrow(data))
      }
      bal_obj <- cobalt::bal.tab(
        as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
        data = data,
        weights = data$weight,
        un = TRUE,
        s.d.denom = "treated"
      )
      samples <- bal_obj$Observations
      if ("Adjusted" %in% rownames(samples)) {
        control_ess <- samples["Adjusted", "Control"]
        treated_ess <- samples["Adjusted", "Treated"]
      } else {
        control_ess <- samples[1, "Control"]
        treated_ess <- samples[1, "Treated"]
      }
      method_name <- paste(weight_method, trim_method, sep = "_")
      data.frame(
        Method = method_name,
        Control = control_ess,
        Treated = treated_ess
      )
    })
    do.call(rbind, res)
  })
  ess_summary <- do.call(rbind, ess_list)
  rownames(ess_summary) <- NULL
  return(ess_summary)
}

### 3.5.3 Visuals
#### plot_comb_overlap()
plot_comb_overlap <- function(comb_meth_cps = NULL, comb_meth_psid = NULL, treat, covar,
                              prefix_cps = "LDW-CPS1", prefix_psid = "LDW-PSID1") {
  all_combined_list <- list()
  if (!is.null(comb_meth_cps)) all_combined_list$CPS <- comb_meth_cps
  if (!is.null(comb_meth_psid)) all_combined_list$PSID <- comb_meth_psid
  for (ds_name in names(all_combined_list)) {
    combined_list <- all_combined_list[[ds_name]]
    plot_list <- list()
    method_names <- character()
    for (weight_method in names(combined_list)) {
      for (trim_method in names(combined_list[[weight_method]])) {
        full_method_name <- paste(weight_method, trim_method, sep = "_")
        plot_list[[full_method_name]] <- list(
          data = combined_list[[weight_method]][[trim_method]],
          method_name = full_method_name
        )
        method_names <- c(method_names, full_method_name)
      }
    }
    total_plots <- length(plot_list)
    plots_per_page <- 4  # 2x2 layout
    for (i in seq(1, total_plots, by = plots_per_page)) {
      page_plots <- plot_list[i:min(i + plots_per_page - 1, total_plots)]
      par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
      for (p in page_plots) {
        assess_overlap(p$data, treat = treat, cov = covar)
        idx <- match(p$method_name, method_names)
        prefix <- if (ds_name == "CPS" && idx <= 25) prefix_cps else
          if (ds_name == "PSID" && idx <= 25) prefix_psid else ""
        
        title(main = paste0(prefix, " - ", p$method_name))
      }
    }
  }
}

#### plot_comb_love_plots()
plot_comb_love_plots <- function(comb_meth_cps = NULL, comb_meth_psid = NULL, treat, covar,
                                 prefix_cps = "LDW-CPS1", prefix_psid = "LDW-PSID1") {
  all_datasets <- list()
  if (!is.null(comb_meth_cps)) all_datasets$CPS <- comb_meth_cps
  if (!is.null(comb_meth_psid)) all_datasets$PSID <- comb_meth_psid
  for (ds_name in names(all_datasets)) {
    method_list <- all_datasets[[ds_name]]
    method_names <- unlist(lapply(names(method_list), function(weighting) {
      trimmed_list <- method_list[[weighting]]
      sapply(names(trimmed_list), function(trim) paste(weighting, trim, sep = "_"))
    }))
    plot_counter <- 0
    for (weighting in names(method_list)) {
      trimmed_list <- method_list[[weighting]]
      for (trim in names(trimmed_list)) {
        df <- trimmed_list[[trim]]
        plot_counter <- plot_counter + 1
        if (!"weight" %in% names(df)) df$weight <- 1
        bal <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = df,
          weights = df$weight,
          un = TRUE,
          s.d.denom = "treated"
        )
        method_name <- paste(weighting, trim, sep = "_")
        prefix <- ""
        if (ds_name == "CPS" && plot_counter <= 25) prefix <- prefix_cps
        if (ds_name == "PSID" && plot_counter <= 25) prefix <- prefix_psid
        title_text <- paste(prefix, method_name, sep = " - ")
        lp <- cobalt::love.plot(
          bal,
          stats = "mean.diffs",
          absolute = TRUE,
          var.order = "unadjusted",
          thresholds = c(m = .1),
          title = title_text
        )
        print(lp)
      }
    }
  }
}

#### save_comb_hist()
save_comb_hist <- function(comb_meth_cps = NULL, comb_meth_psid = NULL, treat, covar,
                           prefix = "model_a",
                           prefix_cps = "LDW-CPS1", prefix_psid = "LDW-PSID1",
                           path = "../graphs/lalonde") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  all_combined_list <- list()
  if (!is.null(comb_meth_cps)) all_combined_list$CPS <- comb_meth_cps
  if (!is.null(comb_meth_psid)) all_combined_list$PSID <- comb_meth_psid
  
  file_index <- 1
  for (ds_name in names(all_combined_list)) {
    combined_list <- all_combined_list[[ds_name]]
    method_names <- character()
    plot_list <- list()
    for (weight_method in names(combined_list)) {
      for (trim_method in names(combined_list[[weight_method]])) {
        full_method_name <- paste(weight_method, trim_method, sep = "_")
        plot_list[[full_method_name]] <- combined_list[[weight_method]][[trim_method]]
        method_names <- c(method_names, full_method_name)
      }
    }
    total_plots <- length(plot_list)
    for (i in seq_len(total_plots)) {
      data <- plot_list[[i]]
      pdf_file <- file.path(path, sprintf("%s_overlap_%d.pdf", prefix, file_index))
      pdf(pdf_file, width = 8, height = 6)
      assess_overlap(data, treat = treat, cov = covar)
      prefix_str <- if (ds_name == "CPS") prefix_cps else prefix_psid
      title(main = paste0(prefix_str, " - ", method_names[i]))
      dev.off()
      file_index <- file_index + 1
    }
  }
}

#### save_comb_loveplots()
save_comb_loveplots <- function(comb_meth_cps = NULL, comb_meth_psid = NULL, treat, covar,
                                prefix = "model_a",
                                prefix_cps = "LDW-CPS1", prefix_psid = "LDW-PSID1",
                                path = "../graphs/lalonde") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  all_datasets <- list()
  if (!is.null(comb_meth_cps)) all_datasets$CPS <- comb_meth_cps
  if (!is.null(comb_meth_psid)) all_datasets$PSID <- comb_meth_psid
  
  file_index <- 1
  for (ds_name in names(all_datasets)) {
    method_list <- all_datasets[[ds_name]]
    method_names <- unlist(lapply(names(method_list), function(weighting) {
      trimmed_list <- method_list[[weighting]]
      sapply(names(trimmed_list), function(trim) paste(weighting, trim, sep = "_"))
    }))
    for (weighting in names(method_list)) {
      trimmed_list <- method_list[[weighting]]
      for (trim in names(trimmed_list)) {
        df <- trimmed_list[[trim]]
        if (!"weight" %in% names(df)) df$weight <- 1
        bal <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = df,
          weights = df$weight,
          un = TRUE,
          s.d.denom = "treated"
        )
        method_name <- paste(weighting, trim, sep = "_")
        prefix_str <- ifelse(ds_name == "CPS", prefix_cps, prefix_psid)
        pdf_file <- file.path(path, sprintf("%s_balance_%d.pdf", prefix, file_index))
        pdf(pdf_file, width = 8, height = 6)
        lp <- cobalt::love.plot(
          bal,
          stats = "mean.diffs",
          absolute = TRUE,
          var.order = "unadjusted",
          thresholds = c(m = .1),
          title = paste(prefix_str, method_name, sep = " - ")
        )
        print(lp)
        dev.off()
        file_index <- file_index + 1
      }
    }
  }
}

## 3.6 Getting top methods and datasets
#### combine_results()
combine_results <- function(dataset_name) {
  dataset_lower <- tolower(dataset_name)
  # Retrieve matching results
  smd_matching <- get(paste0("smd_matchit.", dataset_lower))
  ess_matching <- get(paste0("ess_matchit.", dataset_lower))
  # Retrieve trimming results
  smd_trimming <- get(paste0("smd_trim.", dataset_lower))
  ess_trimming <- get(paste0("ess_trim.", dataset_lower))
  # Retrieve truncation results
  smd_trunc <- get(paste0("smd_trunc.", dataset_lower))
  ess_trunc <- get(paste0("ess_trunc.", dataset_lower))
  # Retrieve weighting results
  smd_weighting <- get(paste0("smd_weight.", dataset_lower))
  ess_weighting <- get(paste0("ess_weight.", dataset_lower))
  # Retrieve combined SMD and ESS
  # Make sure you have these variables in your environment per dataset
  smd_combined_var <- paste0("smd_all_comb_meth.", dataset_lower)
  ess_combined_var <- paste0("ess_all_comb_meth.", dataset_lower)
  smd_combined <- get(smd_combined_var)
  ess_combined <- get(ess_combined_var)
  # Combine all SMD results
  smd_all <- do.call(rbind, list(
    smd_matching,
    smd_trimming,
    smd_trunc,
    smd_weighting,
    smd_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")]
  ))
    # Combine all ESS results
  ess_all <- do.call(rbind, list(
    ess_matching,
    ess_trimming,
    ess_trunc,
    ess_weighting,
    ess_combined[, c("Method", "Control", "Treated")]
  ))
  # Merge SMD and ESS results by Method
  final_df <- merge(smd_all, ess_all, by = "Method", all = TRUE)
  # Remove dataset suffixes like ".cps_plus" or ".psid_plus" from Method names for cleaner labels
  final_df$Method <- gsub("\\.psid_plus", "", final_df$Method, ignore.case = TRUE)
  final_df$Method <- gsub("\\.cps_plus", "", final_df$Method, ignore.case = TRUE)
    # Reset row names
  rownames(final_df) <- NULL
  return(final_df)
}

#### save_csv()
save_csv <- function(data, filename) {
  folder <- "../tables"
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  file_csv <- file.path(folder, paste0(filename, ".csv"))
  write.csv(data, file = file_csv, row.names = FALSE)
}

#### assess_methods
assess_methods <- function(df) {
  df %>%
    mutate(
      # ess score
      ess_balance_ratio = pmin(Control, Treated) / pmax(Control, Treated),
      ess_total = Control + Treated,
      # Normalize balance_ratio and ess_total
      ess_balance_score = (ess_balance_ratio - min(ess_balance_ratio, na.rm = TRUE)) /
        (max(ess_balance_ratio, na.rm = TRUE) - min(ess_balance_ratio, na.rm = TRUE)),
      ess_size_score = (ess_total - min(ess_total, na.rm = TRUE)) /
        (max(ess_total, na.rm = TRUE) - min(ess_total, na.rm = TRUE)),
      # Combine size and balance equally
      ess_score = 0.5 * ess_balance_score + 0.5 * ess_size_score,
      # smd score
      smd_score = 1 - (Mean_Abs_SMD - min(Mean_Abs_SMD, na.rm = TRUE)) /
        (max(Mean_Abs_SMD, na.rm = TRUE) - min(Mean_Abs_SMD, na.rm = TRUE)),
      # final score
      Score = 0.5 * smd_score + 0.5 * ess_score
    ) %>%
    dplyr::select(Method, Score) %>%
    arrange(desc(Score))
}

#### get_top_methods()
get_top_methods <- function(summary_df, top_n = 5) {
  if (!all(c("Method", "Score") %in% names(summary_df))) {
    stop("Data frame must contain columns 'Method' and 'Score'")
  }
  top_methods_df <- summary_df %>%
    arrange(desc(Score)) %>%
    head(top_n)
  print(top_methods_df)
  return(top_methods_df$Method)
}

#### create_top5_datasets()
create_top5_datasets <- function(dataset_list, top5_method_names) {
  lapply(top5_method_names, function(method_name) {
    if (!method_name %in% names(dataset_list)) {
      stop(paste0("Method '", method_name, "' not found in the dataset lookup list"))
    }
    ds <- dataset_list[[method_name]]
    if (inherits(ds, "matchit")) {
      # for matchit class, extract matched data
      return(as.data.frame(match.data(ds)))
    } else if (is.data.frame(ds)) {
      # if already dataframe, return as is
      return(ds)
    } else if (is.vector(ds)) {
      # merge with original data
      if (!"original" %in% names(dataset_list)) {
        stop("Original dataset not found in dataset_list to merge weights")
      }
      original_df <- dataset_list[["original"]]
      if (length(ds) != nrow(original_df)) {
        stop(paste0("Length of weight vector for method '", method_name, "' does not match number of rows in original dataset"))
      }
      # create new data frame with weights appended
      weighted_df <- original_df
      weighted_df[["weights"]] <- ds
      return(weighted_df)
    } else {
      stop(paste("Unsupported data type for method", method_name))
    }
  })
}

#### save_top5_individual_files()
save_top5_individual_files <- function(combined_methods_list, top5_method_names, prefix) {
  for (i in seq_along(top5_method_names)) {
    method_name <- top5_method_names[i]
    if (!method_name %in% names(combined_methods_list)) {
      warning(paste0("Method '", method_name, "' not found in combined methods list"))
      next
    }
    dataset_to_save <- combined_methods_list[[method_name]]
    file_name <- sprintf("data/top%d_%s_method_%s.RData", i, prefix, method_name)
    save(dataset_to_save, file = file_name)
    cat("Saved:", file_name, "\n")
  }
}

# 4. Estimating
## 4.1 ATT
#### estimate_all()
# difference in means
diff <- function(data, Y, treat) {
  fml <- as.formula(paste(Y, "~", treat))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  return(out) # extract coef, se, ci.lower, ci.upper
}

# regression adjustment
reg <- function(data, Y, treat, covar) {
  fml <- as.formula(paste(Y, "~", treat, "+", paste(covar, collapse = " + ")))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  # extract coef, se, ci.lower, ci.upper
  return(out)
}

# matching
#library(Matching)
matching <- function(data, Y, treat, covar) {
  m.out <- Match(Y = data[, Y], Tr = data[, treat], X = data[, covar], Z = data[, covar],
                 estimand = "ATT", M = 5, replace = TRUE, ties = TRUE, BiasAdjust = TRUE)
  out <- c(m.out$est[1], m.out$se[1], m.out$est[1] - 1.96 * m.out$se[1],
           m.out$est[1] + 1.96 * m.out$se[1])
  return(out)
}

# psm
psm <- function(data, Y, treat, covar) {
  ps <- probability_forest(X = data[, covar],
                           Y = as.factor(data[,treat]), seed = 1234, num.trees = 4000)$predictions[,2]
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

# OM (reg)
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

# OM (grf)
#library(grf)
om.grf <- function(data, Y, treat, covar) {
  tr <- which(data[, treat] == 1)
  co <- which(data[, treat] == 0)
  out.co <- regression_forest(X = data[co, covar, drop = FALSE], Y = as.vector(data[co, Y]) )
  Y.tr.hat <- as.vector(unlist(predict(out.co, newdata = data[tr, covar, drop = FALSE])))
  newdata <- cbind.data.frame(Y = c(data[tr, Y], Y.tr.hat), treat = rep(c(1, 0), each = length(tr)))
  out <- summary(lm_robust(Y ~ treat, data = newdata, se_type = "stata"))$coefficients["treat", c(1, 2, 5, 6)]
  return(out)
}

# IPW
ipw <- function(data, Y, treat, covar) {
  ps <- probability_forest(X = data[, covar, drop = FALSE], Y = as.factor(data[, treat]), seed = 1234)$predictions[,2]
  fml <- as.formula(paste(Y, "~", treat))
  weights <- rep(1, nrow(data))
  co <- which(data[, treat] == 0)
  weights[co] <- ps[co]/(1-ps[co])
  out <- summary(lm_robust(fml, data = data, weights = weights, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  # extract coef, se, ci.lower, ci.upper
  return(out)
}

# CBPS
#library("CBPS")
cbps <- function(data, Y, treat, covar) {
  fml <- as.formula(paste(treat, "~", paste(covar, collapse = " + ")))
  ps <- quiet(CBPS(fml, data = data, standardize = TRUE)$fitted.values)
  fml <- as.formula(paste(Y, "~", treat))
  weights <- rep(1, nrow(data))
  co <- which(data[, treat] == 0)
  weights[co] <- ps[co]/(1-ps[co])
  out <- summary(lm_robust(fml, data = data, weights = weights, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  return(out)
}

# ebal
#library(hbal)
ebal <- function(data, Y, treat, covar) {
  ebal.out <- hbal::hbal(Y = Y, Treat = treat, X = covar,  data = data, expand.degree = 1)
  out <- hbal::att(ebal.out, dr = FALSE)[1, c(1, 2, 5, 6)]
  return(out)
}

# hbal
# hbal <- function(data, Y, treat, covar) {
#   hbal.out <- hbal::hbal(Y = Y, Treat = treat, X = covar,  data = data, expand.degree = 2, # cv = TRUE)
#   out <- hbal::att(hbal.out, dr = FALSE)[1, c(1, 2, 5, 6)]
#   return(out)
# }

# AIPW
aipw <- function(data, Y, treat, covar) {
  #library("grf")
  for (var in c(Y, treat, covar)) {
    data[, var] <- as.vector(data[, var])
  }
  c.forest <- causal_forest(X = data[, covar, drop = FALSE], Y = data[, Y],
                            W = data[, treat], seed = 1234)
  att <- average_treatment_effect(c.forest, target.sample = "treated", method = "AIPW")
  att <- c(att, att[1] - 1.96 * att[2], att[1] + 1.96 * att[2])
  return(att)
}

aipw.match <- function(data, Y, treat, covar) { # match on ps
  ps <- probability_forest(X = data[, covar], Y = as.factor(data[, treat]), seed = 1234)$predictions[,2]
  m.out <- Match(Y = data[, Y], Tr = data[, treat], X = ps,
                 estimand = "ATT", M = 1, replace = FALSE, ties = FALSE, BiasAdjust = FALSE)
  mb <- quiet(MatchBalance(treat ~ ps, data = data, match.out = m.out, nboots= 0))
  ks <- mb$AfterMatching[[1]]$ks$ks$statistic
  s <- data[c(m.out$index.treated, m.out$index.control), ]
  out <- aipw(s, Y, treat, covar)
  #return(out)
  return(c(out, ks))
}

aipw_ow <- function(data, Y, treat, covar) {
  for (var in c(Y, treat, covar)) {
    data[, var] <- as.vector(data[, var])
  }
  X <- data[, covar, drop = FALSE]
  Y <- data[, Y]
  W <- data[, treat]
  # run dml_with_smoother with AIPW_ATT
  dml_fit <- dml_with_smoother(Y = Y, D = W, X = X,
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
}

### This script checks for robustness by estimating original model
### using double/debiased machine learning using DoubleML package
dml <-function(data, Y = NULL, treat = NULL, covar = NULL, clust_var = NULL, ml_l = lrn("regr.lm"), ml_m = lrn("regr.lm")){
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
    dml_dat = DoubleMLData$new(dat,
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
    dml_dat = DoubleMLClusterData$new(dat,
                                      y_col = Y,
                                      d_cols = treat,
                                      cluster_cols = clust_var,
                                      use_other_treat_as_covariate = FALSE,
                                      x_cols = covar)
  }
  # Set active treatment treatment
  dml_dat$set_data_model(treat)
  # Estimate with DML
  set.seed(pi)
  dml_mod = DoubleMLPLR$new(dml_dat, ml_l=ml_l, ml_m=ml_m)
  quiet(dml_mod$fit())
  out = c(dml_mod$coef[treat], dml_mod$se[treat], dml_mod$confint()[treat,])
  return(out)
  
}

# execute all estimators
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

#### plot_coef()
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
                      textsize = 1
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

#### plot_att_panels()
plot_att_panels <- function(all_outs, plot_titles, band, est, ylim = c(-15500, 5500), plots_per_page = 4, ylab = "Estimate", textsize = 1) {
  num_pages <- ceiling(length(all_outs) / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx   <- min(page * plots_per_page, length(all_outs))
    par(mfrow = c(plots_per_page, 1), mar = c(3, 4, 3, 2))
    for (i in start_idx:end_idx) {
      out <- all_outs[[i]]
      plot_coef(out, 
                band = band, 
                line = est,
                ylim = ylim,
                main = plot_titles[i],
                ylab = ylab,
                textsize = textsize)
    }
  }
}

#### save_att_panels()
save_att_panels <- function(
    all_outs, plot_titles, band, est, prefix,
    plots_per_page = 4, ylab = "Estimate", textsize = 1) {
  folder <- "../graphs/lalonde"
  ylim   <- c(-15500, 5500)
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  num_pages <- ceiling(length(all_outs) / plots_per_page)
  for (page in seq_len(num_pages)) {
    file_name <- file.path(folder, paste0(prefix, "_", page, ".pdf"))
    pdf(file_name, width = 8, height = 11)
    par(mfrow = c(plots_per_page, 1), mar = c(3,4,3,2))
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(all_outs))
    for (i in start_idx:end_idx) {
      plot_coef(all_outs[[i]],
                band = band, line = est, ylim = ylim,
                main = plot_titles[i], ylab = ylab, textsize = textsize)
    }
    dev.off()
  }
}

#### create_matrix_results()
create_matrix_results <- function(all_outs, sample_names) {
  n_samples <- length(sample_names)
  n_estimators <- nrow(all_outs[[1]])
  result_mat <- matrix("", nrow = n_estimators + 1, ncol = n_samples * 2)
  # Set up alternating column names
  cnames <- character(n_samples * 2)
  for (j in seq_along(sample_names)) {
    cnames[(j-1)*2 + 1] <- sample_names[j]
    cnames[(j-1)*2 + 2] <- "" # SE/parenthesis column
  }
  colnames(result_mat) <- cnames
  estimator_names <- rownames(all_outs[[1]])
  rownames(result_mat) <- c("Experimental Benchmark", estimator_names)
  # Fill values
  for (j in seq_along(all_outs)) {
    out <- all_outs[[j]]
    result_mat[1, (j-1)*2 + 1] <- sprintf("%.2f", out[1, 1])
    result_mat[1, (j-1)*2 + 2] <- paste0("(", sprintf("%.2f", out[1, 2]), ")")
    for (i in 2:(n_estimators+1)) {
      result_mat[i, (j-1)*2 + 1] <- sprintf("%.2f", out[i-1, 1])
      result_mat[i, (j-1)*2 + 2] <- paste0("(", sprintf("%.2f", out[i-1, 2]), ")")
    }
  }
  return(result_mat)
}

####  eval_att()
eval_att <- function(result) {
  data.frame(
    Mean_SE = mean(result[, "SE"], na.rm = TRUE), # mean standard error across all estimator results
    Min_Estimate = min(result[, "Estimate"], na.rm = TRUE), # maximum estimate of all estimator results
    Max_Estimate = max(result[, "Estimate"], na.rm = TRUE), # minimum estimate of all estimator results
    Diff_Estimate = max(result[, "Estimate"], na.rm = TRUE) - min(result[, "Estimate"], na.rm = TRUE)
  )
}

## 4.2 CATT
#### plot_catt_panels()
plot_catt_panels <- function(all_catt, plot_titles, plots_per_page = 4, range = c(-8000, 8000)) {
  num_pages <- ceiling((length(all_catt) - 1) / plots_per_page)
  # Experimental reference (first panel)
  catt_ldw <- all_catt[[1]]$catt
  att_ldw  <- all_catt[[1]]$att[1]
  id_ldw   <- if (!is.null(all_catt[[1]]$id)) all_catt[[1]]$id else seq_along(catt_ldw)
  
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 2 # skip experimental vs itself
    end_idx   <- min(page * plots_per_page + 1, length(all_catt))
    par(mfrow = c(plots_per_page, 1), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
      other <- all_catt[[i]]
      main_label <- plot_titles[i]
      catt2 <- other$catt
      att2  <- other$att[1]
      id2   <- if (!is.null(other$id)) other$id else seq_along(catt2)
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
        ylab  = main_label,
        main  = main_label,
        axes.range = range
      )
    }
  }
}

#### save_catt_panels()
save_catt_panels <- function(all_catt, plot_titles, range = c(-8000, 8000), prefix = "model_a", plots_per_page = 4) {
  dir.create("../graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  catt_ldw <- all_catt[[1]]$catt
  att_ldw  <- all_catt[[1]]$att[1]
  id_ldw   <- if (!is.null(all_catt[[1]]$id)) all_catt[[1]]$id else seq_along(catt_ldw)
  num_panels <- length(all_catt) - 1   # skip experimental on page
  num_pages  <- ceiling(num_panels / plots_per_page)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 2  # always skip all_catt[[1]]
    end_idx   <- min(page * plots_per_page + 1, length(all_catt))
    plots_this_page <- end_idx - start_idx + 1
    file_name <- sprintf("../graphs/lalonde/%s_catt_estimates_%d.pdf", prefix, page)
    pdf(file = file_name, width = 10, height = 12)
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
    for (i in start_idx:end_idx) {
      other <- all_catt[[i]]
      main_label <- plot_titles[i]
      catt2 <- other$catt
      att2  <- other$att[1]
      id2   <- if (!is.null(other$id)) other$id else seq_along(catt2)
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
        ylab  = main_label,
        main  = main_label,
        axes.range = range
      )
    }
    # blanks if fewer than 4 plots on last page
    if (plots_this_page < plots_per_page) {
      for (k in seq_len(plots_per_page - plots_this_page)) plot.new()
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

## 4.3 QTET
#### plot_qte_top()
plot_qte_top <- function(qtet_top, qtet_top0, bm, plot_titles, main_start = 1, ylim = NULL, col = NULL) {
  n <- length(qtet_top)
  for (i in 1:n) {
    main_title <- plot_titles[main_start + i - 1]
    mod <- qtet_top[[i]]
    mod2 <- qtet_top0[[i]]
    plot_qte(mod, mod2, bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
  }
}

#### save_qtet()
save_qtet <- function(plots, plot_titles = NULL, main_start = 1, ylim = NULL, col = NULL, prefix = "model_a") {
  dir.create("../graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  n <- length(plots)
  for (i in seq_len(n)) {
    p <- plots[[i]]
    main_title <- if (is.null(plot_titles)) p$main else plot_titles[main_start + i - 1]
    clean_title <- gsub("[^a-zA-Z0-9]", "_", main_title)
    file_name <- sprintf("../graphs/lalonde/%s_qtet_estimates_%s.pdf", prefix, clean_title)
    pdf(file = file_name, width = 7, height = 5)
    plot_qte(p$mod, p$mod0, p$bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
    dev.off()
  }
}

#### save_qte_top()
save_qte_top <- function(qtet_top, qtet_top0, bm, plot_titles, main_start = 1,
                                ylim = NULL, col = NULL, prefix = "model_a_top") {
  n <- length(qtet_top)
  dir.create("../graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(n)) {
    mod <- qtet_top[[i]]
    mod2 <- qtet_top0[[i]]
    main_title <- plot_titles[main_start + i - 1]
    clean_title <- gsub("[^a-zA-Z0-9]", "_", main_title)
    file_name <- sprintf("../graphs/lalonde/%s_qte_estimates_%s.pdf", prefix, clean_title)
    pdf(file = file_name, width = 7, height = 5)
    plot_qte(mod, mod2, bm, main = main_title, ylim = ylim, col = col)
    legend("bottomleft", legend = c("Experimental", "Unadjusted", "Adjusted"),
           lty = 1, pch = c(16, 17, 16), col = c(4, 2, 1), bty = "n")
    dev.off()
  }
}

## 4.4 Assessing outcome weights (OW)
#### get_res_att()
get_res_att <- function(dataset_list, Y, treat, covar,
                            estimators = "AIPW_ATT",
                            smoother = "honest_forest",
                            n_cf_folds = 5,
                            n_reps = 1) {
  lapply(dataset_list,
      function(data) {
        dml_with_smoother(
          Y = data[[Y]],
          D = data[[treat]],
          X = data[, covar, drop = FALSE],
          estimators = estimators,
          smoother = smoother,
          n_cf_folds = n_cf_folds,
          n_reps = n_reps
      )
    }
  )
}

#### derive_ow()
derive_ow <- function(results_list) {
lapply(results_list, function(res) get_outcome_weights(res))
}

#### plot_ow()
plot_ow <- function(outcome_weights, plot_titles = NULL, breaks = 50, 
                    col = "#ff000080", xlab = "Outcome Weight", 
                    estimand = "AIPW-ATT") {
  N <- length(outcome_weights)
  for (i in seq_len(N)) {
    weights <- outcome_weights[[i]]$omega[estimand, ]
    main_title <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    hist(weights, breaks = breaks, main = main_title, xlab = xlab, col = col)
    mtext(paste("N =", length(weights)), side = 3, line = -1.5, cex = 0.8)
  }
  par(mfrow = c(1, 1))
}

#### eval_ow_prop()
eval_ow_prop <- function(outcome_weights, dataset_list, plot_titles = NULL, treat_var = "treat", estimand = "AIPW-ATT") {
  results <- lapply(seq_along(outcome_weights), function(i) {
    ow <- outcome_weights[[i]]$omega[estimand, ]
    treat <- dataset_list[[i]][[treat_var]]
    method <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    sum_treated <- sum(ow[treat == 1])
    sum_untreated <- sum(ow[treat == 0])
    data.frame(
      Method = method,
      Sum_Treated = sum_treated,
      Sum_Untreated = sum_untreated,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#### save_ow()
save_ow <- function(outcome_weights, plot_titles = NULL,
                    breaks = 50, col = "#ff000080", xlab = "Outcome Weight",
                    prefix = "model_a", estimand = "AIPW-ATT") {
  dir.create("../graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  N <- length(outcome_weights)
  for (i in seq_len(N)) {
    file_name <- sprintf("../graphs/lalonde/%s_outcomewt_%d.pdf", prefix, i)
    pdf(file = file_name, width = 8, height = 6)
    weights <- outcome_weights[[i]]$omega[estimand, ]
    main_title <- if (!is.null(plot_titles)) plot_titles[i] else paste("Dataset", i)
    hist(weights, breaks = breaks, main = main_title, xlab = xlab, col = col)
    mtext(paste("N =", length(weights)), side = 3, line = -1.5, cex = 0.8)
    dev.off()
  }
}

# 5. Sensitivity Analysis
#### check_filter_data()
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
  folder <- "../graphs/lalonde"
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  for (i in seq_along(filtered_datasets)) {
    idx <- i
    file_name <- file.path(folder, paste0(prefix, "_sensitivity_", idx, ".pdf"))
    pdf(file_name, width = 8, height = 8)
    sens_ana(filtered_datasets[[i]], Y, treat, covar, bm, kd = 1:3)
    if (!missing(plot_titles) && length(plot_titles) >= idx) {
      title(main = plot_titles[idx])
    }
    dev.off()
 }
}