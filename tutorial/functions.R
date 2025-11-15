# 1 Set up
# 1.1 Installation
packages <- c("CBPS", "cobalt", "data.table", "dplyr", "DT", "DoubleML", "ebal", "estimatr", "ggplot2", "gridExtra", "grf", "hbal", "highr", "highs", 
              "kableExtra", "MatchIt", "Matching", "mlr3", "mlr3learners", "OutcomeWeights", "optmatch", "optweight", "quickmatch", 
              "readr", "rgenoud", "tidyr", "tidyverse", "ppcor", "WeightIt"
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

# 2. Improving covariate balance and overlap
## 2.3 Trimming
### 2.3.1 Propensity score threshold trimming (similar to tutorial of Imbens & Xu (2024))
# ps_trim() *** 
ps_trim <- function(data, ps = "ps_assoverlap", threshold = 0.9) { 
  sub <- data[which(data[, ps] < threshold), ]
  return(sub)
}

### 2.3.2 Common range trimming
#### common_range_trim ()
common_range_trim <- function(data, ps = "ps_assoverlap", treat = "treat") {
  lower_cut <- min(data[[ps]][data[[treat]] == 1], na.rm = TRUE)
  upper_cut <- max(data[[ps]][data[[treat]] == 0], na.rm = TRUE)
  sub <- data[data[[ps]] >= lower_cut & data[[ps]] <= upper_cut, ]
  return(sub)
}

### 2.3.3 Crump trimming
#### crump_trim ()
crump_trim <- function(data, ps = "ps_assoverlap", lower = 0.1, upper = 0.9) {
  sub <- data[data[[ps]] >= lower & data[[ps]] <= upper, ]
  return(sub)
}

### 2.3.4 Stuermer trimming
#### stuermer_trim ()
stuermer_trim <- function(data, treat = "treat", ps = "ps_assoverlap", 
                          lower_percentile = 0.05, upper_percentile = 0.95) {
  treated_ps   <- data[[ps]][data[[treat]] == 1]
  untreated_ps <- data[[ps]][data[[treat]] == 0]
  lower_cutoff <- quantile(treated_ps, probs = lower_percentile, na.rm = TRUE)
  upper_cutoff <- quantile(untreated_ps, probs = upper_percentile, na.rm = TRUE)
  sub <- data[data[[ps]] >= lower_cutoff & data[[ps]] <= upper_cutoff, ]
  return(sub)
}

### 2.3.5 Walker trimming
#### walker_trim ()
walker_trim <- function(data, treat = "treat", ps = "ps_assoverlap", 
                        lower_cutoff = 0.3, upper_cutoff = 0.7) {
  treat_prevalence  <- mean(data[[treat]], na.rm = TRUE)
  logit_ps          <- log(data[[ps]] / (1 - data[[ps]]))
  logit_prevalence  <- log(treat_prevalence / (1 - treat_prevalence))
  preference_score  <- 1 / (1 + exp(-(logit_ps - logit_prevalence)))
  sub <- data[preference_score >= lower_cutoff & preference_score <= upper_cutoff, ]
  return(sub)
}

### 2.3.6 Reestimate propensity scores
#### ps_estimate()***
ps_estimate <- function(data, treat, cov, odds = TRUE, num.trees = NULL, seed = 1234, breaks = 50, xlim = NULL, ylim = NULL) {
  if(is.null(num.trees))
  {p.forest1 <- probability_forest(X = data[, cov],
                                    Y = as.factor(data[,treat]), seed = seed)}
  else
  {p.forest1 <- probability_forest(X = data[, cov],
                                    Y = as.factor(data[,treat]), seed = seed, num.trees = num.trees)}
  data$ps_assoverlap <- p.forest1$predictions[,2]
  data$ps_assoverlap[which(abs(data$ps_assoverlap) <= 1e-7)] <- 1e-7
  return(data)
}

## 2.5 Integrated methods
### 2.5.1 Trimming and matching
#### attach_matchit()
attach_matchit <- function(model, data_list, ..., verbose = FALSE) {
  matchit_results <- vector("list", length(data_list))
  names(matchit_results) <- names(data_list)
  for (i in seq_along(data_list)) {
    res <- tryCatch(
      matchit(model, data = data_list[[i]], ...),
      error = function(e) {
        return(NULL)
      }
    )
    matchit_results[[i]] <- res
  }
  return(matchit_results)
}

# 3. Reassessing methods
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
      s.d.denom = "treated"
    )
    smds <- bal$Balance$Diff.Adj
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

## 3.3 Trimming and matching
### 3.3.1 SMD
#### compute_abs_smd_matchit()
compute_abs_smd_matchit <- function(match_list, data_list) {
  match_sizes <- vapply(match_list, length, integer(1))
  smd_results <- list()
  for (method_name in names(match_list)) {
    method_sublist <- match_list[[method_name]]
    n_match <- length(method_sublist)
    for (i in seq_len(n_match)) {
      match_obj <- method_sublist[[i]]
      data <- data_list[[i]] 
      if (is.null(match_obj)) {
        warning(sprintf(
          "compute_abs_smd_matchit(): skipping '%s' match object %d because it is NULL; check attach_matchit() failures. match_list lengths: %s",
          method_name,
          i,
          paste(sprintf("%s=%d", names(match_sizes), match_sizes), collapse = ", ")
        ), call. = FALSE)
        next
      }
      if (is.null(data)) {
        warning(sprintf(
          "compute_abs_smd_matchit(): skipping '%s' source data %d because it is NULL; verify trimming outputs. match_list lengths: %s",
          method_name,
          i,
          paste(sprintf("%s=%d", names(match_sizes), match_sizes), collapse = ", ")
        ), call. = FALSE)
        next
      }
      bal <- bal.tab(
        match_obj,
        data = data,
        stats = "mean.diffs",
        un = TRUE,
        s.d.denom = "treated"
      )
      smds <- bal$Balance$Diff.Adj
      smd_vals <- abs(smds)
      mean_smd <- mean(smd_vals, na.rm = TRUE)
      max_smd <- max(smd_vals, na.rm = TRUE)
      smd_results[[paste(method_name, i, sep = "_")]] <- data.frame(
        Method = method_name,
        MatchIndex = i,
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

### 3.3.2 OVL
#### compute_ovl_matchit()
compute_ovl_matchit <- function(match_list, data_list, ps = "ps_assoverlap", treat = "treat", 
                                covar = NULL, num.trees = NULL, seed = 42, n_points = 512) {
  if (!is.list(match_list) || !length(match_list)) {
    stop("compute_ovl_matchit(): 'match_list' must be a non-empty list of MatchIt outputs")
  }
  if (!is.list(data_list) || !length(data_list)) {
    stop("compute_ovl_matchit(): 'data_list' must be a non-empty list of trimming datasets")
  }

  match_sizes <- vapply(match_list, length, integer(1))
  data_count <- length(data_list)
  if (!all(match_sizes == data_count)) {
    warning(sprintf(
      "compute_ovl_matchit(): not all match method lists align with data_list length (%d). sizes: %s",
      data_count,
      paste(sprintf("%s=%d", names(match_sizes), match_sizes), collapse = ", ")
    ), call. = FALSE)
  }

  ovl_results <- list()

  for (method_name in names(match_list)) {
    method_sublist <- match_list[[method_name]]
    if (!is.list(method_sublist)) {
      warning(sprintf("compute_ovl_matchit(): skipping method '%s' because entry is not a list", method_name), call. = FALSE)
      next
    }
    for (i in seq_along(method_sublist)) {
      match_obj <- method_sublist[[i]]

      if (is.null(match_obj)) {
        warning(sprintf(
          "compute_ovl_matchit(): skipping '%s' match object %d because it is NULL; check attach_matchit() failures",
          method_name,
          i
        ), call. = FALSE)
        next
      }

      if (!inherits(match_obj, "matchit")) {
        warning(sprintf(
          "compute_ovl_matchit(): skipping '%s' match object %d because it is class '%s' not 'matchit'",
          method_name,
          i,
          paste(class(match_obj), collapse = "/")
        ), call. = FALSE)
        next
      }

      if (i > length(data_list) || is.null(data_list[[i]])) {
        warning(sprintf(
          "compute_ovl_matchit(): trimming dataset %d missing for method '%s'; verify data_list ordering",
          i,
          method_name
        ), call. = FALSE)
        next
      }

      data <- data_list[[i]]
      if (!is.data.frame(data)) {
        warning(sprintf(
          "compute_ovl_matchit(): data_list[[%d]] for method '%s' is type '%s' not data.frame",
          i,
          method_name,
          paste(class(data), collapse = "/")
        ), call. = FALSE)
        next
      }
      if (!nrow(data)) {
        warning(sprintf(
          "compute_ovl_matchit(): data_list[[%d]] for method '%s' has zero rows",
          i,
          method_name
        ), call. = FALSE)
        next
      }

      matched_data <- match.data(match_obj, data = data)

      if (!is.data.frame(matched_data) || !nrow(matched_data)) {
        warning(sprintf(
          "compute_ovl_matchit(): matched data empty for method '%s' (index %d)",
          method_name,
          i
        ), call. = FALSE)
        next
      }

      if (ps %in% names(matched_data)) {
        ps_vals <- as.numeric(matched_data[[ps]])
      } else {
        if (is.null(covar)) {
          stop(paste0("Propensity scores not found for '", method_name, "' on iteration ", i,
                      " and 'covar' not provided to estimate them"))
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
          "compute_ovl_matchit(): insufficient treated/control overlap after cleaning for method '%s' (index %d)",
          method_name,
          i
        ), call. = FALSE)
        next
      }

      treated_scores <- ps_vals[treat_vals == 1]
      control_scores <- ps_vals[treat_vals == 0]

      dens_treated <- density(treated_scores, from = 0, to = 1, n = n_points)
      dens_control <- density(control_scores, from = 0, to = 1, n = n_points)
      x_grid <- dens_treated$x
      min_density <- pmin(dens_treated$y, dens_control$y)
      bin_width <- x_grid[2] - x_grid[1]
      ovl <- sum(min_density) * bin_width

      ovl_results[[paste(method_name, i, sep = "_")]] <- data.frame(Method = paste(method_name, i, sep = "_"), OVL = ovl)
    }
  }

  if (!length(ovl_results)) {
    return(data.frame(Method = character(0), OVL = numeric(0)))
  }

  final_df <- do.call(rbind, ovl_results)
  rownames(final_df) <- NULL
  return(final_df)
}


## 3.6 Getting top methods and datasets
#### combine_results()
combine_results <- function(dataset_name) {
  dataset_lower <- tolower(dataset_name)
  # retrieve individual method results
  smd_trimming <- get(paste0("smd_trim.", dataset_lower))
  ovl_trimming <- get(paste0("ovl_trim.", dataset_lower))
  smd_trunc <- get(paste0("smd_trunc.", dataset_lower))
  ovl_trunc <- get(paste0("ovl_trunc.", dataset_lower))
  smd_trim_match_combined <- get(paste0("smd_trim_match_comb.", dataset_lower))
  ovl_trim_match_combined <- get(paste0("ovl_trim_match_comb.", dataset_lower))
  smd_trim_weight_combined <- get(paste0("smd_trim_weight_comb.", dataset_lower))
  ovl_trim_weight_combined <- get(paste0("ovl_trim_weight_comb.", dataset_lower))
  smd_all <- do.call(rbind, list(
    smd_trimming,
    smd_trim_match_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")],
    smd_trim_weight_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")],
    smd_trunc_match_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")],
    smd_trunc_weight_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")]
  ))
  # combine all OVL results
  ovl_all <- do.call(rbind, list(
    ovl_trimming,
    ovl_trim_match_combined[, c("Method", "OVL")],
    ovl_trim_weight_combined[, c("Method", "OVL")],
    ovl_trunc_match_combined[, c("Method", "OVL")],
    ovl_trunc_weight_combined[, c("Method", "OVL")]
  ))
  # merge absolute SMD and OVL results by method
  final_df <- merge(smd_all, ovl_all, by = "Method", all = TRUE)
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
    dplyr::select(Method, OVL) %>%
    arrange(desc(OVL))
}

#### rank_methods()


#### get_top_methods()
get_top_methods <- function(summary_df, top_n = 5) {
  if (!all(c("Method", "Score") %in% names(summary_df))) {
    stop("Data frame must contain columns 'Method' and 'Score'")
  }
  top_methods_df <- summary_df %>%
    arrange(desc(Score)) %>%
    head(top_n)
  return(top_methods_df$Method)
}

#### create_top5_datasets()
create_top5_datasets <- function(dataset_list, top5_method_names, weight_list, trunc_list, trim_list) {
  lapply(top5_method_names, function(method_name) {
      if (!method_name %in% names(dataset_list)) {
        stop(paste0("Method '", method_name, "' not found in the dataset lookup list"))
      }
      ds <- dataset_list[[method_name]]
      # matching 
      if (inherits(ds, "matchit")) {
        return(as.data.frame(match.data(ds)))
      # dataframe (trimming and truncation)
      } else if (is.data.frame(ds)) {
        return(ds)
      } else {
        stop(paste("Unsupported data type for method", method_name))
      }
    })
}  

#### save_top5_datasets()
save_top5_datasets <- function(combined_methods_list, top5_method_names, prefix) {
  for (i in seq_along(top5_method_names)) {
    method_name <- top5_method_names[i]
    if (!method_name %in% names(combined_methods_list)) {
      warning(paste0("Method '", method_name, "' not found in combined methods list"))
      next
    }
    dataset_to_save <- combined_methods_list[[method_name]]
    file_name <- sprintf("tutorial/data/top%d_%s_method_%s.RData", i, prefix, method_name)
    save(dataset_to_save, file = file_name)
  }
}

# 4. Estimating
## 4.1 ATT
# difference in means
# diff()
diff <- function(data, Y, treat) {
  fml <- as.formula(paste(Y, "~", treat))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  return(out) # extract coef, se, ci.lower, ci.upper
}

# regression adjustment
# reg()
reg <- function(data, Y, treat, covar) {
  fml <- as.formula(paste(Y, "~", treat, "+", paste(covar, collapse = " + ")))
  out <- summary(lm_robust(fml, data = data, se_type = "stata"))$coefficients[treat, c(1, 2, 5, 6)]
  # extract coef, se, ci.lower, ci.upper
  return(out)
}

# matching
# library(Matching)
# matching()
matching <- function(data, Y, treat, covar) {
  m.out <- Match(Y = data[, Y], Tr = data[, treat], X = data[, covar], Z = data[, covar],
                 estimand = "ATT", M = 5, replace = TRUE, ties = TRUE, BiasAdjust = TRUE)
  out <- c(m.out$est[1], m.out$se[1], m.out$est[1] - 1.96 * m.out$se[1],
           m.out$est[1] + 1.96 * m.out$se[1])
  return(out)
}

# psm()
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

# OM (reg)
# om.reg()
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
# library(grf)
# om.grf()
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
# ipw()
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

# CBPS
# library("CBPS")
# cbps()
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

# ebal()
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

# aipw_grf()
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

# aipw.match()
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

# aipw_ow() ***
aipw_ow <- function(data, Y, treat, covar) {
  tryCatch({
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
  }, error = function(e) {
    cat("Error in aipw_ow method:", e$message, "\n")
    return(c(NA, NA, NA, NA))
  })
}

### this script checks for robustness by estimating original model
### using double/debiased machine learning using DoubleML package
dml <-function(data, Y = NULL, treat = NULL, covar = NULL, clust_var = NULL, ml_l = lrn("regr.lm"), ml_m = lrn("regr.lm")){
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

# estimate_all ***
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

#### plot_coef() ***
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
save_att_panels <- function(out_list, plot_titles, band_list, est_list,
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
create_matrix_results <- function(all_outs, sample_names) {
  n_samples <- length(sample_names)
  n_estimators <- nrow(all_outs[[1]])
  # initialize empty matrix
  result_mat <- matrix("", nrow = n_estimators + 1, ncol = n_samples * 2)
  # alternating column names (estimate / SE)
  cnames <- character(n_samples * 2)
  for (j in seq_along(sample_names)) {
    cnames[(j - 1) * 2 + 1] <- sample_names[j]
    cnames[(j - 1) * 2 + 2] <- ""
  }
  colnames(result_mat) <- cnames
  estimator_names <- rownames(all_outs[[1]])
  rownames(result_mat) <- c("Experimental Benchmark", estimator_names)
  # extract benchmark estimates 
  bench.exp  <- all_outs[[1]][1, 1:2]   # experimental LDW
  bench.cps  <- all_outs[[2]][1, 1:2]   # LDW-CPS1 trimmed benchmark
  bench.psid <- all_outs[[3]][1, 1:2]   # LDW-PSID1 trimmed benchmark
  # fill results column by column
  for (j in seq_along(all_outs)) {
    out <- all_outs[[j]]
    # assign benchmark row depending on dataset type
    if (j <= 3) {  
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.exp[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.exp[2]), ")")
    } 
    else if (grepl("CPS", sample_names[j], ignore.case = TRUE)) {
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.cps[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.cps[2]), ")")
    } 
    else if (grepl("PSID", sample_names[j], ignore.case = TRUE)) {
      result_mat[1, (j - 1) * 2 + 1] <- sprintf("%.2f", bench.psid[1])
      result_mat[1, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", bench.psid[2]), ")")
    }
    # fill in estimates + SEs
    for (i in 2:(n_estimators + 1)) {
      result_mat[i, (j - 1) * 2 + 1] <- sprintf("%.2f", out[i - 1, 1])
      result_mat[i, (j - 1) * 2 + 2] <- paste0("(", sprintf("%.2f", out[i - 1, 2]), ")")
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

## 4.2 CATT
#### plot_catt_panels()
plot_catt_panels <- function(exp_catt, catt_list, plot_titles, plots_per_page = 4, range = c(-8000, 8000)) {
  n <- length(catt_list)
  num_pages <- ceiling(n / plots_per_page)
  catt_ldw <- exp_catt$catt
  att_ldw  <- exp_catt$att[1]
  id_ldw   <- if (!is.null(exp_catt$id)) exp_catt$id else seq_along(catt_ldw)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx   <- min(page * plots_per_page, n)
    par(mfrow = c(2, 2), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
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

#### save_catt_panels()
save_catt_panels <- function(
    exp_catt, 
    catt_list, 
    plot_titles, 
    prefix = "catt_top5", 
    plots_per_page = 4, 
    range = c(-8000, 8000), 
    folder = "graphs/lalonde"
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  n <- length(catt_list)
  num_pages <- ceiling(n / plots_per_page)
  catt_ldw <- exp_catt$catt
  att_ldw  <- exp_catt$att[1]
  id_ldw   <- if (!is.null(exp_catt$id)) exp_catt$id else seq_along(catt_ldw)
  for (page in seq_len(num_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx   <- min(page * plots_per_page, n)
    file_name <- file.path(folder, paste0(prefix, "_", page, ".pdf"))
    pdf(file = file_name, width = 10, height = 12)
    par(mfrow = c(plots_per_page, 1), mar = c(4.5, 5, 3, 2))
    for (i in start_idx:end_idx) {
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

#### save_main_catt_panels()
save_main_catt_panels <- function(
    catt_refs,
    catt_comps, 
    ylabels,
    prefix = "catt_main", 
    plots_per_page = 4, 
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

#### save_plus_catt_panels()
save_plus_catt_panels <- function(
    catt1_list, 
    catt2_list, 
    ylabels, 
    prefix = "catt_plus", 
    plots_per_page = 4, 
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

## 4.3 QTET
#### plot_qte_top()
plot_qte_top <- function(qtet_top, qtet_top0, bm, plot_titles, main_start = 1, 
                         ylim = c(-25000, 15000), col = NULL) {
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

#### save_qte_top()
save_qte_top <- function(qtet_top, qtet_top0, bm, plot_titles, main_start = 1,
                                ylim = c(-25000, 15000), col = NULL, prefix = "ldw_top") {
  n <- length(qtet_top)
  dir.create("graphs/lalonde", showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(n)) {
    mod <- qtet_top[[i]]
    mod2 <- qtet_top0[[i]]
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

## 4.4 Assessing outcome weights (OW)
#### get_res_att()
get_res_att <- function(dataset_list, Y, treat, covar,
                            estimator = "AIPW_ATT",
                            smoother = "honest_forest",
                            n_cf_folds = 5,
                            n_reps = 1) {
  lapply(dataset_list,
      function(data) {
        dml_with_smoother(
          Y = data[[Y]],
          D = data[[treat]],
          X = data[, covar, drop = FALSE],
          estimator = estimator,
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

# 5. Sensitivity Analysis
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