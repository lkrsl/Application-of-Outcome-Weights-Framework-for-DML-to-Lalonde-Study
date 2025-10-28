# 1 Set up
# 1.1 Installation
packages <- c("cobalt", "data.table", "dplyr", "ebal", "ggplot2", "gridExtra", "highr", "highs", 
              "kableExtra", "MatchIt", "OutcomeWeights", "optmatch", "optweight", "quickmatch", 
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
library(cobalt)
library(data.table)
library(dplyr)
library(ebal)
library(ggplot2)
library(gridExtra)
library(highr)
library(highs)
library(kableExtra)
library(MatchIt)
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
# 2.1 Matching
# NA

# 2.2 Weighting
# NA

## 2.3 Truncation
### 2.3.1 Fixed maximum value truncation
#### truncate_weights_fixed()
truncate_weights_fixed <- function(data, weight_col, lower = 0.025, upper = 0.975) {
  data[[weight_col]] <- pmax(data[[weight_col]], lower)
  data[[weight_col]] <- pmin(data[[weight_col]], upper)
  return(data)
}

### 2.3.2 At percentile truncation
#### truncate_weights_percentile()
truncate_weights_percentile <- function(data, weight_col, lower = 0.05, upper = 0.95) {
  quantiles <- quantile(data[[weight_col]], probs = c(lower, upper), na.rm = TRUE)
  data[[weight_col]] <- pmax(data[[weight_col]], quantiles[1])
  data[[weight_col]] <- pmin(data[[weight_col]], quantiles[2])
  return(data)
} 

### 2.3.3 Adaptive weight truncation
#### check_weights()
check_weights <- function(data, weight_col = "weight") {
  w <- data[[weight_col]]
  variance <- var(w, na.rm = TRUE)
  data.frame(
    Weight_Column = weight_col,
    Variance = variance
  )
}

#### truncate_weights_adaptive()
truncate_weights_adaptive <- function(data, weight_col, c = 3) {
  w <- data[[weight_col]]
  # only apply if variance is above zero
  if (var(w, na.rm = TRUE) > 0) {
    cutoff <- mean(w, na.rm = TRUE) + c * sd(w, na.rm = TRUE)
    data[[weight_col]] <- pmin(w, cutoff)
  }
  return(data)
}

## 2.4 Trimming
### ps_estimate() ***
ps_estimate <- function(data, Y, treat, covar, num.trees = 4000, seed = 42) {
  set.seed(seed)
  data$ps_estimate <- probability_forest(
    X = data[, covar],
    Y = as.factor(data[, treat]),
    seed = seed, num.trees = num.trees
  )$predictions[, 2]
  return(data)
}

### 2.4.1 Propensity score threshold trimming (similar to tutorial of Imbens & Xu (2024))
# ps_trim() *** 
ps_trim <- function(data, ps = "ps_assoverlap", threshold = 0.9) { 
  sub <- data[which(data[, ps] < threshold), ]
  return(sub)
}

### 2.4.2 Common range trimming
#### common_range_trim ()
common_range_trim <- function(data, ps = "ps_assoverlap", treat = "treat") {
  lower_cut <- max(
    min(data[[ps]][data[[treat]] == 1], na.rm = TRUE),
    min(data[[ps]][data[[treat]] == 0], na.rm = TRUE)
  )
  upper_cut <- min(
    max(data[[ps]][data[[treat]] == 1], na.rm = TRUE),
    max(data[[ps]][data[[treat]] == 0], na.rm = TRUE)
  )
  sub <- data[data[[ps]] >= lower_cut & data[[ps]] <= upper_cut, ]
  return(sub)
}

### 2.4.3 Crump trimming
#### crump_trim ()
crump_trim <- function(data, ps = "ps_assoverlap", lower = 0.1, upper = 0.9) {
  sub <- data[data[[ps]] >= lower & data[[ps]] <= upper, ]
  return(sub)
}

### 2.4.4 Stuermer trimming
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

### 2.4.5 Walker trimming
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

## 2.5 Combination of methods
#### trim_attach_weights()
trim_attach_weights <- function(trimmed_list, model, weight_type = c("ipw_weight", "opt_weight", "cbps_weight", "ebal_weight"), estimand = "ATT") {
  weight_type <- match.arg(weight_type)
  weight_fun <- function(data) {
    if (weight_type == "ipw_weight") {
      w.obj <- WeightIt::weightit(model, data = data, estimand = estimand, method = "glm")
      data$ipw_weight <- w.obj$weights
    } else if (weight_type == "opt_weight") {
      w.obj <- optweight::optweight(model, data = data, estimand = estimand)
      data$opt_weight <- w.obj$weights
    } else if (weight_type == "cbps_weight") {
      w.obj <- WeightIt::weightit(model, data = data, estimand = estimand, method = "cbps")
      data$cbps_weight <- w.obj$weights
    } else if (weight_type == "ebal_weight") {
      w.obj <- WeightIt::weightit(model, data = data, estimand = estimand, method = "ebal")
      data$ebal_weight <- w.obj$weights
    }
    return(data)
  }
  weighted_list <- lapply(trimmed_list, weight_fun)
  return(weighted_list)
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
#### plot_matching_balance()***
plot_matching_balance <- function(matchit_objects, threshold = 0.1, title = NULL) {
  if (!is.list(matchit_objects)) {
    matchit_objects <- list(single = matchit_objects)
  }
  plot_list <- list()
  for (name in names(matchit_objects)) {
    matchit_object <- matchit_objects[[name]]
    bal <- cobalt::bal.tab(matchit_object, un = TRUE)
    smd_df <- data.frame(
      Variable = rownames(bal$Balance),
      Pre = abs(bal$Balance[, "Diff.Un"]),
      Post = abs(bal$Balance[, "Diff.Adj"])
    )
    # remove propensity score or non-covariate rows if present
    smd_df <- smd_df[!grepl("^distance$|^prop.score$|distance|prop.score", smd_df$Variable), ]
    # convert to long format
    smd_long <- pivot_longer(
      smd_df,
      cols = c("Pre", "Post"),
      names_to = "Matching",
      values_to = "Std_Diff"
    )
    smd_long$Matching <- factor(
      smd_long$Matching, 
      levels = c("Pre", "Post"), 
      labels = c("Pre-Matching", "Post-Matching")
    )
    # build titles
    plot_title <- if (length(matchit_objects) > 1) paste(title, "-", name) else title
    # create ggplot 
    p <- ggplot(smd_long, aes(x = Variable, y = Std_Diff, color = Matching)) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = -threshold, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = plot_title, x = "Covariates", y = "Absolute Standardized Mean Differences") +
      theme_minimal() +
      theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
            plot.title = element_text(face = "plain", size = 11, hjust = 0.5)) +
      scale_color_manual(values = c("Pre-Matching" = "#00CFC1", "Post-Matching" = "#FC766A"))
    plot_list[[name]] <- p
  }
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}

## 3.2 Weighting
### 3.2.1 SMD
#### compute_abs_smd_weight()
compute_abs_smd_weight <- function(data, treat, covar, weights_list) {
  smd_list <- lapply(names(weights_list), function(method) {
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = weights_list[[method]],
      un = TRUE,
      s.d.denom = "treated"
    )
    smd_df <- as.data.frame(bal_obj$Balance)
    smd_vals <- abs(smd_df$Diff.Adj)
    mean_smd <- mean(smd_vals, na.rm = TRUE)
    max_smd  <- max(smd_vals, na.rm = TRUE)
    return(data.frame(
      Method = method,
      Mean_Abs_SMD = mean_smd,
      Max_Abs_SMD  = max_smd
    ))
  })
  do.call(rbind, smd_list)
}

### 3.2.2 ESS
#### compute_ess_weight()
compute_ess_weight <- function(data, treat, covar, weights_list) {
  ess_list <- lapply(names(weights_list), function(method) {
    bal_obj <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = weights_list[[method]],
      un = FALSE
    )
    samples <- bal_obj$Observations
    if ("Adjusted" %in% rownames(samples)) {
      df <- samples["Adjusted", c("Control", "Treated"), drop = FALSE]
    } else {
      df <- samples[1, c("Control", "Treated"), drop = FALSE]
    }
    df <- cbind(Method = method, df)
    rownames(df) <- NULL
    return(df)
  })
  do.call(rbind, ess_list)
}

### 3.2.3 Visuals
#### plot_weighting_balance()***
plot_weighting_balance <- function(data, treat, covar, weight_list, title = NULL) {
  plot_list <- list()
  for (wname in names(weight_list)) {
    invisible(capture.output({
    bal <- cobalt::bal.tab(
      as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
      data = data,
      weights = weight_list[[wname]],
      un = TRUE,
      s.d.denom = "treated"
    )
    smd_df <- data.frame(
      Variable = rownames(bal$Balance),
      Pre  = abs(bal$Balance[, "Diff.Un"]),
      Post = abs(bal$Balance[, "Diff.Adj"])
    )
    smd_df <- smd_df[!grepl("^distance$|^prop.score$|distance|prop.score", smd_df$Variable), ]
    smd_long <- pivot_longer(
      smd_df,
      cols = c("Pre", "Post"),
      names_to = "Weighting",
      values_to = "Std_Diff"
    )
    smd_long$Weighting <- factor(
      smd_long$Weighting,
      levels = c("Pre", "Post"),
      labels = c("Pre-Weighting", "Post-Weighting")
    )
    plot_title <- if (length(weight_list) > 1) paste(title, "-", wname) else title
    p <- ggplot(smd_long, aes(x = Variable, y = Std_Diff, color = Weighting)) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
      geom_hline(yintercept = -0.1, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = plot_title, x = "Covariates", y = "Absolute Standardized Mean Differences") +
      theme_minimal() +
      theme(panel.border = element_rect(color = "black", fill = NA, size = 1) ,
            plot.title = element_text(face = "plain", size = 11, hjust = 0.5)) +
      scale_color_manual(values = c("Pre-Weighting" = "#00CFC1", "Post-Weighting" = "#FC766A"))
    plot_list[[wname]] <- p
   }))
  }  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
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
#### plot_trunc_balance()
plot_trunc_balance <- function(trunc_list, treat, covar, weight_cols, dataset_name = NULL, threshold = 0.1) {
  plot_list <- list()
  for (trunc_name in names(trunc_list)) {
    dataset <- trunc_list[[trunc_name]]
    for (wcol in weight_cols) 
      invisible(capture.output({
      if (wcol %in% names(dataset)) {
        bal_obj <- cobalt::bal.tab(
          as.formula(paste(treat, "~", paste(covar, collapse = " + "))),
          data = dataset,
          weights = dataset[[wcol]],
          un = TRUE,
          s.d.denom = "treated"
        )
        smd_df <- data.frame(
          Variable = rownames(bal_obj$Balance),
          Pre = abs(bal_obj$Balance[,"Diff.Un"]),
          Post = abs(bal_obj$Balance[,"Diff.Adj"])
        )
        smd_df <- smd_df[!grepl("^distance$|^prop.score$|distance|prop.score", smd_df$Variable), ]
        smd_long <- pivot_longer(
          smd_df,
          cols = c("Pre", "Post"),
          names_to = "Truncation",
          values_to = "Std_Diff"
        )
        smd_long$Truncation <- factor(
          smd_long$Truncation,
          levels = c("Pre", "Post"),
          labels = c("Pre-Truncation", "Post-Truncation")
        )
        plot_title <- paste0(
          if (!is.null(dataset_name)) paste0(dataset_name, " - ") else "",
          trunc_name, "_", wcol, " truncation"
        )
        p <- ggplot(smd_long, aes(x = Variable, y = Std_Diff, color = Truncation)) +
          geom_point(size = 3) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
          geom_hline(yintercept = -threshold, linetype = "dashed", color = "red") +
          coord_flip() +
          labs(title = plot_title, x = "Covariates", y = "Absolute Standardized Mean Differences") +
          theme_minimal() +
          theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                plot.title = element_text(face = "plain", size = 11, hjust = 0.5)) +
          scale_color_manual(values = c("Pre-Truncation" = "#00CFC1", "Post-Truncation" = "#FC766A"))
        plot_list[[paste(trunc_name, wcol, sep = "_")]] <- p
      } else {
        message(sprintf("Column %s not found in %s", wcol, trunc_name))
      }
    }))
  }
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
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
#### plot_trim_overlap()
plot_trim_overlap <- function(data_list, treat, covar, prefix = NULL, 
                              main.size = 1.1, lab.size = 1, axis.size = 1) {
  par(mfrow = c(2, 3), mar = c(4, 4, 1, 1),
      cex.lab = lab.size, font.lab = 1,
      cex.axis = axis.size, font.axis = 1)
  for (name in names(data_list)) {
    data_obj <- data_list[[name]]
    invisible(capture.output(
      assess_overlap(data_obj, treat = treat, cov = covar)
    ))
    title(main = paste0(prefix, " - ", name), cex.main = main.size, font.main = 1)
  }
}

## 3.5 Combined methods
### 3.5.1 SMD
#### compute_smd_comb()
compute_abs_smd_comb <- function(combined_list, treat, covar) {
  smd_list <- lapply(names(combined_list), function(weight_method) {
    method_list <- combined_list[[weight_method]]
    res <- lapply(names(method_list), function(trim_method) {
      data <- method_list[[trim_method]]
      # use equal weights if weight column missing
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
#### compute_ess_comb()
compute_ess_comb <- function(combined_list, treat, covar) {
  ess_list <- lapply(names(combined_list), function(weight_method) {
    method_list <- combined_list[[weight_method]]
    res <- lapply(names(method_list), function(trim_method) {
      data <- method_list[[trim_method]]
      # use equal weights if weight column missing
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
                              prefix_cps = NULL, prefix_psid = NULL,
                              main.size = 1.1, lab.size = 1, axis.size = 1) {
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
        plot_list[[full_method_name]] <- combined_list[[weight_method]][[trim_method]]
        method_names <- c(method_names, full_method_name)
      }
    }
    total_plots <- length(plot_list)
    plots_per_page <- 6
    for (i in seq(1, total_plots, by = plots_per_page)) {
      start_idx <- i
      end_idx <- min(i + plots_per_page - 1, total_plots)
      par(mfrow = c(2, 3), mar = c(4, 4, 1, 1), 
          cex.lab = lab.size, font.lab = 1,
          cex.axis = axis.size, font.axis = 1)
      for (j in start_idx:end_idx) {
        data <- plot_list[[j]]
        prefix <- if (ds_name == "CPS") prefix_cps else prefix_psid
        invisible(capture.output(
          assess_overlap(data, treat = treat, cov = covar) 
        ))
        title(main = paste0(prefix, " - ", method_names[j]), cex.main = main.size, font.main = 1)
      }
      if ((end_idx - start_idx + 1) < plots_per_page) {
        for (k in seq_len(plots_per_page - (end_idx - start_idx + 1))) plot.new()
      }
    }
  }
}

#### save_comb_hist()
save_comb_hist <- function(comb_meth_cps = NULL, comb_meth_psid = NULL, treat, covar,
                           prefix = NULL,
                           prefix_cps = NULL, prefix_psid = NULL,
                           path = "graphs/lalonde") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  print(getwd())
  all_combined_list <- list()
  if (!is.null(comb_meth_cps)) all_combined_list$CPS <- comb_meth_cps
  if (!is.null(comb_meth_psid)) all_combined_list$PSID <- comb_meth_psid
  
  file_index <- 1
  for (ds_name in names(all_combined_list)) {
    combined_list <- all_combined_list[[ds_name]]
    plot_list <- list()
    method_names <- character()
    for (weight_method in names(combined_list)) {
      for (trim_method in names(combined_list[[weight_method]])) {
        full_method_name <- paste(weight_method, trim_method, sep = "_")
        plot_list[[full_method_name]] <- combined_list[[weight_method]][[trim_method]]
        method_names <- c(method_names, full_method_name)
      }
    }
    total_plots <- length(plot_list)
    plots_per_page <- 4
    for (i in seq(1, total_plots, by = plots_per_page)) {
      start_idx <- i
      end_idx <- min(i + plots_per_page - 1, total_plots)
      pdf_file <- file.path(path, sprintf("%s_ov_%d.pdf", prefix, file_index))
      pdf(pdf_file, width = 10, height = 8)
      par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), cex.lab = 1, font.lab = 1, cex.axis = 1, font.axis = 1)
      for (j in start_idx:end_idx) {
        data <- plot_list[[j]]
        prefix_str <- if (ds_name == "CPS") prefix_cps else prefix_psid
        assess_overlap(data, treat = treat, cov = covar)
        title(main = paste0(prefix_str, " - ", method_names[j]), cex.main = 1, font.main = 1)
      }
      if ((end_idx - start_idx + 1) < 4) {
        for (k in seq_len(4 - (end_idx - start_idx + 1))) plot.new()
      }
      dev.off()
      file_index <- file_index + 1
    }
  }
}

## 3.6 Getting top methods and datasets
#### combine_results()
combine_results <- function(dataset_name) {
  dataset_lower <- tolower(dataset_name)
  # retrieve matching results
  smd_matching <- get(paste0("smd_matchit.", dataset_lower))
  ess_matching <- get(paste0("ess_matchit.", dataset_lower))
  # retrieve trimming results
  smd_trimming <- get(paste0("smd_trim.", dataset_lower))
  ess_trimming <- get(paste0("ess_trim.", dataset_lower))
  # retrieve truncation results
  smd_trunc <- get(paste0("smd_trunc.", dataset_lower))
  ess_trunc <- get(paste0("ess_trunc.", dataset_lower))
  # retrieve weighting results
  smd_weighting <- get(paste0("smd_weight.", dataset_lower))
  ess_weighting <- get(paste0("ess_weight.", dataset_lower))
  # retrieve combined SMD and ESS
  smd_combined_var <- paste0("smd_all_comb_meth.", dataset_lower)
  ess_combined_var <- paste0("ess_all_comb_meth.", dataset_lower)
  smd_combined <- get(smd_combined_var)
  ess_combined <- get(ess_combined_var)
  # combine all SMD results
  smd_all <- do.call(rbind, list(
    smd_matching,
    smd_trimming,
    smd_trunc,
    smd_weighting,
    smd_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")]
  ))
  # combine all ESS results
  ess_all <- do.call(rbind, list(
    ess_matching,
    ess_trimming,
    ess_trunc,
    ess_weighting,
    ess_combined[, c("Method", "Control", "Treated")]
  ))
  # merge SMD and ESS results by Method
  final_df <- merge(smd_all, ess_all, by = "Method", all = TRUE)
  # remove dataset suffixes for clean labels
  final_df$Method <- gsub("\\.psid", "", final_df$Method, ignore.case = TRUE)
  final_df$Method <- gsub("\\.cps", "", final_df$Method, ignore.case = TRUE)
  # reset row names
  rownames(final_df) <- NULL
  return(final_df)
}

#### combine_results_plus()
combine_results_plus <- function(dataset_name) {
  dataset_lower <- tolower(dataset_name)
  # retrieve standard results
  smd_matching <- get(paste0("smd_matchit.", dataset_lower))
  ess_matching <- get(paste0("ess_matchit.", dataset_lower))
  smd_trimming <- get(paste0("smd_trim.", dataset_lower))
  ess_trimming <- get(paste0("ess_trim.", dataset_lower))
  smd_trunc <- get(paste0("smd_trunc.", dataset_lower))
  ess_trunc <- get(paste0("ess_trunc.", dataset_lower))
  smd_weighting <- get(paste0("smd_weight.", dataset_lower))
  ess_weighting <- get(paste0("ess_weight.", dataset_lower))
  smd_combined <- get(paste0("smd_all_comb_meth.", dataset_lower))
  ess_combined <- get(paste0("ess_all_comb_meth.", dataset_lower))
  # optionally retrieve plus dataset results if they exist
  smd_combined_plus_name <- paste0("smd_all_comb_meth.", dataset_lower, "_plus")
  ess_combined_plus_name <- paste0("ess_all_comb_meth.", dataset_lower, "_plus")
  smd_combined_plus <- if (exists(smd_combined_plus_name)) get(smd_combined_plus_name) else NULL
  ess_combined_plus <- if (exists(ess_combined_plus_name)) get(ess_combined_plus_name) else NULL
  # combine all SMD results
  smd_list <- list(
    smd_matching,
    smd_trimming,
    smd_trunc,
    smd_weighting,
    smd_combined[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")]
  )
  if (!is.null(smd_combined_plus)) {
    smd_list <- append(smd_list, list(smd_combined_plus[, c("Method", "Mean_Abs_SMD", "Max_Abs_SMD")]))
  }
  smd_all <- do.call(rbind, smd_list)
  # combine all ESS results
  ess_list <- list(
    ess_matching,
    ess_trimming,
    ess_trunc,
    ess_weighting,
    ess_combined[, c("Method", "Control", "Treated")]
  )
  if (!is.null(ess_combined_plus)) {
    ess_list <- append(ess_list, list(ess_combined_plus[, c("Method", "Control", "Treated")]))
  }
  ess_all <- do.call(rbind, ess_list)
  # merge SMD and ESS results
  final_df <- merge(smd_all, ess_all, by = "Method", all = TRUE)
  # clean method names
  final_df$Method <- gsub("\\.psid|\\.cps", "", final_df$Method, ignore.case = TRUE)
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

#### assess_methods
assess_methods <- function(data) {
  data %>%
    mutate(
      # ess score
      ess_balance_ratio = pmin(Control, Treated) / pmax(Control, Treated),
      ess_total = Control + Treated,
      # normalize balance_ratio and ess_total
      ess_balance_score = (ess_balance_ratio - min(ess_balance_ratio, na.rm = TRUE)) /
        (max(ess_balance_ratio, na.rm = TRUE) - min(ess_balance_ratio, na.rm = TRUE)),
      ess_size_score = (ess_total - min(ess_total, na.rm = TRUE)) /
        (max(ess_total, na.rm = TRUE) - min(ess_total, na.rm = TRUE)),
      # combine size and balance equally
      ess_score = 0.5 * ess_balance_score + 0.5 * ess_size_score,
      # smd score
      smd_score = 1 - (Mean_Abs_SMD - min(Mean_Abs_SMD, na.rm = TRUE)) /
        (max(Mean_Abs_SMD, na.rm = TRUE) - min(Mean_Abs_SMD, na.rm = TRUE)),
      # composite score
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