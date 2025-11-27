# AGENTS.md - R Bookdown Project

## Build & Test Commands
- Render bookdown: `Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"`
- Render single chapter: `Rscript -e "rmarkdown::render('tutorial/02-LDW-Model-A.Rmd')"`
- Test OutcomeWeights package: `Rscript -e "testthat::test_file('package/OutcomeWeights/tests/test-get_outcome_weights.R')"`
- Install packages: `Rscript -e "source('tutorial/functions.R')"`

## Architecture
- **Bookdown**: Multi-chapter tutorial built from index.Rmd + tutorial/*.Rmd files; output to `docs/`
- **Data**: Three datasets (LDW, NSW, LCS) in `data/` with subfolders for model variants
- **Package**: Local `OutcomeWeights` R package in `package/OutcomeWeights/` with custom DML estimators
- **Code**: Data preparation scripts in `code/` (e.g., lcs_prepare.R)
- **Outputs**: Graphics in `graphs/`, tables in `tables/` (CSV format)

## Code Style
- **Imports**: Load all packages explicitly via library() at start; no :: operators mid-code
- **Naming**: Snake_case for functions and variables (e.g., `compute_abs_smd_matchit()`, `lcs_psid`)
- **Data frames**: Use base R data.frame; apply() for operations; minimal tidyverse except dplyr/ggplot2
- **Seed**: Always set seed before randomized operations (e.g., `set.seed(1234)`)
- **Paths**: Use forward slashes in file.path(); data loaded from `data/`, tables saved to `tables/`, graphs to `graphs/`
- **Comments**: Section headers with `####`; inline comments explain statistical choices (e.g., threshold = 0.9)
