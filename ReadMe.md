Replication Files for Kreisel (2025)
================

This Github repository contains the data and code necessary to replicate
**Kreisel (2025)**: "Application of the Outcome Weights Framework for Double Machine Learning to the LaLonde Study (2025)"

**Tutorial**: For a detailed bookdown tutorial, see
[here](https://github.com/lkrsl/Application-of-Outcome-Weights-Framework-for-DML-to-Lalonde-Study). 

## Folder Structure

The folder structure of this repository is as follows:

| folder | usage                                                      |
|:-------|:-----------------------------------------------------------|
| bookdown| Bookdown files                                             |
| code    | R scripts                                                  |
| data    | Data files, including four subfolders “ldw_model_a”, “ldw_model_b”, “nsw” & “lcs” |
| docs    | HTML documentation and website files from bookdown rendering |
| graphs  | Graphics including one subfolder "lalonde"                   |
| package | Local R package including "OutcomeWeights" package source code |
| tables  | CSV files                                           |
| tutorial  | Bookdown RMD files, functions, bibliography, data and style sheets   |

## Data Files

Kreisel (2025) uses the following datasets, which are based on
LaLonde (1986), Dehejia and Wahba (1999) and Calónico and Smith (2017).

| Data.files        | Details                                                          | File_Type | Experimental |
|:------------------|:-----------------------------------------------------------------|:----------|:-------------|
| nsw.dta           | NSW experimental data, used in LaLonde (1986)                    | Stata     | Yes          |
| nsw_dw.dta        | Subset of NSW experimental data, used in Dehejia & Wahba (1999)  | Stata     | Yes          |
| cps_controls.dta  | CPS-SSA-1 controls, used in both papers                          | Stata     | No           |
| psid_controls.dta | PSID-1 controls, used in both papers                             | Stata     | No           |
| NSW_AFDC_CS.dta   | Reconstructed NSW AFDC female samples                            | Stata     | Both         |

## Install R Packages

Successful replication requires installation of the following R packages:

``` r
# required packages
packages <- c("CBPS", "cobalt", "data.table", "dplyr", "DT", "DoubleML", "ebal", "estimatr", "ggplot2", "gridExtra", "grf",
              "hbal", "highr", "highs", "kableExtra", "MatchIt", "Matching", "mlr3", "mlr3learners", "OutcomeWeights",
              "optmatch", "optweight", "quickmatch", "readr", "rgenoud", "tidyr", "tidyverse", "patchwork", "ppcor", "WeightIt"
)

# install packages
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}

install_all(packages)
```
## Report Errors

To report errors, please contact <laura.kreisel@student.uni-tuebingen.de>. Any suggestions or comments are welcome.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-dehejiawahba" class="csl-entry">

Dehejia, Rajeev H, and Sadek Wahba. 1999. “Causal Effects in
Nonexperimental Studies: Reevaluating the Evaluation of Training
Programs.” *Journal of the American Statistical Association* 94 (448):
1053–62.

</div>

<div id="ref-imbensrubinsacerdote" class="csl-entry">

Imbens, Guido W, Donald B Rubin, and Bruce I Sacerdote. 2001.
“Estimating the Effect of Unearned Income on Labor Earnings, Savings,
and Consumption: Evidence from a Survey of Lottery Players.” *American
Economic Review* 91 (4): 778–94.

</div>

<div id="ref-imbensxu" class="csl-entry">

Imbens, Guido W, and Yiqing Xu. 2024. “LaLonde (1986) After Nearly Four
Decades: Lessons Learned.” arXiv:2406.00827.

</div>

<div id="ref-LaLonde" class="csl-entry">

LaLonde, Robert J. 1986. “Evaluating the Econometric Evaluations of
Training Programs with Experimental Data.” *The American Economic
Review* 76 (4): 604–20.

</div>

<div id="ref-imbensxu" class="csl-entry">
  
Xie, Yihui. 2015. Dynamic Documents with R and knitr. 2nd ed. Boca Raton, Florida: Chapman and Hall/CRC. ISBN 978-1498716963. http://yihui.org/knitr/.

</div>

<div id="ref-outcomeweightspackage" class="csl-entry">
  
Knaus, Michael C., and Henri Pfleiderer. 2024. “Outcome Weights.” Retrieved from https://cran.r-project.org/web/packages/OutcomeWeights/OutcomeWeights.pdf.
</div>

<div id="ref-greifer2025" class="csl-entry">
  
Greifer, Noah. 2025. “Assessing Balance.” Retrieved from https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html.

</div>

<div id="ref-greifer2_2025" class="csl-entry">
  
Greifer, Noah. 2025. “Estimating Effects After Matching.” Retrieved from https://cran.r-project.org/web/packages/MatchIt/vignettes/estimating-effects.html.

</div>

<div id="ref-greifer3_2025" class="csl-entry">
  
Greifer, Noah. 2025. “Matching Methods.” Retrieved from https://cran.r-project.org/web/packages/MatchIt/vignettes/matching-methods.html.

</div>

<div id="ref-sekhon2011" class="csl-entry">
  
Sekhon, Jasjeet S. 2011. “Multivariate and Propensity Score Matching Software with Automated Balance Optimization: The Matching Package for R.” Journal of Statistical Software.

</div>

<div id="ref-sharon2024" class="csl-entry">
  
Green, Sharon. 2024. “Propensity Score Matching for Causal Inference: Creating Data Visualizations to Assess Covariate Balance in R.” Retrieved from https://dlab.berkeley.edu/news/propensity-score-matching-causal-inference-creating-data-visualizations-assess-covariate.

</div>

<div id="ref-bryer2025" class="csl-entry">
  
Bryer, Jason. 2025. “Applied Propensity Score Analysis with R.” Retrieved from https://psa.bryer.org/chapter-weighting.html.

</div>

<div id="ref-larson2023" class="csl-entry">

Larson, Dirk R., Isabella Zaniletti, David G. Lewallen, Daniel J. Berry, and Hilal Maradit Kremers. 2023. “Propensity Scores: Confounder Adjustment When Comparing Nonrandomized Groups in Orthopaedic Surgery.” The Journal of Arthroplasty 38: 622–26.

</div>

<div id="ref-crump2009" class="csl-entry">

Crump, Richard K., V. Joseph Hotz, Guido W. Imbens, and Oscar A. Mitnik. 2009. “Dealing with Limited Overlap in Estimation of Average Treatment Effects.” Biometrika 96 (1): 187–99.

</div>

<div id="ref-stuermer2021" class="csl-entry">

Stürmer, Til, Michael Webster-Clark, Jennifer L. Lund, Richard Wyss, Alan R. Ellis, Mark Lunt, Kenneth J. Rothman, and Robert J. Glynn. 2021. “Propensity Score Weighting and Trimming Strategies for Reducing Variance and Bias of Treatment Effect Estimates: A Simulation Study.” American Journal of Epidemiology.

</div>

<div id="ref-stuermer2010" class="csl-entry">

Stürmer, Til, Kenneth J. Rothman, Jerry Avorn, and Robert J. Glynn. 2010. “Treatment Effects in the Presence of Unmeasured Confounding: Dealing With Observations in the Tails of the Propensity Score Distribution—A Simulation Study.” American Journal of Epidemiology.

</div>

<div id="ref-matsouaka2023" class="csl-entry">

Matsouaka, Roland A., and Yunji Zhou. 2023. “Causal Inference in the Absence of Positivity: The Role of Overlap Weights.” Biometrical Journal.

</div>

<div id="ref-Firpo2007" class="csl-entry">

Firpo, Sergio. 2007. “Efficient Semiparametric Estimation of Quantile Treatment Effects.” Econometrica 75 (1): 259-276 

</div>
