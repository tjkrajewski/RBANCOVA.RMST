RBANCOVA.RMST: RBANCOVA for RMST estimates with Categorized
time-to-event data
================

``` r
# During development, you can uncomment the next line in an interactive session
# to use the functions in R/ without installing the package:
# devtools::load_all("..")
devtools::install()
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> WARNING: Rtools is required to build R packages, but is not currently installed.
#> 
#> Please download and install Rtools 4.5 from https://cran.r-project.org/bin/windows/Rtools/.
#>          checking for file 'C:\Users\tjk42\Documents\GitHub\RBANCOVA.RMST/DESCRIPTION' ...  ✔  checking for file 'C:\Users\tjk42\Documents\GitHub\RBANCOVA.RMST/DESCRIPTION' (500ms)
#>       ─  preparing 'RBANCOVA.RMST':
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts (459ms)
#>   ─  checking for empty or unneeded directories
#>       ─  building 'RBANCOVA.RMST_0.0.0.9000.tar.gz'
#>      
#> 
#> WARNING: Rtools is required to build R packages, but is not currently installed.
#> 
#> Please download and install Rtools 4.5 from https://cran.r-project.org/bin/windows/Rtools/.
#> Running "C:/PROGRA~1/R/R-45~1.1/bin/x64/Rcmd.exe" INSTALL \
#>   "C:\Users\tjk42\AppData\Local\Temp\Rtmp2b4iyf/RBANCOVA.RMST_0.0.0.9000.tar.gz" \
#>   --install-tests 
#> * installing to library 'C:/Users/tjk42/AppData/Local/R/win-library/4.5'
#> * installing *source* package 'RBANCOVA.RMST' ...
#> ** this is package 'RBANCOVA.RMST' version '0.0.0.9000'
#> ** using staged installation
#> ** R
#> ** byte-compile and prepare package for lazy loading
#> ** help
#> *** installing help indices
#> ** building package indices
#> ** installing vignettes
#> ** testing if installed package can be loaded from temporary location
#> ** testing if installed package can be loaded from final location
#> ** testing if installed package keeps a record of temporary installation path
#> * DONE (RBANCOVA.RMST)

library(RBANCOVA.RMST)
```

# Overview

This vignette shows how to implement Randmoziation-Based Covariance
Analysis (RB_ANCOVA) for Restricted Mean Survival Time (RMST) contrasts
using RBANCOVA.RMST. We’ll demonstrate via application to estimating the
effect of D-Penicillamine for Primary Biliary Cholangitis (PBC).

## APPLICATION: D-PENICILLAMINE FOR PRIMARY BILIARY CHOLANGITIS

Primary biliary cholangitis (PBC) is a chronic autoimmune disease
characterized by progressive inflammation and destruction of the small
bile ducts within the liver, ultimately leading to cirrhosis (permanent
scarring of liver tissues) and liver failure in advanced stages . The
progression of PBC is inevitable, but it can be prolonged with
treatment. Due to its clinical importance and the need for effective
therapeutic options, PBC has been a focus of clinical trials for several
decades. A dataset often used to study the natural history and treatment
response in PBC patients originates from a Mayo Clinic trial conducted
between 1974 and 1984.

The dataset includes 424 patients diagnosed with PBC who were referred
to the Mayo Clinic between 1974 and 1984 . Of these, 312 patients were
enrolled in a randomized, placebo-controlled trial of the drug
D-penicillamine, while the remaining 112 cases did not participate in
the clinical trial but consented to basic measurements and follow-up for
survival. The dataset encompasses clinical variables, including age,
sex, and a range of laboratory biomarkers such as serum bilirubin, serum
albumin, and alkaline phosphatase levels, alongside clinical indicators
like the presence of ascites, hepatomegaly, blood vessel malformations
(spider angiomata), and histologic stage of disease at entry.
Participants’ survival outcomes are tracked until death, loss to
follow-up (censored), liver transplantation (considered censored at time
of transplant), or censoring at the study endpoint in July 1986. The PBC
dataset can be accessed by installing the package in R with dataset .

This analysis focuses on the 312 trial participants (158 in the
D-penicillamine arm, 154 in the placebo arm). We implement
randomization-based covariate adjustment for RMST estimation, as
described in Section , to estimate treatment effects when adjusting for
a selection of covariates.

The methodology described in Section was applied to the PBC dataset. The
Kaplan-Meier curves illustrate the time to event for participants in the
placebo and D-penicillamine groups over a 12.5-year follow-up period
(Figure ). Over the first five years, there is minimal separation in
survival probabilites between the two arms, with the D-penicillamine arm
showing slightly higher survival probabilities compared to the placebo
group. However, starting around five years, the Kaplan-Meier survival
curves demonstrate that participants in the placebo arm generally had a
higher probability of survival compared to those in the D-penicillamine
group until approximately year nine or ten. %The separation between the
curves becomes more noticeable around 5 years, indicating a potential
survival disadvantage for the D-penicillamine group.

<!-- \begin{figure} -->

<!-- \centering -->

<!-- \includegraphics[width=.8\linewidth]{Figures/pbc_survival_KM.png} -->

<!-- \caption{Kaplan-Meier (KM) estimate for time-to-death of each treatment group (D-penicillamine or placebo) for patients enrolled in the PBC Mayo Clinic trial. Time represents years from trial registration to death or last date known alive.}\label{fig:pbc_km} -->

<!-- \end{figure} -->

## Single Timepoint

To quantify the observed difference, we calculated the RMST for each
group from the start of the trial until year ten using ten yearly
intervals. Time was restricted to ten years of follow-up as many
individuals at risk after year ten were censored before the end of
follow-up. This yielded an RMST of 7.13 years (SE 0.28) for the
D-penicillamine group and 7.28 years (SE 0.30) for the placebo group.
%These results estimate that the average survival time during the first
ten years of follow-up was 7.13 years in the D-penicillamine group and
7.28 years in the placebo group. The unadjusted difference in RMST
(D-penicillamine - placebo) was estimated to be -0.15 years with a
confidence interval of \[-0.95, 0.65\] (Figure ).

To further refine the estimated difference in RMST, we employed
RB-ANCOVA to adjust for covariates selected via stepwise selection. The
covariates retained for adjustment were histologic stage of disease,
copper (ug/day), serum bilirubin (mg/dl), and serum albumin (g/dl),
reflecting clinically relevant factors that might impact survival
outcomes. Baseline balance between the trial arms was similar for most
of these covariates (Table ): the average histoloic stage was slightly
higher in the placebo group (3.09 vs. 2.97), while serum bilirubin
levels were lower in the D-penicillamine group (2.87 vs. 3.65 mg/dl).
Copper levels and serum albumin were comparable between groups,
indicating minimal baseline differences. Assessing the extent of random
imbalance between the two treatment arms, as described in Section gave
no evidence of unexpected imbalance among the covariates (p-value =
0.23).

<!-- \begin{table} -->

<!-- \caption{Mean (SE) of selected variables by treatment arm} -->

<!-- \centering -->

<!-- \begin{tabular}{lcccc} -->

<!-- \toprule -->

<!-- \textbf{Trial Arm} & \textbf{Histologic Stage} & \textbf{Copper} & \textbf{Serum Bilirubin} & \textbf{Serum Albumin}\\ -->

<!-- \midrule -->

<!-- Placebo & 3.09 (0.07) & 97.65 (6.46) & 3.65 (0.43) & 3.52 (0.03)\\ -->

<!-- D-penicillamine & 2.97 (0.07) & 97.64 (7.18) & 2.87 (0.29) & 3.52 (0.04)\\ -->

<!-- \hiderowcolors -->

<!-- \hline  % Please only put a hline at the end of the table -->

<!-- \end{tabular}\label{table:pbc_covs} -->

<!-- \end{table} -->

After adjusting for these covariates using RB-ANCOVA, the adjusted RMST
difference was estimated to be -0.44 years with a confidence interval of
\[-1.05, 0.17\].  
While the RB-ANCOVA results are consistent with the unadjusted results,
the covariate adjustment improved the precision of the estimate. The CI
half-width for the unadjusted estimate was 0.80 versus 0.61 for the
adjusted estimate.

Figure also shows the results of implementing the continuous
methodology, unadjusted, as well as adjusted for covariates via ’s
ANCOVA-type method. Results when implementing the continuous methodology
are consistent with those found with our methodology. Adjusting for
covariates via RB-ANCOVA results in a similar estimate and confidence
interval to that of ’s method, but the RB-ANCOVA adjustment does not
require the assumptions of to hold.

<!-- \begin{figure} -->

<!-- \centering -->

<!-- \includegraphics[width=0.8\linewidth]{Figures/estimates_CIs_pbc_4covs.png} -->

<!-- \caption{Estimates and corresponding confidence intervals for difference in RMST of D-penicillamine vs. Placebo. Adjusted estimates adjust for histolic stage of disease, serum bilirubin (mg/dl), serum albumin (g/dl), and urine copper (ug/day). The Continuous-Adjusted results are shown from implementation of the continuous methodology, adjusting for covariates with \cite{Tian2018}'s ANCOVA-type method.}\label{fig:pbc_ci_1} -->

<!-- \end{figure} -->

## Multiple Timepoints

<!-- rmarkdown::render("vignettes/RMST-RBANCOVA.Rmd") -->
