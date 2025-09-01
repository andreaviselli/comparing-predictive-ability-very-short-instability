# Replication Package for: Comparing predictive ability in presence of instability over a very short time

**September 01, 2025**

This replication package accompanies Iacone F., Rossini L. and Viselli A. (2025). Comparing predictive ability in presence of instability over a very short time. *The Econometrics Journal*.

https://academic.oup.com/ectj/advance-article-abstract/doi/10.1093/ectj/utaf018/8238631?redirectedFrom=fulltext

## Data availability and provenance statements

### Statement about rights

The author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

### Summary of availability

All data are publicly available.

### Details on each data source

Data on the US Nominal Gross Domestic Product (NGDP) Survey of Professional Forecasters (SPF) were downloaded from the Federal Reserve Bank of Philadelphia. We use the Median Responses.

Data can be downloaded from https://www.philadelphiafed.org/surveys-and-data/ngdp, under "Billions of dollars. Seasonally adjusted. Annual rate. GNP prior to 1992. GDP 1992 to present", then download "Median responses".

Data can also be directly downloaded using https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/survey-of-professional-forecasters/data-files/files/Median_NGDP_Level.xlsx?sc_lang=en&hash=3D3DCDF99D4C8EB44EC5C7AD41B6BE47.

A copy of the data is provided as part of this archive. The data are in the public domain. Datafile: `SPF.mat` (Matlab, proprietary) and `SPF.csv` (CSV, non-proprietary).

### Details on the dataset

We use column (4) of the SPF data for Nominal GDP, labeled `NGDP2`, as specified in the SPF documentation available at: https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/spf-documentation.pdf?la=en&hash=F2D73A2CE0C3EA90E71A363719588D25

## Description of programs/code

Run the following programs/code in Matlab_R2023A (recommended).
The replication files are compatible with macOS, Windows, and Linux, provided that a compatible version of MATLAB is installed.
The code does not include any OS- or Windows-specific paths or functions, so no platform-dependent adjustments should be necessary.

*   `dmtest.m` includes the function to perform the Diebold and Mariano (DM, 1995) test.
*   `esitest.m` includes the function to perform the Andrews (2003) test.
*   `grtest.m` includes the function to perform the Giacomini and Rossi (GR, 2010) "Fluctuation" test.
*   `longrunvariance.m` includes the function to estimate the long-run variance for the DM and GR tests.
*   `maxtest.m` includes the function to perform the MAX test by Harvey et al. (2020).
*   `ReplicationCode_MC_Table_1.m` generates Table 1 in Section 3 (MONTE CARLO STUDIES).
*   `ReplicationCode_MC_Table_2.m` generates Table 2 in Section 3 (MONTE CARLO STUDIES).
*   `ReplicationCode_MC_Table_3.m` generates Table 3 in Section 3 (MONTE CARLO STUDIES).
*   `ReplicationCode_MC_Table_4.m` generates Table 4 in Section 3 (MONTE CARLO STUDIES).
*   `ReplicationCode_MC_Table_5.m` generates Table 5 in Section 3 (MONTE CARLO STUDIES).
*   `ReplicationCode_Table_6_7.m` generates Figure 1-2 and Tables 6-7 in Section 4 (EVALUATING THE NOWCAST OF US NOMINAL GDP GROWTH).

The replication codes for the Monte Carlo studies include a random seed.

## License for Code

The code is licensed under a MIT license. See license.txt for details.

## List of tables and figures

The provided code reproduces all tables and figures in the paper.

| Figure/Table # | Program | Page Number |
| :--- | :--- | :--- |
| Table 1 | `ReplicationCode_MC_Table_1.m` | 12 |
| Table 2 | `ReplicationCode_MC_Table_2.m` | 13 |
| Table 3 | `ReplicationCode_MC_Table_3.m` | 14 |
| Table 4 | `ReplicationCode_MC_Table_4.m` | 15 |
| Table 5 | `ReplicationCode_MC_Table_5.m` | 15 |
| Table 6 | `ReplicationCode_Table_6_7.m` | 17 |
| Table 7 | `ReplicationCode_Table_6_7.m` | 18 |
| Figure 1 | `ReplicationCode_Table_6_7.m` | 16 |
| Figure 2 | `ReplicationCode_Table_6_7.m` | 17 |

## Instructions to replicators

*   Set the working directory to the folder containing all MATLAB scripts and data files.
*   Execute a replication script (e.g. `ReplicationCode_MC_Table_#.m`) to reproduce the results.
*   Tables and figures will be automatically generated and displayed in the **MATLAB Command Window**.
*   The expected running time is approximately a few minutes for the Monte Carlo simulations and only a few seconds for `ReplicationCode_Table_6_7.m`.

## References

*   Andrews, D. W. (2003). End-of-sample instability tests. *Econometrica* 71, 1661–1694.
*   Diebold, F. X. and R. S. Mariano (1995). Comparing predictive accuracy. *Journal of Business & Economic Statistics* 20, 134–144.
*   Giacomini, R. and B. Rossi (2010). Forecast comparisons in unstable environments. *Journal of Applied Econometrics* 25, 595–620.
*   Harvey, D. I., S. J. Leybourne, R. Sollis, and A. R. Taylor (2021). Real-time detection of regimes of predictability in the US equity premium. *Journal of Applied Econometrics* 36, 45–70.
*   Survey of Professional Forecasters Documentation.

## Authors

*   Fabrizio Iacone
*   Luca Rossini
*   Andrea Viselli
