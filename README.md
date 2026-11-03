# BOAT

**Authors:** Yong Zang and Saijun Zhao.

**Bayesian Organoid-Assisted Treatment evaluation**

BOAT is an R package for randomized umbrella trials that use patient-derived organoid (PDO)
data. It jointly models clinical and PDO responses, adaptively borrows information across
subgroups, and uses each patient's PDO responses to guide treatment assignment.

This repository contains the package source and the accompanying code for the simulation study
and illustrative trial. The illustrative data are simulated; no empirical patient data are included.

## PDF documentation

- [BOAT reference manual](docs/BOAT-manual.pdf): function descriptions, arguments, return values and examples.
- [BOAT user guide](docs/BOAT-guide.pdf): a worked example of trial design, conduct and simulation.

Open either link to view or download the PDF on GitHub. The original R Markdown vignette is in
[vignettes/BOAT.Rmd](vignettes/BOAT.Rmd).

## Installation

Requires **R 4.1 or later** and the **msm** package. The reproduction code was prepared with R 4.5.3.
Download or clone this repository, open R in the repository root, and run:

```r
install.packages("msm")
install.packages("reproduce/BOAT_0.1.0.tar.gz", repos = NULL, type = "source")
library(BOAT)
```

The bundled source archive contains the implementation supplied with the manuscript, with updated
author metadata. To install the source at the repository root instead, run `R CMD INSTALL .`
in a terminal from this folder after installing `msm`.

## Quick start

```r
library(BOAT)

design <- boat_design(
  arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)),
  N = 150,
  interim = c(0.5, 0.75),
  lambda = 0.70
)
boat_boundary(design)
```

See `help(package = "BOAT")` for the function index. After installing the bundled archive,
use `vignette("BOAT")` for the worked example. The PDF links above can be read without installing R.

To rebuild both PDFs, run `Rscript tools/build-docs.R` from the repository root.
This requires `knitr`, `rmarkdown`, Pandoc and a LaTeX installation.

## Functions

| Task | Functions |
| --- | --- |
| Specify a design and its cutoffs | `boat_design()`, `boat_boundary()` |
| Fit the joint model and make decisions | `boat_fit()`, `boat_decide()` |
| Assign a new patient | `boat_rar_model()`, `boat_randomize()` |
| Specify scenarios and generate data | `boat_scenario()`, `boat_dgm()`, `boat_generate_data()` |
| Simulate trials and summarize performance | `boat_simulate_trial()`, `boat_get_oc()` |
| Calibrate the tuning parameter | `boat_calibrate()` |

## Reproduce the simulation study and illustrative trial

The [reproduction guide](reproduce/README.md) maps each script to the manuscript and supplementary
tables. It includes the package source archive, all five designs, calibration,
the main simulation and five sensitivity analyses.

| Analysis | Entry point |
| --- | --- |
| Illustrative trial | [illustration.R](reproduce/illustration/illustration.R) |
| Main simulation | [main.R](reproduce/simulation/main.R) |
| Calibration | [calibrate.R](reproduce/simulation/calibrate.R) |
| Sensitivity analyses | [simulation directory](reproduce/simulation/) |

For example, from a terminal:

```sh
cd reproduce/illustration
Rscript illustration.R
```

The main simulation uses 1,000 replications per setting and scenario. Full simulations and
calibration are computationally intensive; see the reproduction guide for runtime estimates and
parallel execution instructions. The scripts generate numerical results and summary CSV files;
they do not automatically typeset the manuscript or render all publication figures.

For each new set of simulation settings, use a fresh working copy or a separate results folder.
The supplied scripts reuse existing results when rerun.

## Repository contents

```text
DESCRIPTION, NAMESPACE   R package metadata
R/                       package implementation
man/                     R help files
vignettes/               worked example source
docs/                    PDF reference manual and user guide
reproduce/               manuscript reproduction code and package source archive
```

## License

The BOAT package is licensed under **GPL (>= 3)**, as stated in `DESCRIPTION`.
See [LICENSE](LICENSE) for the GNU General Public License, version 3.
