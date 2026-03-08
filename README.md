# ak

R package for event-study designs with discrete outcomes using transition independence. The package estimates treatment effects when outcomes are discrete and categorical, replacing the parallel trends assumption with transition independence — the assumption that, absent treatment, transition dynamics conditional on pre-treatment outcomes are identical across groups. Accommodates unobserved heterogeneity via latent group structure and provides flow decomposition of treatment effects.

## Installation

```r
# From the package/ directory:
devtools::install()

# Or load for development without installing:
devtools::load_all()
```

## Replication

The package includes three empirical applications and one pedagogical illustration. All examples use cached estimated models by default, so plots can be regenerated without re-estimation.

### Quick start (all examples)

From the `package/` directory in R:

```r
devtools::load_all()
source("examples/regenerate-all-plots.R")
```

### Hvide & Jones (2018) — Norwegian Patent Reform

Effect of a university patent rights reform on researcher patenting rates in Norway (1995--2010). Binary outcome: whether an inventor filed at least one patent in a year. Compares university researchers (treated) to non-university inventors (control).

**Data source:** AER replication files (`data/raw/113125-V1.zip`).

```r
devtools::load_all()
rmarkdown::render("examples/hvide-jones-18/hvide-jones-18.Rmd", envir = new.env())
```

### ADA on Employment (1990)

Employment effects of the Americans with Disabilities Act on individuals with disabilities. Three-state outcome: employed, unemployed, or out of the labor force. Data from the SIPP 1990 Panel, comparing individuals with work-limiting disabilities (treated) to those without (control).

**Data source:** SIPP 1990 Panel. Run the data preparation notebook first if `data/generated/ada.RData` does not exist.

```r
devtools::load_all()

# If data/generated/ada.RData does not exist, run data prep first:
rmarkdown::render("examples/ada-on-employment/ada-on-employment-data-prep.Rmd", envir = new.env())

# Then run the analysis:
rmarkdown::render("examples/ada-on-employment/ada-on-employment-analysis.Rmd", envir = new.env())
```

### Charoenwong et al. (2019) — Dodd-Frank Act

Effect of the Dodd-Frank regulatory jurisdiction shift on investment adviser complaint rates (2006--2016). Binary outcome: whether a registered investment adviser received a customer complaint in a year. The 2012 reform transferred midsize advisers from SEC to state oversight.

**Data source:** AER replication files (`data/raw/116210-V1.zip`).

```r
devtools::load_all()
rmarkdown::render("examples/charoenwong-19/charoenwong-19.Rmd", envir = new.env())
```

### Illustration — Parallel Trends vs Transition Independence

Simulation-based pedagogical example comparing the parallel trends assumption to the transition independence assumption.

```r
devtools::load_all()
source("examples/illustration/plot-illustration-tranformation-invariance.R")
```

## Re-estimation

Cached models are stored in `models/` directories within each example folder. To re-estimate from scratch, delete the relevant `models/` directory and rerun the example. Note that re-estimation with bootstrap (500 replications) is computationally intensive.

## Running Tests

```r
devtools::test()
```
