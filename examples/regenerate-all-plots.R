# Regenerate all plots for the three empirical applications and illustration
# Uses cached models (no re-estimation). Skips did event study plots.
# Run from: package/ directory

devtools::load_all()

examples_dir <- file.path(getwd(), "examples")

# ============================================================
# 1. Charoenwong et al. (2019) - Dodd-Frank Act
# ============================================================
message("=== Charoenwong et al. (2019) ===")
rmarkdown::render(
  file.path(examples_dir, "charoenwong-19", "charoenwong-19.Rmd"),
  envir = new.env()
)
message("Charoenwong plots done.")
gc()

# ============================================================
# 2. Hvide & Jones (2018) - Norwegian Patent Reform
# ============================================================
message("\n=== Hvide & Jones (2018) ===")
rmarkdown::render(
  file.path(examples_dir, "hvide-jones-18", "hvide-jones-18.Rmd"),
  envir = new.env()
)
message("Hvide-Jones plots done.")
gc()

# ============================================================
# 3. ADA on Employment
# ============================================================
message("\n=== ADA on Employment ===")
ada_dir <- file.path(examples_dir, "ada-on-employment")
if (!file.exists(file.path(ada_dir, "data", "generated", "ada.RData"))) {
  message("ADA data not found. Running data prep...")
  rmarkdown::render(
    file.path(ada_dir, "ada-on-employment-data-prep.Rmd"),
    envir = new.env()
  )
}
rmarkdown::render(
  file.path(ada_dir, "ada-on-employment-analysis.Rmd"),
  envir = new.env()
)
message("ADA plots done.")
gc()

# ============================================================
# 4. Illustration
# ============================================================
message("\n=== Illustration ===")
source(file.path(examples_dir, "illustration",
                 "plot-illustration-tranformation-invariance.R"),
       chdir = TRUE)
message("Illustration plots done.")

message("\nAll plots regenerated successfully!")
