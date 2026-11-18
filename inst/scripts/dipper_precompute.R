# Precompute the DiPPER results shown in the differential abundance chapter
# (section "DPA with DiPPER", @sec-dpa-dipper).
#
# The book does not run DiPPER, because it needs the CmdStan backend, which is
# installed separately (https://github.com/jepelt/DiPPER#installation).
# Instead, this script is run by hand, and the chapter shows the stored
# results. Rerun it, and update the version and date in the chapter, when the
# data preparation in the chapter or the DiPPER package changes.
#
# Run from the repository root:
#   Rscript inst/scripts/dipper_precompute.R
#
# Outputs:
#   inst/extdata/dipper_tengeler2020.csv       summary(fit), one row per genus
#   inst/pages/figures/dipper_tengeler2020.png  plot(fit), the forest plot

library(mia)
library(DiPPER)
library(ggplot2)

# Prepare the data as in the chapter (chunks load_packages_and_data and
# prepare_data)
data("Tengeler2020", package = "mia")
tse <- Tengeler2020
tse <- subsetByPrevalent(tse, rank = "Genus", prevalence = 0.10)
tse$patient_status <- factor(tse$patient_status, levels = c("Control", "ADHD"))
tse$cohort <- factor(tse$cohort, levels = c("Cohort_1", "Cohort_2", "Cohort_3"))

# Run DiPPER as in the chapter (chunk run_dpa_dipper)
set.seed(123)
fit_dipper <- dipper(
    tse = tse,
    assay.type = "counts",
    data.type = "counts",
    formula = ~ patient_status + cohort,
    read.depth = TRUE,
    seed = 123
)

res_dipper <- summary(fit_dipper)
write.csv(res_dipper, "inst/extdata/dipper_tengeler2020.csv", row.names = FALSE)

p <- plot(fit_dipper)
ggsave("inst/pages/figures/dipper_tengeler2020.png", p,
    width = 7, height = 4, dpi = 150, bg = "white"
)

cat("DiPPER", as.character(packageVersion("DiPPER")), "on", format(Sys.Date()), "\n")
print(fit_dipper)
