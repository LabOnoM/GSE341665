#!/usr/bin/env Rscript
# 08_mendelian_randomization.R  --  OPTIONAL human-genetics support (template)
#
# Purpose
# -------
# Two-sample Mendelian randomization (MR) testing whether genetically-proxied
# cholesterol lowers/raises risk of DENTAL disease outcomes. This provides
# human-genetic PLAUSIBILITY for the cholesterol-metabolism angle discovered in
# the Smad4 dcKO ameloblast transcriptome.
#
# IMPORTANT SCOPE CAVEATS (state these in the manuscript if used):
#   1. Exposure mismatch: MR on SERUM cholesterol is not the same as intracellular
#      ameloblast STEROL BIOSYNTHESIS. A cis-HMGCR instrument (statin proxy) is the
#      closest available proxy for perturbing the biosynthetic pathway.
#   2. Outcome mismatch: no GWAS of amelogenesis imperfecta / enamel formation
#      exists (AI is monogenic/rare). Available dental outcomes are caries,
#      periodontitis, and tooth loss/number of teeth -- diseases of MATURE teeth.
#      A positive result is supportive plausibility; a null does NOT refute the
#      developmental, tissue-local hypothesis.
#
# PREREQUISITES (this script does NOT run in the analysis session -- it needs a token):
#   install.packages("remotes")
#   remotes::install_github(c("MRCIEU/TwoSampleMR","MRCIEU/ieugwasr","rondolab/MR-PRESSO"))
#   Register a free OpenGWAS token at https://api.opengwas.io and paste below.
#
# USAGE: edit the TOKEN and the *_ID vectors (verify current IDs with
#   ieugwasr::gwasinfo() / TwoSampleMR::available_outcomes()), then:
#   Rscript scripts/JDR_reanalysis/08_mendelian_randomization.R

# --- one-time bootstrap: install packages if missing ---
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github(c("MRCIEU/ieugwasr", "MRCIEU/TwoSampleMR"))
}
suppressMessages({library(TwoSampleMR)})

# --- OpenGWAS token (register free at https://api.opengwas.io) ---
# Preferred: pass it via the environment so it is NEVER written to disk:
#   OPENGWAS_JWT="your_token" Rscript scripts/JDR_reanalysis/08_mendelian_randomization.R
# Fallback: replace the placeholder below. Never commit a real token; regenerate if exposed.
if (Sys.getenv("OPENGWAS_JWT") == "")
  Sys.setenv(OPENGWAS_JWT = "PASTE_YOUR_OPENGWAS_TOKEN_HERE")
OUT <- "scripts/JDR_reanalysis/output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Exposures (VERIFIED OpenGWAS IDs; confirm current with ieugwasr::gwasinfo) ----
# Systemic lipid exposures:
EXPOSURES <- c(
  LDL_cholesterol_GLGC = "ieu-a-300",   # GLGC 2013 (Willer), 173,082 samples
  LDL_cholesterol_UKB  = "ieu-b-110",   # larger LDL GWAS (replication of the LDL arm)
  total_cholesterol    = "ieu-a-301",   # GLGC total cholesterol
  HDL_cholesterol      = "ieu-a-299"    # negative-control-ish lipid arm
)
# PREFERRED mechanistic instrument: cis-HMGCR "statin proxy" = LDL-associated SNPs
# within +/-100 kb of HMGCR (chr5:74,632,154-74,657,929, GRCh37), LD-clumped.
# Build with: extract_instruments("ieu-a-300") then subset to the HMGCR window.
HMGCR_WINDOW <- list(chr = 5, start = 74532154, end = 74757929)

## ---- 2. Outcomes: dental disease ----
# Route A (OpenGWAS, FinnGen R5 'finn-b-' prefix). VERIFY exact IDs with
# ieugwasr::gwasinfo() / https://opengwas.io -- FinnGen endpoint codes change by release.
OUTCOMES <- c(
  periodontitis_chronic = "finn-b-K11_PERIODON_CHRONIC", # chronic periodontitis (FinnGen)
  dental_problems_UKB   = "ukb-b-19191"                  # UKB "Mouth/teeth dental problems: Toothache"
  # dental_caries       = "<verify finn-b caries code>"  # confirm FinnGen caries endpoint id
)
# Route B (recommended for CARIES/PERIODONTITIS): the GLIDE consortium summary stats
# (Shungin et al. 2019, Nat Commun; PMID 31235808) are downloadable from the Univ. of
# Bristol data repository. Load locally and use instead of / alongside Route A:
#   out <- read_outcome_data("GLIDE_dental_caries.txt", snps = exp$SNP,
#            snp_col="SNP", beta_col="beta", se_col="se", effect_allele_col="ea",
#            other_allele_col="oa", eaf_col="eaf", pval_col="pval")
GLIDE_LOCAL <- NULL   # set to a file path to enable Route B

## ---- 3. Run MR per exposure x outcome with full sensitivity analyses ----
run_pair <- function(exp_id, out_id, exp_name, out_name) {
  message(sprintf("MR: %s -> %s", exp_name, out_name))
  exp <- clump_data(extract_instruments(exp_id))               # genome-wide instruments, LD-clumped
  out <- extract_outcome_data(exp$SNP, out_id)
  dat <- harmonise_data(exp, out)
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  het  <- mr_heterogeneity(dat)                                # Cochran's Q
  pleio<- mr_pleiotropy_test(dat)                              # MR-Egger intercept
  list(exposure = exp_name, outcome = out_name, res = res, het = het, pleio = pleio,
       n_snp = nrow(dat))
}

results <- list()
for (en in names(EXPOSURES)) for (on in names(OUTCOMES)) {
  if (grepl("^FILL_", OUTCOMES[[on]])) next                    # skip un-filled outcome IDs
  results[[paste(en, on, sep = "__")]] <- tryCatch(
    run_pair(EXPOSURES[[en]], OUTCOMES[[on]], en, on),
    error = function(e) { message("  skipped: ", conditionMessage(e)); NULL })
}

## ---- 4. Report ----
tab <- do.call(rbind, lapply(results, function(r) if (!is.null(r))
  cbind(exposure = r$exposure, outcome = r$outcome, n_snp = r$n_snp, r$res)))
if (!is.null(tab)) {
  write.csv(tab, file.path(OUT, "MR_results.csv"), row.names = FALSE)
  print(tab)
  cat("\nInterpretation guide: concordant IVW + weighted-median effects with a\n",
      "non-significant MR-Egger intercept (no directional pleiotropy) support a\n",
      "causal cholesterol->dental-disease relationship. Report OR/beta, 95% CI, and\n",
      "the sensitivity analyses; frame as human-genetic plausibility, not proof of\n",
      "the ameloblast sterol-biosynthesis mechanism.\n")
} else {
  cat("No results -- fill in valid outcome IDs (OUTCOMES) and a token first.\n")
}
