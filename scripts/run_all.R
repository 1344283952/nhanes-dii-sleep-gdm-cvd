# ============================================
# run_all.R (003_v2 — full pipeline)
# 用法: 在 projects/003_v2_dii_sleep_gdm/ 目录下执行
#   cd projects/003_v2_dii_sleep_gdm
#   Rscript scripts/run_all.R
# 所有相对路径 (scripts/, data/, output/) 以当前 cwd 为基准
# ============================================

cat("╔══════════════════════════════════════════╗\n")
cat("║  003_v2 DII × Sleep × GDM — 全流程       ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

start_time <- Sys.time()

# Sequenced for the 003_v2 pipeline. Comment out 01 if raw data already in
# data/raw/ to skip re-download (~5-20 min depending on bandwidth).
scripts <- c(
  "scripts/01_download_data.R",     # raw NHANES + mortality download
  "scripts/02_merge_data.R",        # 6 cycles + linked mortality
  "scripts/03_clean_data.R",        # variable construction + analytic sample
  "scripts/04_survey_design.R",     # svydesign with WTDR_6YR / SDMVPSU / SDMVSTRA
  "scripts/05_table1.R",            # Table 1 baseline characteristics
  "scripts/06_cox_main.R",          # Table 2 Cox + Schoenfeld PH tests
  "scripts/07_logistic_cvd.R",      # Table 3 logistic + BH FDR
  "scripts/08_subgroup.R",          # Table 4 subgroup × interactions
  "scripts/09_rcs.R",               # RCS dose-response (Supplementary Fig S4/S5)
  "scripts/11_mediation.R",         # Table 5 CMAverse mediation
  "scripts/12_sensitivity.R",       # Table 6 S1-S6 (incl. mice m=20, IPTW)
  "scripts/13_predictive.R",        # Table 7 C-index + NRI/IDI + E-values
  "scripts/14_consort.R",           # Fig 1 CONSORT
  "scripts/15_dag.R",               # Fig 2 DAG
  "scripts/16_forest.R"             # Fig 3A-D subgroup forest plots
)

n_ok <- 0
n_fail <- 0
fail_list <- character()
for (s in scripts) {
  cat(paste0("\n>>> 正在执行: ", s, "\n"))
  cat(paste0(rep("-", 50), collapse = ""), "\n")

  tryCatch({
    source(s, echo = FALSE)
    cat(paste0(">>> ", s, " 执行成功 ✓\n"))
    n_ok <- n_ok + 1
  }, error = function(e) {
    cat(paste0(">>> ", s, " 执行失败 ✗\n"))
    cat(paste0("    错误: ", e$message, "\n"))
    n_fail <<- n_fail + 1
    fail_list <<- c(fail_list, s)
  })
}

end_time <- Sys.time()
elapsed <- round(difftime(end_time, start_time, units = "mins"), 1)

cat(paste0("\n╔══════════════════════════════════════════╗\n"))
cat(sprintf("║  全部完成: 成功 %d / 失败 %d / 耗时 %s 分钟\n", n_ok, n_fail, elapsed))
cat(paste0("╚══════════════════════════════════════════╝\n"))

if (length(fail_list) > 0) {
  cat("\n失败脚本:\n")
  for (f in fail_list) cat(paste0("  ", f, "\n"))
}

cat("\n关键输出 (post-pipeline):\n")
cat("  output/tables/table1.xlsx                — baseline characteristics\n")
cat("  output/tables/table2_cox.xlsx            — Cox main + PH + GDM × DII\n")
cat("  output/tables/table_schoenfeld_ph.csv    — Schoenfeld per-variable + global\n")
cat("  output/tables/table3_logistic_cvd.xlsx   — logistic + BH FDR\n")
cat("  output/tables/table3c_logistic_fdr.csv   — FDR q-values\n")
cat("  output/tables/table4_dii_gdm_interaction.csv\n")
cat("  output/tables/table5_mediation.xlsx      — CMAverse mediation\n")
cat("  output/tables/table6_subgroup.xlsx       — subgroup-stratified\n")
cat("  output/tables/table7_sensitivity.xlsx    — S1-S6 (incl. mice m=20)\n")
cat("  output/tables/table7b_mice_per_imp.csv   — per-imputation Cox\n")
cat("  output/tables/table8_predictive.xlsx     — C-index + NRI/IDI + E-values\n")
cat("  output/figures/fig1_consort.png\n")
cat("  output/figures/fig2_dag.png\n")
cat("  output/figures/fig3{a-d}_forest_*.png\n")
