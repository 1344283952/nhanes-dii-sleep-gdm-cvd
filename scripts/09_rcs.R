# ============================================
# 09_rcs.R (003_v2 - W12)
# RCS 限制性立方样条非线性
# DII / sleep_score → 8 outcomes (2 mortality + 6 CVD)
# 用 rms::cph + 加权（同 002 模板；rms 不支持 PSU/Strata，SE 偏乐观）
#
# 输出:
#   output/tables/rcs_pvalues.csv (p-overall + p-nonlinear)
#   output/figures/rcs_*.png (16 张 + overlay)
# ============================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(dplyr)
  library(ggplot2)
})

cat("========================================\n")
cat("W12 RCS 非线性剂量-反应 (003_v2)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)

# M2 协变量 (Wang 2024 同款 lifestyle-only)
m2_covs <- c("age","race","education","pir","bmi",
             "smoke_status","alcohol_any","parity")

# ---- datadist ----
dd_cols <- unique(c("DII", "sleep_score", "WTDR_6YR",
                    m2_covs,
                    "permth", "mort_allcause", "mort_cvd",
                    "cvd_composite", "cvd_chf", "cvd_chd",
                    "cvd_angina", "cvd_mi", "cvd_stroke"))
nhanes_dd <- nhanes_final[, intersect(dd_cols, names(nhanes_final))]
ddist <- datadist(nhanes_dd)
options(datadist = "ddist")

# ---- outcomes & exposures ----
outcomes <- list(
  list(name="all_mortality", type="cox",   var="mort_allcause"),
  list(name="cvd_mortality", type="cox",   var="mort_cvd"),
  list(name="total_cvd",     type="logit", var="cvd_composite"),
  list(name="chf",           type="logit", var="cvd_chf"),
  list(name="chd",           type="logit", var="cvd_chd"),
  list(name="angina",        type="logit", var="cvd_angina"),
  list(name="mi",            type="logit", var="cvd_mi"),
  list(name="stroke",        type="logit", var="cvd_stroke")
)
exposures <- c("DII", "sleep_score")

# ---- 模型 ----
run_rcs_cox <- function(out_var, exp_name) {
  rhs <- paste(c(sprintf("rcs(%s, 4)", exp_name), m2_covs), collapse = " + ")
  f <- as.formula(sprintf("Surv(permth, %s) ~ %s", out_var, rhs))
  tryCatch(cph(f, data = nhanes_final, x = TRUE, y = TRUE,
               weights = nhanes_final$WTDR_6YR),
           error = function(e) { cat("cph ERR:", e$message, "\n"); NULL })
}
run_rcs_logit <- function(out_var, exp_name) {
  rhs <- paste(c(sprintf("rcs(%s, 4)", exp_name), m2_covs), collapse = " + ")
  f <- as.formula(sprintf("%s ~ %s", out_var, rhs))
  tryCatch(lrm(f, data = nhanes_final, x = TRUE, y = TRUE,
               weights = nhanes_final$WTDR_6YR),
           error = function(e) { cat("lrm ERR:", e$message, "\n"); NULL })
}

# ---- P-overall + P-nonlinear（按 rms::anova 输出格式提取）----
extract_p <- function(fit, exposure) {
  if (is.null(fit)) return(c(p_overall = NA, p_nonlin = NA))
  an <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(an)) return(c(p_overall = NA, p_nonlin = NA))
  rn <- rownames(an)
  i_overall <- which(rn == exposure)
  i_nonlin  <- which(grepl("Nonlinear", rn))
  if (length(i_nonlin) > 1 && length(i_overall) > 0) {
    i_nonlin <- i_nonlin[i_nonlin > i_overall][1]
  } else if (length(i_nonlin) >= 1) {
    i_nonlin <- i_nonlin[1]
  } else {
    i_nonlin <- integer(0)
  }
  c(
    p_overall = if (length(i_overall) > 0) an[i_overall, "P"] else NA,
    p_nonlin  = if (length(i_nonlin)  > 0) an[i_nonlin,  "P"] else NA
  )
}

# ---- 绘 RCS 图 ----
plot_rcs <- function(fit, exposure, outcome_label, fname) {
  if (is.null(fit)) return(invisible(NULL))
  pred <- tryCatch(Predict(fit, name = exposure, fun = exp),
                   error = function(e) NULL)
  if (is.null(pred)) return(invisible(NULL))
  df <- as.data.frame(pred)
  names(df)[names(df) == exposure] <- "x"
  p <- ggplot(df, aes(x = x, y = yhat)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#5DADE2", alpha = 0.3) +
    geom_line(color = "#1A5276", linewidth = 1.1) +
    labs(x = exposure, y = "HR / OR (95% CI)",
         title = sprintf("RCS: %s → %s", exposure, outcome_label)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  ggsave(fname, plot = p, width = 5.5, height = 4, dpi = 200)
}

# ---- 主循环 ----
results <- list()
for (out in outcomes) {
  for (exp in exposures) {
    cat(sprintf("--- %s × %s (%s) ---\n", exp, out$name, out$type))
    fit <- if (out$type == "cox") run_rcs_cox(out$var, exp)
           else                   run_rcs_logit(out$var, exp)
    ps  <- extract_p(fit, exp)
    results[[paste(exp, out$name, sep="_")]] <- data.frame(
      exposure       = exp,
      outcome        = out$name,
      outcome_type   = out$type,
      p_overall      = as.numeric(ps["p_overall"]),
      p_nonlinear    = as.numeric(ps["p_nonlin"]),
      stringsAsFactors = FALSE
    )
    if (!is.null(fit)) {
      fname <- sprintf("output/figures/rcs_%s_%s.png", exp, out$name)
      plot_rcs(fit, exp, out$name, fname)
    }
  }
}

res_df <- bind_rows(results)
cat("\n--- RCS p-values ---\n")
print(res_df)

write.csv(res_df, "output/tables/rcs_pvalues.csv", row.names = FALSE)
cat("\n已保存:\n")
cat("  output/tables/rcs_pvalues.csv\n")
cat("  output/figures/rcs_DII_*.png (8 张)\n")
cat("  output/figures/rcs_sleep_score_*.png (8 张)\n")
cat("========================================\n")
