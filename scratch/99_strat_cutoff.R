# =====================================================================
# find_optimal_ml_thresholds()
# =====================================================================
#' Find operationally- and statistically-feasible ML risk thresholds
#'
#' Given a precision-recall curve and the historical distribution of model
#' scores, this function identifies two "safe zone" thresholds that bound
#' the proportion of the clinic population that should be flagged by the
#' model:
#'
#'   - a FLOOR threshold: the minimum proportion needed to reach your
#'     target trial sample size (n_req) within the trial window (w weeks)
#'   - a CEILING threshold: the maximum proportion the care team can
#'     operationally handle, given call capacity (c_max) per week
#'
#' Any threshold you choose between these two values should simultaneously
#' satisfy statistical power requirements and stay within staff capacity.
#'
#' @param pr_curve_data A tibble/data frame from yardstick::pr_curve(),
#'   containing (at minimum) the columns `.threshold` and `precision`.
#'   Typically produced via:
#'   collect_predictions(...) %>% pr_curve(truth = ..., .pred_DNA)
#'
#' @param historical_scores Numeric vector of historical ML risk scores
#'   (predicted probabilities) from your cohort. Used to empirically
#'   estimate, for each candidate threshold, what proportion of real
#'   patients would be flagged. Must be a plain numeric vector — if
#'   extracting from collect_predictions(), remember to pull() the score
#'   column first, e.g. collect_predictions(...) %>% pull(.pred_DNA).
#'
#' @param v_week Numeric. Expected weekly clinic volume (patients/week)
#'   eligible for enrollment. Default 308.
#'
#' @param w Numeric. Planned trial duration in weeks. Default 26.
#'
#' @param n_req Numeric. Target total sample size required for adequate
#'   statistical power. Default 800.
#'
#' @param a_int Numeric in (0, 1]. Allocation ratio to the intervention
#'   arm (e.g. 0.50 = 50% of enrolled patients receive the intervention).
#'   Default 0.50.
#'
#' @param r_call Numeric in (0, 1]. Expected non-response rate among
#'   intervention-arm patients that will require a follow-up phone call
#'   from a coordinator. Default 0.70.
#'
#' @param c_max Numeric. Maximum number of calls coordinators can
#'   realistically handle per week (operational capacity constraint).
#'   Default 30.
#'
#' @return An (invisible) list with:
#'   \describe{
#'     \item{floor_threshold}{Model score threshold (>=) that most closely
#'       achieves the statistical floor proportion.}
#'     \item{ceiling_threshold}{Model score threshold (>=) that most
#'       closely achieves the operational ceiling proportion.}
#'     \item{floor_prop}{Empirical proportion of historical_scores flagged
#'       at floor_threshold.}
#'     \item{ceiling_prop}{Empirical proportion of historical_scores
#'       flagged at ceiling_threshold.}
#'     \item{floor_target_prop}{The theoretical target proportion derived
#'       from n_req, v_week, and w (before matching to an actual
#'       threshold).}
#'     \item{ceiling_target_prop}{The theoretical target proportion
#'       derived from c_max, v_week, a_int, and r_call (before matching).}
#'     \item{feasible}{Logical. TRUE if floor_prop <= ceiling_prop, i.e.
#'       a valid "safe zone" of thresholds exists between the two bounds.}
#'     \item{enriched_pr_curve}{The full pr_curve_data with an added
#'       prop_flagged column, useful for downstream plotting or
#'       sensitivity analysis.}
#'   }
#'
#' @examples
#' \dontrun{
#' optimal_thresholds <- find_optimal_ml_thresholds(
#'   pr_curve_data     = collect_predictions(model_tune$tune_res,
#'                          parameters = model_tune$best_params) %>%
#'                        pr_curve(truth = dna_outcome, .pred_DNA),
#'   historical_scores = collect_predictions(model_tune$tune_res,
#'                          parameters = model_tune$best_params) %>%
#'                        pull(.pred_DNA),
#'   v_week = 308,
#'   w      = 26,
#'   n_req  = 800
#' )
#' }
find_optimal_ml_thresholds <- function(
    pr_curve_data,
    historical_scores,
    v_week = 308,
    w = 26,
    n_req = 800,
    a_int = 0.50,
    r_call = 0.70,
    c_max = 30
) {
  
  # ---------------------------------------------------------
  # 1. Calculate Operational & Statistical Bounds (Proportions)
  # ---------------------------------------------------------
  floor_prop   <- n_req / (v_week * w)
  ceiling_prop <- c_max / (v_week * a_int * r_call)
  is_feasible  <- floor_prop <= ceiling_prop
  
  cat("====================================================\n")
  cat("      ML THRESHOLD SAFE ZONE EXTRACTOR              \n")
  cat("====================================================\n")
  cat(sprintf("Statistical Floor Target:   %.2f%% of clinic (Min N = %d over %d weeks)\n",
              floor_prop * 100, n_req, w))
  cat(sprintf("Operational Ceiling Target: %.2f%% of clinic (Max Calls = %d/week)\n",
              ceiling_prop * 100, c_max))
  cat("----------------------------------------------------\n")
  
  if (!is_feasible) {
    cat("[!] WARNING: Statistical floor exceeds operational ceiling!\n")
    cat("    No single threshold can satisfy both constraints -\n")
    cat("    consider raising c_max, extending w, or lowering n_req.\n\n")
  }
  
  # ---------------------------------------------------------
  # 2. Enrich PR Curve with Empirical Proportion Flagged (VECTORIZED)
  # ---------------------------------------------------------
  sorted_scores <- sort(historical_scores)
  n_hist <- length(sorted_scores)
  
  enriched_pr <- pr_curve_data %>%
    filter(.threshold != Inf & !is.na(.threshold)) %>%
    mutate(
      n_less       = findInterval(.threshold, sorted_scores, left.open = TRUE),
      prop_flagged = (n_hist - n_less) / n_hist
    ) %>%
    select(-n_less)
  
  # ---------------------------------------------------------
  # 3. Match Target Proportions to Closest PR Curve Thresholds
  # ---------------------------------------------------------
  floor_match <- enriched_pr %>%
    slice_min(abs(prop_flagged - floor_prop), n = 1, with_ties = FALSE)
  
  ceiling_match <- enriched_pr %>%
    slice_min(abs(prop_flagged - ceiling_prop), n = 1, with_ties = FALSE)
  
  # Estimated weekly/total volumes at each matched threshold, for context
  floor_n_per_week   <- floor_match$prop_flagged * v_week
  ceiling_n_per_week <- ceiling_match$prop_flagged * v_week
  floor_total_n      <- floor_n_per_week * w
  ceiling_est_calls  <- ceiling_n_per_week * a_int * r_call
  
  # ---------------------------------------------------------
  # 4. Console Report
  # ---------------------------------------------------------
  cat("RECOMMENDED MODEL THRESHOLD RANGE:\n\n")
  cat(sprintf("-> Power Floor Threshold:      >= %.4f\n", floor_match$.threshold))
  cat(sprintf("   Flags top %.1f%% of clinic (~%.0f pts/week, ~%.0f over %d weeks) | Precision: %.3f\n",
              floor_match$prop_flagged * 100, floor_n_per_week, floor_total_n, w, floor_match$precision))
  cat(sprintf("-> Capacity Ceiling Threshold:  >= %.4f\n", ceiling_match$.threshold))
  cat(sprintf("   Flags top %.1f%% of clinic (~%.0f pts/week, ~%.0f est. calls/week) | Precision: %.3f\n",
              ceiling_match$prop_flagged * 100, ceiling_n_per_week, ceiling_est_calls, ceiling_match$precision))
  cat("----------------------------------------------------\n")
  cat(sprintf("Safe zone: %s\n", if (is_feasible) "VALID (floor <= ceiling)" else "INVALID (floor > ceiling)"))
  cat("====================================================\n\n")
  
  # ---------------------------------------------------------
  # 5. Return structured list for downstream pipeline scripts
  # ---------------------------------------------------------
  return(invisible(list(
    floor_threshold     = floor_match$.threshold,
    ceiling_threshold   = ceiling_match$.threshold,
    floor_prop          = floor_match$prop_flagged,
    ceiling_prop        = ceiling_match$prop_flagged,
    floor_target_prop   = floor_prop,
    ceiling_target_prop = ceiling_prop,
    feasible            = is_feasible,
    enriched_pr_curve   = enriched_pr
  )))
}

optimal_thresholds <- find_optimal_ml_thresholds(
  pr_curve_data     = collect_predictions(model_tune$tune_res, parameters = model_tune$best_params) %>%
    pr_curve(truth = dna_outcome, .pred_DNA),
  historical_scores = collect_predictions(model_tune$tune_res, parameters = model_tune$best_params) %>%
    pull(.pred_DNA),
  v_week            = 308,
  w                 = 26,
  n_req             = 800
)