* The authors' Table 3 weight model fitted on observed baselines only (no imputation).
* kg45 = the posted weight_val column (lb x 0.45); kg = weight_lb_val / 2.20462.
foreach month in 6 12 {
  foreach unit in kg45 kg {
    insheet using "analysis_sample.csv", comma clear
    unab agebins: age5_*
    unab ageomit: age5_10
    local agebins: list agebins - ageomit
    if "`unit'" == "kg" {
        foreach m in 0 6 12 {
            replace weight_val_`m'mo = weight_lb_val_`m'mo / 2.20462
        }
    }
    local var weight
    drop if mi(`var'_val_`month'mo) | mi(`var'_val_0mo)
    keep subject_id treat_now rural hba1c_cat_control age5_* male race_nonhisp_black hispanic consent_date `var'_val_0mo `var'_val_`month'mo `var'_date_0mo `var'_date_`month'mo
    gen `var'_baseline = `var'_val_0mo
    reshape long `var'_val_@mo `var'_date_@mo, i(subject_id) j(months)
    rename `var'_val*mo `var'_panel
    rename `var'_date*mo `var'_time
    replace `var'_time = `var'_time - consent_date
    gen post_period = (months == `month')
    gen post_treat_now = post_period * treat_now
    gen post_baseline = post_period * `var'_baseline
    egen n_subject_id = group(subject_id)
    xtset n_subject_id post_period
    quietly xtreg `var'_panel treat_now post_period post_treat_now post_baseline `var'_time rural hba1c_cat_control `agebins' male race_nonhisp_black hispanic, fe i(n_subject_id) robust
    di as result "RESULT month=`month' unit=`unit' N_pat=" e(N_g) " b=" %6.3f _b[post_treat_now] " lo=" %6.3f r(table)["ll","post_treat_now"] " hi=" %6.3f r(table)["ul","post_treat_now"] " p=" %6.4f r(table)["pvalue","post_treat_now"]
  }
}
