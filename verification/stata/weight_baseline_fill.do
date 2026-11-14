* Mirrors the authors' Table 3 steps for weight at 6 months (01_published_tables.do lines 318-403)
* and shows which baseline weight the fill step (line 394) puts into the model for the
* 5 patients with no observed baseline. Usage: stata-se -b do weight_baseline_fill.do [sortseed]
args sortseed
insheet using "analysis_sample.csv", comma clear
if "`sortseed'" != "" set sortseed `sortseed'
local var weight
local month 6
drop if mi(`var'_val_`month'mo)
unab todrop: `var'_*mo
unab tokeep: `var'_*_`month'mo `var'_*_0mo
local todrop: list todrop - tokeep
drop `todrop'
gen `var'_baseline = `var'_val_0mo
reshape long `var'_val_@mo `var'_date_@mo, i(subject_id) j(months)
rename `var'_val*mo `var'_panel
rename `var'_date*mo `var'_time
replace `var'_time = `var'_time - consent_date
gen post_period = (months == `month')
gen post_treat_now = post_period * treat_now
egen n_subject_id = group(subject_id)
xtset n_subject_id post_period
keep `var'_panel post_period treat_now post_treat_now rural `var'_baseline `var'_time hba1c_cat_control age5_* male race_nonhisp_black hispanic n_subject_id
mi set mlong
mi register imputed `var'_panel `var'_time
mi impute mvn `var'_panel `var'_time = treat_now post_period post_treat_now rural hba1c_cat_control male race_nonhisp_black hispanic age5_*, add(10) rseed(9112023)
gen byte miss_base = mi(`var'_baseline)
bysort n_subject_id: egen anymiss = max(miss_base)
sort n_subject_id post_period
bysort n_subject_id: replace `var'_baseline = `var'_panel[1] if mi(`var'_baseline)
gen post_baseline = post_period * `var'_baseline
* compact summary per patient
keep if anymiss
bysort n_subject_id: egen fu_weight = max(cond(post_period==1, weight_panel, .))
bysort n_subject_id: egen imp_min = min(cond(post_period==0 & _mi_m>0, weight_panel, .))
bysort n_subject_id: egen imp_max = max(cond(post_period==0 & _mi_m>0, weight_panel, .))
gen used_in_model = post_baseline if post_period==1
bysort n_subject_id: egen used = max(used_in_model)
gen fu_row_dropped = mi(used)
bysort n_subject_id: keep if _n==1
list n_subject_id treat_now fu_weight imp_min imp_max used fu_row_dropped, noobs clean abbreviate(14)
