*********************************
* NELSON & PLOSSER REPLICATION
*********************************

	use "NelsonPlosserData.dta", clear
	
	set more off
	summarize

	* The variables are presented in their raw form, and then once again in their logarithms.
	* The logged variabes are denoted with an ``l'' prefix.

	keep l* year bnd


foreach v of var bnd-lsp500 {
	dfuller `v'
}
*

* Generate the time variables (each variable gets it's own time,
	* because the numbering starts with "1" for each variable in
	* Nelson and Plosser's paper
	foreach v of var bnd-lsp500{
		gen time_`v' = .
		replace time_`v' = year if `v'~=.
		quietly summarize time_`v'
		*replace time_`v' = time_`v' - r(min) + 1 /* if time started at year=1 */
		replace  time_`v' = time_`v' - r(min)     /* How N&P do it (time=0)    */
		* browse year `v'  time_`v'
	}






********************************
* N&P's Table 2
********************************

* Replicating the results
foreach v of var bnd-lsp500 {
	display ""
	display "-----------------------------------------"
	display "`v'"
	corrgram `v', noplot lags(6)
}








*****************************************
* N&P's Table 3: Sample autocorrelations 
* of the first differences of the logs
*****************************************

foreach v of var bnd-lsp500 {
	display ""
	display "-----------------------------------------"
	display "D.`v'"
	corrgram D.`v', noplot lags(6)
}






*****************************************
* N&P's Table 4: Sample autocorrelations 
* of the deviations from time trend
*****************************************

foreach v of var bnd-lsp500 {
	display ""
	display "-----------------------------------------"
	display "Deviations from trend of: `v'"
	quietly reg `v' year
	quietly predict `v'devs, resid
	corrgram `v'devs, noplot lags(6)
	drop `v'devs
}
*







*****************************************
* N&P's Table 5: Tests for Autoregressive 
* Unit Roots
*****************************************


*use "NelsonPlosserData.dta", clear


* Real GNP
local var = "lrgnp"
dfuller `var', trend reg lags(1)
local lags = r(lags)
reg D.`var' L1.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time]/_se[time]
display "rho_1 (coeff on lagged term) = " 1 + _b[L.l]
display "tau(rho_1) = " _b[L.l]/_se[L.l]
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
correlate errors l.errors
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Nominal GNP
local var = "lgnp"
dfuller `var' , trend reg lags(1)
reg D.`var' L1.`var' LD.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time]/_se[time]
display "rho_1 (coeff on lagged term) = " 1 + _b[L.l]
display "tau(rho_1) = " _b[L.l]/_se[L.l]
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Real per capita GNP
local var = "lpcrgnp"
dfuller `var'  , trend reg lags(1)
local lags = r(lags)
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time]/_se[time]
display "rho_1 (coeff on lagged term) = " 1 + _b[L.l]
display "tau(rho_1) = " _b[L.l]/_se[L.l]
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Industrial Production
local var = "lip"
dfuller `var'  , trend reg lags(5)
local lags = r(lags)
reg D.`var' L1.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time]/_se[time]
display "rho_1 (coeff on lagged term) = " 1 + _b[L.l]
display "tau(rho_1) = " _b[L.l]/_se[L.l]
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors




* Employment
local var = "lemp"
dfuller `var'  , trend reg lags(2)
local lags = r(lags)
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_lemp]/_se[time]
display "rho_1 (coeff on lagged term) = " 1 + _b[L.l]
display "tau(rho_1) = " _b[L.l]/_se[L.l]
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
* The residual autocorrelation on employment is wrong in N&P
display "r1 = " r(ac1)
drop errors




* Unemployment
local var  "lun"
local lags "3"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
local r1 = r(ac1)
display "r1 = " `r1'
drop errors



* GNP Deflator
local var  "lprgnp"
local lags "1"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* CPI
local var  "lcpi"
local lags "3"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
* The correlation coeff is a little off on this one
display "r1 = " r(ac1)
drop errors



* Wages
local var  "lwg"
local lags "2"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Real Wages
local var  "lrwg"
local lags "1"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Money stock
local var  "lm"
local lags "1"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* Velocity 
local var  "lvel"
local lags "3"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors
 



* Interest rate
local var  "bnd"
local lags "2"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors



* S&P 500 
local var  "lsp500"
local lags "2"
dfuller `var'  , trend reg lags(`lags')
reg D.`var' L.`var' L(1/`lags')D.`var' time_`var'
local k = `lags'+1
display "K = `k'"
display "mu (the constant) = " _b[_cons]
display "t(mu) = " _b[_cons]/_se[_cons]
display "gamma (the coeff on time) = " _b[time]
display "t(gamma) = " _b[time_`var']/_se[time_`var']
display "rho_1 (coeff on lagged term) = " 1 + _b[L.`var']
display "tau(rho_1) = " _b[L.`var']/_se[L.`var']
display "s(u) the std err of the regression  = " e(rmse)
predict errors, resid
corrgram errors, lags(1)
display "r1 = " r(ac1)
drop errors














*********************************
* DF-GLS (ERS)
*********************************


use "NelsonPlosserData.dta", clear
keep l* year bnd
order lrgnp lgnp lpcrgnp lip lemp lun lprgnp lcpi lwg lrwg lm lvel bnd lsp500	

set more off


* DFGLS-ERS
	foreach v of var lrgnp-lsp500 {
	display ""
	display ""
	display ""
	display ""
	display ""
	display "-----------------------------------------"
	display "`v'"
	dfgls `v' , ers
	}
* In NO cases do we reject the Ho of a Unit root ==> we have unit roots
* This is true for all lags (up to the max that we tested, lag=10)











