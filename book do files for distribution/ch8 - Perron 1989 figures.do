

set scheme s1mono

************************************************************************
* PERRON (1989) replication
* "The Great Crash, the Oil Price Shock, and the Unit Root Hypothesis"
************************************************************************


********************************
* Getting the data in order
********************************

cd "C:\Users\johnlevendis\Desktop\book do files\Perron\"
use "NelsonPlosserData.dta", clear

set more off
keep l* year bnd

foreach v of var bnd-lsp500 {
	dfuller `v'
}


* Generate the time variables (each variable gets it's own time,
	* because the numbering starts with "1" for each variable in
	* Nelson and Plosser's paper
	foreach v of var bnd-lsp500{
		gen time_`v' = .
		replace time_`v' = year if `v'~=.
		quietly summarize time_`v'
		*replace time_`v' = time_`v' - r(min) + 1 /* if time started at year=1 */
		replace  time_`v' = time_`v' - r(min)     /* How N&P do it (time=0)    */
	}




* Dummy variables indicating the one period after the shock
	gen DGD = 0
	replace DGD = 1 if year == 1930
	*gen DOS = 0
	*replace DOS = 1 if year == 1974







********************************
* Figure 8-5
********************************

* Figure 8.5.a
	gen DU = 0
	replace DU = 1 if year > 1929
	quietly reg lwg DU year
	predict yhat
	graph twoway (connected lwg year if (year>=1900 & year <=1970), msymbol(i)) (connected yhat year if (year>=1900 & year <=1970), msymbol(i)) 
	drop yhat


* Figure 8.5.b
	gen DT = 0
	replace DT = year if year > 1929
	quietly reg lsp500 DU year DT
	predict yhat
	graph twoway (connected lsp500 year if (year>=1870 & year <=1970), msymbol(i)) (connected yhat year if (year>=1870 & year <=1970), msymbol(i))
	drop yhat









**************************************************
* Figure 8.1
* Structural Breaks and Unit Roots
**************************************************

	drop _all
	clear all
	set obs 200
	set seed 1234
	gen time = _n
	tsset time
	gen error = rnormal()
	gen X = error in 1
	replace X = 0.50 * L.X + error in 2/L
	replace X = 20+X in 101/200
	graph twoway (connected X time, msymbol(i) xline(100))
	graph save sb1, replace
	graph export sb1.pdf, replace


* Same trend, different intercepts
	drop _all
	clear all
	set obs 200
	set seed 1234
	gen time = _n
	tsset time
	gen DU = 0
	replace DU = 1 if time > 100
	gen error = rnormal()
	local beta1 = "0.10"
	local mu1 = "1"
	local mu2 = "10" 
	gen X = error + `mu1' in 1
	replace X = `mu1'+ `beta1'*time + (`mu2' - `mu1')*DU + error in 2/L
	graph twoway (connected X time, msymbol(i) xline(100)) 

	

* Two trend-stationary processes, with different trends (slopes)
	drop _all
	clear all
	set obs 200
	set seed 1234
	gen time = _n
	tsset time
	gen TimeSince = 0
	replace TimeSince = time - 100 if time > 100
	gen error = rnormal()
	local beta1 = "0.10"
	local beta2 = "0.60"
	local mu1 = "0"
	gen X = error + `mu1' in 1
	replace X = `mu1'+ `beta1'*time + (`beta2' - `beta1')*TimeSince + error in 2/L
	graph twoway (connected X time, msymbol(i) xline(100)) 
	

	
* Two trend-stationary processes, with different interepts and trends (slopes)
	drop _all
	clear all
	set obs 200
	set seed 1234
	gen time = _n
	tsset time
	gen DT = 0
	replace DT = time if time > 100
	gen DU = 0
	replace DU = 1 if time > 100
	gen error = rnormal()
	local beta1 = "0.10"
	local beta2 = "0.60"
	local mu1 = "1"
	local mu2 = "10" 
	gen X = error + `mu1' in 1
	replace X = `mu1'+ `beta1'*time + (`mu2' - `mu1')*DU + (`beta2' - `beta1')*DT + error in 2/L
	graph twoway (connected X time, msymbol(i) xline(100)) 




	




**********************************
* Figures 8.2, 8.3, and 8.4
**********************************


* Deterministic
	drop _all
	clear all
	set obs 100
	set seed 3456
	gen t = _n -1
	tsset t
	gen DL = 0
	replace DL = 1 if t > 50
	gen DP = 0
	replace DP = 1 if t == 51
	gen error = 0
	gen Ynull = 1 in 1
	gen Yalt  = 1 in 1
	replace Ynull = 1 + L.Ynull + 10*DP + error in 2/L
	replace Yalt = 1 + t + 10*DL + error in 2/L
	graph twoway (connected Ynull Yalt t, msymbol(i i) xline(50))



* Stochastic
	drop _all
	clear all
	set obs 100
	set seed 3456
	gen t = _n-1
	tsset t
	gen DL = 0
	replace DL = 1 if t > 50
	gen DP = 0
	replace DP = 1 if t == 51
	gen error = rnormal()
	gen Ynull = 1 in 1
	gen Yalt  = 1 in 1
	replace Ynull = 1 + L.Ynull + 10*DP + error in 2/L
	replace Yalt = 1 + t + 10*DL + error in 2/L

	graph twoway (connected Yalt t, msymbol(i) xline(50) ylabel(0(40)120) yscale(range(0 120)) )
	graph twoway (connected Ynull t, msymbol(i) xline(50) ylabel(0(40)120) yscale(range(0 120)) ) 

