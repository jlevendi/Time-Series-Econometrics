**********************************************
* ARMA
**********************************************

set scheme s1mono

**********************************************
* Stationarity
**********************************************

* Figure 2.1
		drop _all
		local observations 300
		set obs `observations'
		set seed 12345
		gen time = _n
		tsset time

		gen sigma = 1
		replace sigma = 10 if (time < 20) | (time > 50 & time < 60) | (time > 90 & time < 130) | (time > 190 & time < 250)

		gen X = 0
		forvalues t = 2/`observations'{
			quietly sum sigma in `t'
			quietly local sigma = r(mean)
			quietly replace X = rnormal(0,`sigma') in `t'
		}
		graph twoway (line X time, ylabel(none))
		
			
			
			
* Figure 2.2
		drop _all
		set obs 300
		set seed 12345
		gen time = _n
		tsset time
		
		gen error = rnormal(0,1)
		gen X = 0
		replace X = 0.5*L.X + error in 2/L
		replace X = 0  in 100/200
		
		graph twoway (line X time, ylabel(none))
		


**********************************************
* Generating fake data for illustration purposes
**********************************************

*AR(1), AR(2) AND AR(3) PROCESSES
	drop _all
	clear all
	set seed 1234

	sim_arma X, nobs(3000) ar(0.50)           spin(2000) time(time)
	sim_arma Y, nobs(3000) ar(0.70 0.20)      spin(2000) time(time)
	sim_arma Z, nobs(3000) ar(0.60 0.20 0.10) spin(2000) time(time) 

	arima X, ar(1) nocons nolog
	arima Y, ar(1/2) nocons nolog
	arima Z, ar(1/3) nocons nolog
	
	order time X Y Z
	save ".\data\ARexamples.dta", replace

	
* MA(1), MA(2) AND MA(3) PROCESSES
	drop _all
	clear all
	set seed 1234
	
	sim_arma X, nobs(3000) ma(0.50)           spin(2000) time(time) 
	sim_arma Y, nobs(3000) ma(0.70 0.20)      spin(2000) time(time) 
	sim_arma Z, nobs(3000) ma(0.60 0.20 0.10) spin(2000) time(time) 
	
	arima X, ma(1) nocons nolog
	arima X, ma(1/2) nocons nolog
	arima X, ma(1/3) nocons nolog
	
	ac X, lags(10)
	ac Y, lags(10)
	ac Z, lags(10)
	
	pac X, lags(10)
	pac Y, lags(10)
	pac Z, lags(10)

	order time X Y Z
	save ".\data\MAexamples.dta", replace   

	
**********************************************
* AR(1) processes
**********************************************

***********************************
* 2.2.1 Estimating an AR(1) model
***********************************

* Estimating an AR(1) process
	use ".\data\ARexamples.dta", clear
	reg X L.X, nocons
	arima X, ar(1) nocons nolog

***********************************
* 2.2.2 Impulse Responses	
***********************************
	
* Figure 2.5
	drop _all
	clear all
	local observations = 10		
	set obs `observations'
	gen time = _n 
	tsset time
	gen X= 0 
	replace X = 1 in 1
	replace X= 0.75*L.X in 2/L
	*browse
	graph twoway line X time

	set obs `=_N+1'
	replace time = .999 in L
	replace X = 0 in L
	sort time
	graph twoway line X time

	set obs `=_N+1'
	replace time = 0 in L
	replace X = 0 in L
	sort time
	graph twoway line X time

	
* Figure 2.6
	use ".\data\ARexamples.dta", clear
	arima X, ar(1) nocons
	irf create AR1, step(10) set(name)
	irf graph irf
	irf drop AR1



	
	
**********************************************
* 2.3.1 Estimating an AR(p) Model
**********************************************	
	
* Estimating an AR(p) process
	regress X L.X L2.X L3.X L4.X, nocons
	arima X, ar(1/4) nocons
	
* Estimation Exercises (3), (4) and (5):
	use ".\data\ARexamples.dta", clear
	arima X, ar(1) nocons nolog
	arima Y, ar(1/2) nocons nolog
	arima Z, ar(1/3) nocons nolog
	arima Z, ar(1 3) nocons nolog

**********************
* 2.3.3 Forecasting
**********************
	use ".\data\ARexamples.dta", clear
	arima Z, ar(1/3) nocons
	list Z in 2998/L /* The data in the last five periods */
	tsappend , add(4)
	list Z in 2998/L
	predict Zhat
	list Z Zhat in 2998/L	
	
* Exercise (2)
	display _b[ARMA:L1.ar]*(.58345208) +  _b[ARMA:L2.ar]*( 1.1054195) +_b[ARMA:L3.ar]*(.68069001) 
	display _b[ARMA:L1.ar]*(.63622778) +  _b[ARMA:L2.ar]*(.58345208) +_b[ARMA:L3.ar]*( 1.1054195) 
	display _b[ARMA:L1.ar]*(.61675564) +  _b[ARMA:L2.ar]*(.63622778) +_b[ARMA:L3.ar]*(.58345208) 
	display _b[ARMA:L1.ar]*(.55970004) +  _b[ARMA:L2.ar]*(.61675564) +_b[ARMA:L3.ar]*(.63622778) 


	
	


**********************************************
* 2.4 MA(1) Models
* 2.4.3 Forecasting
**********************************************	

* Forecasting with MA(1) models
	use ".\data\MAexamples.dta", clear
	*keep X time
	arima X, ma(1) nocons nolog
	predict Xhat         /*Calculate the fitted values*/
	predict r, residuals /*Calculate the residuals*/
	list X Xhat r in 2996/3000 /* Listing the last 5 periods*/
	tsappend, add(1)     /*Adding a blank observation*/
	list X Xhat r in 2996/3001 /* Listing the last 6 periods*/
	drop Xhat
	predict Xhat         /*Re-calculating the fitted values*/
	list X Xhat r in 2996/3001 /* Listing the last 6 periods */

	
	
	
**********************************************
* 2.5 MA(q) Models
* 2.5.1 Estimation
**********************************************	
	
	
* Estimation of MA(q) models
	use ".\data\MAexamples.dta", clear
	arima Y, ma(1/3) nocons nolog

* Exercises (1) and (2)
	use ".\data\MAexamples.dta", clear
	arima Y, ma(1/3) nocons nolog
	arima Y, ma(1/2) nocons nolog
	arima Z, ma(1/3) nocons nolog
		
		
	


*************************************
* 2.6 Non-zero ARMA processes
*************************************

* Non-zero AR(1) process
	set more off
	drop _all
	clear all
	set obs 10000
	set seed 1234
	gen time = _n
	tsset time
	gen double error = rnormal()
	label var e "errors"

	gen double X = error in 1
	replace X = 10 + 0.50*L.X + error in 2/L
	drop if _n <= 1000
	label var X "X = 10 + 0.50*L.X + e"
	save ".\data\AR1nonzero.dta", replace

	graph twoway line X time if time> 9900  /* Figure 2.7 */

	
* Two different approaches to dealing with non-zero processes
	* (1) Just estimate the constant
		arima X, ar(1) nolog     // We leave out the ``nocons'' option
	* (2) Manually de-mean
		sum X
		local mean = r(mean)
		gen X_demeaned = X - `mean'
		arima X_demeaned, ar(1) nolog nocons  //De-mean, but include ``nocons''



	
	