**********************************************
* ARMA
**********************************************

cd "C:\Users\johnlevendis\Google Drive\papers\book\ARMA\"
set scheme s1mono

**********************************************
* Stationarity
**********************************************

* Returns are, roughly, mean stationary but not vaiance stationary.

	* Initialize
		drop _all
		local observations 300
		set obs `observations'
		set seed 12345
		gen time = _n
		tsset time
	
	* (1) Generate the stdev of X as a simple regime-switching type model of high and low volatility
		gen sigma = 1
		replace sigma = 10 if (time < 20) | (time > 50 & time < 60) | (time > 90 & time < 130) | (time > 190 & time < 250)
	
	* (2) Generate the data, drawing errors from N(0,sigma)
		gen X = 0
		forvalues t = 2/`observations'{
			quietly sum sigma in `t'
			quietly local sigma = r(mean)
			quietly replace X = rnormal(0,`sigma') in `t'
		}
		graph twoway (line X time, ylabel(none))
		graph save   "graphs\mean_stationary_but_not_variance_stationary"    , replace
		graph export "graphs\mean_stationary_but_not_variance_stationary.pdf", replace
			
			
			
* Example of time-series that is not covariance-stationary

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
		graph save   "graphs\not_covariance_stationary"    , replace
		graph export "graphs\not_covariance_stationary.pdf", replace
		


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

* Estimating an AR(1) process
	use ".\data\ARexamples.dta", clear
	reg X L.X, nocons
	arima X, ar(1) nocons nolog


* IRFs of an AR(1) process:
	drop _all
	clear all
	set obs 10
	gen time = _n 	// Create a time variable
	tsset time		// Tell Stata that it is a time variable
	gen X= 0 		// Starting from an initial value of zero
	replace X = 1 in 1	// X takes a shock of 1 unit at t=1
	replace X= 0.75*L.X in 2/L
	graph twoway line X time

	
* IRFs of an AR(1) process (looks like step-ladder)
	drop _all
	clear all
	local observations = 10		//Enter the number of observations here
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
	graph export ".\graphs\ar1IRF.pdf", replace


* Automatic IRFs 
	use ".\data\ARexamples.dta", clear
	arima X, ar(1) nocons
	irf create AR1, step(10) set(name)
	irf graph irf
	irf drop AR1
	graph export ".\graphs\ar1irfauto.pdf", replace


	
	
**********************************************
* AR(p) processes
**********************************************	
	
* Estimating an AR(p) process
	regress X L.X L2.X L3.X L4.X, nocons
	arima X, ar(1/4) nocons
	
* Estimation Exercises:
	use ".\data\ARexamples.dta", clear
	arima X, ar(1) nocons nolog
	arima Y, ar(1/2) nocons nolog
	arima Z, ar(1/3) nocons nolog
	arima Z, ar(1 3) nocons nolog


* Forecasting with AR(p) models
	use ".\data\ARexamples.dta", clear
	arima Z, ar(1/3) nocons
	list Z in 2998/L /* The data in the last five periods */
	tsappend , add(4)
	list Z in 2998/L
	predict Zhat
	list Z Zhat in 2998/L	
	* Verifying:
	*display _b[ARMA:L1.ar]*(.58345208) +  _b[ARMA:L2.ar]*( 1.1054195) +_b[ARMA:L3.ar]*(.68069001) 
	*display _b[ARMA:L1.ar]*(.63622778) +  _b[ARMA:L2.ar]*(.58345208) +_b[ARMA:L3.ar]*( 1.1054195) 
	*display _b[ARMA:L1.ar]*(.61675564) +  _b[ARMA:L2.ar]*(.63622778) +_b[ARMA:L3.ar]*(.58345208) 
	*display _b[ARMA:L1.ar]*(.55970004) +  _b[ARMA:L2.ar]*(.61675564) +_b[ARMA:L3.ar]*(.63622778) 


	
	


**********************************************
* MA processes
**********************************************	

* Forecasting with MA(1) models
	drop _all
	clear all
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

* Estimation of MA(q) models
	use ".\data\MAexamples.dta", clear
	arima Y, ma(1/3) nocons nolog

* Estimation exercises
	use ".\data\MAexamples.dta", clear
	arima Y, ma(1/3) nocons nolog
	arima Y, ma(1/2) nocons nolog
	arima Z, ma(1/3) nocons nolog
		
		
		


		
				
* It is often difficult to tell visually whether a time-series is an AR or an MA process. 
* Consider the following:
	drop _all
	clear all
	set seed 123
	sim_arma X1, nobs(30) ar(0.50) spin(2000) time(time)
	sim_arma X2, nobs(30) ma(0.50) spin(2000) time(time)
	sim_arma X3, nobs(30) ar(0.25) ma(0.25) spin(2000) time(time)
	sim_arma X4, nobs(30) ar(0.50) ma(0.50) spin(2000) time(time)

	label var X1 "AR(1): beta = 0.50"
	label var X2 "MA(1): beta = 0.50"
	label var X3 "ARMA(1,1): beta1 = 0.25, gamma1 = 0.25"
	label var X4 "ARMA(1,1): beta1 = 0.50, gamma1 = 0.50"

	line X1 time
		graph save ".\graphs\x1.gph", replace
	line X2 time
		graph save ".\graphs\x2.gph", replace
	line X3 time
		graph save ".\graphs\x3.gph", replace
	line X4 time
		graph save ".\graphs\x4.gph", replace
			
	graph combine ".\graphs\x1.gph" ".\graphs\x2.gph" ".\graphs\x3.gph" ".\graphs\x4.gph", rows(2) cols(2) xcommon 
	erase ".\graphs\x1.gph"
	erase ".\graphs\x2.gph"
	erase ".\graphs\x3.gph"
	erase ".\graphs\x4.gph"
	graph export ".\graphs\x1x2x3x4.pdf", replace





*************************************
* GENERATING NON-ZERO MEAN AR(p) DATA
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

	graph twoway line X time if time> 9900
	graph export ".\graphs\AR1nonzero.pdf", replace

	
* Two different approaches to dealing with non-zero processes
	* (1) Just estimate the constant
		arima X, ar(1) nolog     // We leave out the ``nocons'' option
	* (2) Manually de-mean
		sum X
		local mean = r(mean)
		gen X_demeaned = X - `mean'
		arima X_demeaned, ar(1) nolog nocons  //De-mean, but include ``nocons''



	
	

**********************************************
* Model selection by Information criteria
**********************************************	

	use ".\data\ARexamples.dta", clear
	* we know that it is an AR(1) model

	quietly arima X, ar(1/3) nocons
	estat ic
	quietly arima X, ar(1/2) nocons
	estat ic
	quietly arima X, ar(1) nocons
	estat ic
	
	quietly arima X, ma(1/3) nocons
	estat ic
	quietly arima X, ma(1/2) nocons
	estat ic
	quietly arima X, ma(1) nocons
	estat ic
	








************************************
* Calculating Empirical Partial ACFs
************************************

* Empirical PACFs of AR processes
		*cd "C:\Documents and Settings\johnlevendis\Desktop\"
		drop _all
		clear all
		use ".\data\ARexamples.dta", clear

	* By hand
		reg X L.X
		reg X L.X L2.X
		reg X L.X L2.X L3.X
		reg X L.X L2.X L3.X L4.X
		
	* Making the table
		qui reg X L.X
		estimates store temp1
		qui reg X L.X L2.X
		estimates store temp2
		qui reg X L.X L2.X L3.X
		estimates store temp3
		qui reg X L.X L2.X L3.X L4.X
		estimates store temp4
		est table temp1 temp2 temp3 temp4
	* Automatically
		corrgram X, lags(4)
		ac X
			graph export acfXar.pdf, replace
		pac X
			graph export pacfXar.pdf, replace
		
	* By hand
		reg Y L.Y
		reg Y L.Y L2.Y
		reg Y L.Y L2.Y L3.Y
		reg Y L.Y L2.Y L3.Y L4.Y
	* Automatically 
		corrgram Y, lags(4)
		pac Y
			graph export pacfYar.pdf, replace
		ac Y
			graph export acfYar.pdf, replace

			
			
* Empirical PACFs of MA processes
		*cd "C:\Documents and Settings\johnlevendis\Desktop\"
		drop _all
		clear all
		use ".\data\MAexamples.dta", clear

		reg X L.X
		reg X L.X L2.X
		reg X L.X L2.X L3.X
		reg X L.X L2.X L3.X L4.X

		reg Y L.Y
		reg Y L.Y L2.Y
		reg Y L.Y L2.Y L3.Y
		reg Y L.Y L2.Y L3.Y L4.Y

		corrgram X, lags(4)
		corrgram Y, lags(4)

		pac X
			graph export pacfXma.pdf, replace
		pac Y
			graph export pacfYma.pdf, replace
		ac X
			graph export acfXma.pdf, replace
		ac Y
			graph export acfYma.pdf, replace

	

		reg X L.X
		reg X L.X L2.X
		reg X L.X L2.X L3.X
		reg X L.X L2.X L3.X L4.X

		reg Y L.Y
		reg Y L.Y L2.Y
		reg Y L.Y L2.Y L3.Y
		reg Y L.Y L2.Y L3.Y L4.Y

		corrgram X, lags(4)
		corrgram Y, lags(4)

		ac X
			graph export acfXar.pdf, replace
		pac X
			graph export pacfXar.pdf, replace
		ac Y
			graph export acfYar.pdf, replace
		pac Y
			graph export pacfYar.pdf, replace

			
			
* Empirical PACFs of MA processes
		*cd "C:\Documents and Settings\johnlevendis\Desktop\"
		drop _all
		clear all
		use ".\data\MAexamples.dta", clear

		reg X L.X
		reg X L.X L2.X
		reg X L.X L2.X L3.X
		reg X L.X L2.X L3.X L4.X

		reg Y L.Y
		reg Y L.Y L2.Y
		reg Y L.Y L2.Y L3.Y
		reg Y L.Y L2.Y L3.Y L4.Y

		corrgram X, lags(4)
		corrgram Y, lags(4)

		ac X
			graph export acfXma.pdf, replace
		pac X
			graph export pacfXma.pdf, replace
		ac Y
			graph export acfYma.pdf, replace
		pac Y
			graph export pacfYma.pdf, replace

	
	










************************************
* Putting it all together
************************************

* Example 1: a toy example
	drop _all
	set seed 2468
	sim_arma X, ma(0.70 0.10) nobs(5000) spin(1000) time(time)
	ac X, lags(20)
			graph export ".\graphs\example_1_ma_acf.pdf", replace
	pac X, lags(20)
			graph export ".\graphs\example_1_ma_pacf.pdf", replace
	arima X, ma(1/2) nolog
	arima X, ma(1/2) nolog nocons  

	
* Example 2: GDP growth rate

	*freduse A191RL1Q225SBEA, clear 
	* Real GDP growth rate, quarterly and sesonally adjusted
	* 1947-04-01 through 2017-04-01
	*rename A19 rGDPgr
	label var rGDPgr "Real GDP growth rate"
	*save ".\data\rGDPgr.dta", replace
	use ".\data\rGDPgr.dta", clear
	gen time = _n
	tsset time
	ac rGDPgr, lags(20)
		graph export ".\graphs\rGDPgr_ACF.pdf", replace
	pac rGDPgr, lags(20)
		graph export ".\graphs\rGDPgr_PACF.pdf", replace
	arima rGDPgr, ar(1) nolog
	tsappend, add(5)
	predict rGDPgrhat
	twoway (connected rGDPgr time if time >=250, msymbol(none)) ///
	(connected rGDPgrhat time if time>= 250 & rGDPgr==., msymbol(none) lpattern(dash))
	graph export ".\graphs\rGDPgr_AR1_forecast.pdf", replace

	
	

************************************
* Putting it all together - EXERCISES
************************************

* Generate some pretend data:
	drop _all
	set seed 9876
	sim_arma X1, ar(0.70) nobs(5000) spin(1000) time(time)
		sort time
		replace X1 = X1 + 2
		save ".\data\ARMAexercises1a.dta", replace
	sim_arma X2, ma(-0.60 0.10) nobs(5000) spin(1000) time(time)
		sort time
		save ".\data\ARMAexercises1b.dta", replace
	sim_arma X3, ar(0.70 0.20) nobs(5000) spin(1000) time(time)
		sort time
		save ".\data\ARMAexercises1c.dta", replace
	sim_arma X4, ma(-0.60 0.10) nobs(5000) spin(1000) time(time)
		sort time
		replace X4 = X4 + 3
		save ".\data\ARMAexercises1d.dta", replace
* Merge the datasets
	use ".\data\ARMAexercises1a.dta", clear
	merge 1:1 time using ".\data\ARMAexercises1b.dta"
		drop _merge
	merge 1:1 time using ".\data\ARMAexercises1c.dta"
		drop _merge
	merge 1:1 time using ".\data\ARMAexercises1d.dta"
		drop _merge
	order time X1 X2 X3 X4
	save ".\data\ARMAexercises.dta", replace
* Clean up
	erase ".\data\ARMAexercises1a.dta"
	erase ".\data\ARMAexercises1b.dta"
	erase ".\data\ARMAexercises1c.dta"
	erase ".\data\ARMAexercises1d.dta"
* Solution
	arima X1, ar(1) nocons
	arima X2, ma(1/2)
	arima X3, ar(1/2)
	arima X4, ma(1/2) nocons















	
	
	
	
	
	
	
	
	
	
	
