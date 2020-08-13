
****************************************************
* Chapter 3: Model Selection in ARMA(p,q) processes
****************************************************
	

	
* Figure 3.1				
	drop _all
	clear all
	set seed 123
	sim_arma X1, nobs(30) ar(0.50) spin(2000) time(time)
	sim_arma X2, nobs(30) ma(0.50) spin(2000) time(time)
	sim_arma X3, nobs(30) ar(0.25) ma(0.25) spin(2000) time(time)
	sim_arma X4, nobs(30) ar(0.50) ma(0.50) spin(2000) time(time)

	label var X1 "AR(1): beta = 0.50"
	line X1 time

	label var X2 "MA(1): beta = 0.50"
	line X2 time
		
	label var X3 "ARMA(1,1): beta1 = 0.25, gamma1 = 0.25"
	line X3 time
	
	label var X4 "ARMA(1,1): beta1 = 0.50, gamma1 = 0.50"
	line X4 time






************************************
* 3.2.2 Calculating Empirical PACFs
************************************

* Empirical PACFs of AR processes
		use ".\data\ARexamples.dta", clear

	* By hand
		reg X L.X
		reg X L.X L2.X
		reg X L.X L2.X L3.X
		reg X L.X L2.X L3.X L4.X
		
	* Table 3.1
		qui reg X L.X
		estimates store reg1
		qui reg X L.X L2.X
		estimates store reg2
		qui reg X L.X L2.X L3.X
		estimates store reg3
		qui reg X L.X L2.X L3.X L4.X
		estimates store reg4
		est table reg1 reg2 reg3 reg4
	* Automatically
		corrgram X, lags(4)
		ac X   /* Figure 3.16a */
		pac X  /* Figure 3.16b */

* Next, we complete the same type of exercise, but with data from an MA process.		
	* By hand
		reg Y L.Y
		reg Y L.Y L2.Y
		reg Y L.Y L2.Y L3.Y
		reg Y L.Y L2.Y L3.Y L4.Y
	* Automatically 
		corrgram Y, lags(4)
		ac Y 	/* Figure 3.17a */
		pac Y 	/* Figure 3.17b */

			

		
			
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

		ac X 	/* Figure 3.18a */
		pac X 	/* Figure 3.18b */
		
		ac Y 	/* Figure 3.19a */
		pac Y 	/* Figure 3.19b */
		

	
	



************************************
* 3.3 Putting it all together
************************************

* Example
	drop _all
	set seed 2468
	sim_arma X, ma(0.70 0.10) nobs(5000) spin(1000) time(time)
	ac X, lags(20)		/* Figure 3.20a */
	pac X, lags(20)		/* Figure 3.20b */
	arima X, ma(1/2) nolog
	arima X, ma(1/2) nolog nocons  

	
* Example
	*freduse A191RL1Q225SBEA, clear 
	*rename A19 rGDPgr
	*label var rGDPgr "Real GDP growth rate"
	*gen time = _n
	*save ".\data\rGDPgr.dta", replace
	use ".\data\rGDPgr.dta", clear
	tsset time
	ac rGDPgr, lags(20)		/* Figure 3.21a */
	pac rGDPgr, lags(20)	/* Figure 3.21b */
	arima rGDPgr, ar(1) nolog
	tsappend, add(5)
	predict rGDPgrhat
	* Figure 3.22:
	twoway (connected rGDPgr time if time >=250, msymbol(none)) ///
		(connected rGDPgrhat time if time>= 250 & rGDPgr==., msymbol(none) lpattern(dash))
	

* 3.3 Exercises
	use ".\data\ARMAexercises.dta", replace
	arima X1, ar(1) nocons
	arima X2, ma(1/2)
	arima X3, ar(1/2)
	arima X4, ma(1/2) nocons






**********************************************
* 3.4 Information Criteria
**********************************************	

	use ".\data\ARexamples.dta", clear
	* we know that it is an AR(1) model

	quietly arima X, ar(1/3) nocons
	estat ic
	quietly arima X, ar(1/2) nocons
	estat ic
	quietly arima X, ar(1) nocons
	estat ic
	* How do these compare to the MA models?
	quietly arima X, ma(1/3) nocons
	estat ic
	quietly arima X, ma(1/2) nocons
	estat ic
	quietly arima X, ma(1) nocons
	estat ic
	













	
	
	
	
	
	
	
	
	
	
	
