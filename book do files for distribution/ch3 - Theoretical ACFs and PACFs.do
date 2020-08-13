

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





**********************************************
* 3.1 Theoretical ACFs and PACFs
**********************************************

cd ".\graphs\"
set scheme s1mono

**********************************************
* ACF of AR(p)
**********************************************

set more off
drop _all
clear all
set obs 1000
set seed 1234
gen time = _n
tsset time
gen X = rnormal()
save "temp_data.dta", replace

* ACF of AR(1)
	use "temp_data.dta", clear
	constraint define 1 _b[ARMA:L1.ar] =  0.50
	constraint define 2 _b[ARMA:L1.ar] = -0.50
	
	quietly arima X , ar(1) noconstant nolog constraints(1)
		estat acplot, lags(10) recast(spike) saving("temp", replace) 
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar1_a.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	quietly arima X , ar(1) noconstant nolog constraints(2)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar1_b.pdf", replace
		erase "temp.dta"


* ACF of AR(2)
	constraint define 1 _b[ARMA:L1.ar] =  0.50
	constraint define 2 _b[ARMA:L1.ar] = -0.50
	constraint define 3 _b[ARMA:L2.ar] =  0.20
	constraint define 4 _b[ARMA:L2.ar] = -0.20
	use "temp_data.dta", clear

	quietly arima X, ar(1/2) noconstant nolog constraints(1 3)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0		
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar2_a.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	quietly arima X, ar(1/2) noconstant nolog constraints(1 4)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar2_b.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	quietly arima X, ar(1/2) noconstant nolog constraints(2 4)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar2_c.pdf", replace
		erase "temp.dta"


	
**********************************************
* ACF of MA(q)
**********************************************

* ACF of MA(1)
	constraint define 1 _b[ARMA:L1.ma] =  0.50
	constraint define 2 _b[ARMA:L1.ma] = -0.50
	use "temp_data.dta", clear
	
	arima X, ma(1) noconstant nolog constraints(1)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma1_a.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	arima X, ma(1) noconstant nolog constraints(2)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma1_b.pdf", replace
		erase "temp.dta"

	
* ACF of MA(2)
	constraint define 1 _b[ARMA:L1.ma] =  0.50
	constraint define 2 _b[ARMA:L1.ma] = -0.50
	constraint define 3 _b[ARMA:L2.ma] =  0.20
	constraint define 4 _b[ARMA:L2.ma] = -0.20
	use "temp_data.dta", clear
	
	quietly arima X, ma(1/2) noconstant nolog constraints(1 3)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma2_a.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	quietly arima X, ma(1/2) noconstant nolog constraints(1 4)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma2_b.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	quietly arima X, ma(1/2) noconstant nolog constraints(2 4)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma2_c.pdf", replace
		erase "temp.dta"
		


**********************************************
* ACF of ARMA(1,1)
**********************************************

	constraint define 1 _b[ARMA:L1.ar] =  0.40
	constraint define 2 _b[ARMA:L1.ma] =  0.40
	constraint define 3 _b[ARMA:L1.ar] = -0.40
	constraint define 4 _b[ARMA:L1.ma] = -0.40
	use "temp_data.dta", clear
	
	arima X, ar(1) ma(1) noconstant nolog constraints(1 2)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar1ma1_a.pdf", replace
		erase "temp.dta"
		use "temp_data.dta", clear
	arima X, ar(1) ma(1) noconstant nolog constraints(3 4)
		estat acplot, lags(10) recast(spike) saving("temp", replace)
		use "temp", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ar1ma1_b.pdf", replace
		erase "temp.dta"
	

**********************************************
* Summary of above
**********************************************

* Theoretical ACF of AR(1)
	* X_t =  0.50 X_{t-1} + e_t ==> "theoretical_acf_ar1_a.pdf"
	* X_t = -0.50 X_{t-1} + e_t ==> "theoretical_acf_ar1_b.pdf"
	
	
* Theoretical ACF of AR(2)
	* X_t =  0.50 X_{t-1} +  0.20 X_{t-2} + e_t ==> "theoretical_acf_ar2_a.pdf"
	* X_t =  0.50 X_{t-1} + -0.20 X_{t-2} + e_t ==> "theoretical_acf_ar2_b.pdf"
	* X_t = -0.50 X_{t-1} + -0.20 X_{t-2} + e_t ==> "theoretical_acf_ar2_c.pdf"


* Theoretical ACF of MA(1)
	* X_t = u_t +  0.50 u_{t-1} ==> "theoretical_acf_ma1_a.pdf"
	* X_t = u_t + -0.50 u_{t-1}	==> "theoretical_acf_ma1_b.pdf"


* Theoretical ACF of MA(2)
	*X_t = u_t +  0.50 u_{t-1} +  0.20 u_{t-2} ==> "theoretical_acf_ma2_a.pdf"
	*X_t = u_t +  0.50 u_{t-1} + -0.20 u_{t-2} ==> "theoretical_acf_ma2_b.pdf"
	*X_t = u_t + -0.50 u_{t-1} + -0.20 u_{t-2} ==> "theoretical_acf_ma2_c.pdf"


* Theoretical ACF of ARMA(1,1)
	* X_t =  0.40 X_{t-1} + u_t +  0.40 u_{t-1} ==> "theoretical_acf_ar1ma1_a.pdf"
	* X_t = -0.40 X_{t-1} + u_t + -0.40 u_{t-1} ==> "theoretical_acf_ar1ma1_b.pdf"
	

	
	
	
	
	
	
	

**********************************************
* ACF of MA(q) (by hand)
**********************************************
* Verifying that what I calculated by hand in the text is correct
* X_t = u + beta1*L.u + beta2*LL.u + beta3*LLL.u
*     = u +  0.40*L.u +  0.20*LL.u +  0.10*LLL.u

	use "temp_data.dta", clear
	constraint define 1 _b[ARMA:L1.ma] = 0.40
	constraint define 2 _b[ARMA:L2.ma] = 0.20
	constraint define 3 _b[ARMA:L3.ma] = 0.10
	
	arima X, ma(1/3) noconstant nolog constraints(1 2 3)
		estat acplot, lags(10) saving("theoretical_acf_ma3", replace)
		use "theoretical_acf_ma3.dta", clear
		replace lag = lag-1
		drop if lag == 0
		twoway dropline ac lag, msymbol(O) yline(0)
		graph export "theoretical_acf_ma3.pdf", replace
		erase "theoretical_acf_ma3.dta"




	
	




**********************************************
* Theoretical PACFs
**********************************************

/*

* Approach: Generate a huge AR(1) dataset, so that the empirical PACF is spot-on with the Theoretical PACF


* Theoretical PACF of AR(1)
	drop _all
	local obs = 15000000	// 15 million

	sim_arma X, nobs(`obs') ar(0.50) spin(2000) time(time)
	corrgram X, lags(10)
	drop X

	sim_arma X, nobs(`obs') ar(-0.50) spin(2000) time(time)
	corrgram X, lags(10)
	drop X

	* The Th-PACF of an AR(1) process is equal to \beta_1, and zeros everywhere else.

	

* Theoretical PACF of AR(2)
	sim_arma X, nobs(`obs') ar(0.50 0.20) spin(2000) time(time)
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ar(0.50 -0.20) spin(2000) time(time)
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ar(-0.50 -0.20) spin(2000) time(time)
		corrgram X, lags(10)
		drop X
		
	
* Theoretical PACF of MA(1)
	sim_arma X, nobs(`obs') ma(0.50) spin(2000) time(time)
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ma(-0.50) spin(2000) time(time)	
		corrgram X, lags(10)
		drop X
		
		
		
* Theoretical PACF of MA(2)
	sim_arma X, nobs(`obs') ma(0.50 0.20) spin(2000) time(time)
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ma(0.50 -0.20) spin(2000) time(time)	
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ma(-0.50 -0.20) spin(2000) time(time)	
		corrgram X, lags(10)
		drop X


* Theoretical PACF of ARMA(1,1)
	sim_arma X, nobs(`obs') ar(0.40) ma(0.40) spin(2000) time(time)	
		corrgram X, lags(10)
		drop X
	sim_arma X, nobs(`obs') ar(-0.40) ma(-0.40) spin(2000) time(time)	
		corrgram X, lags(10)
		drop X

*/		
		
		
		
		
		
**********************************************
* Drawing Theoretical PACFs
**********************************************
local obs = 15000000	// 15 million

*sim_arma X, nobs(`nobs') ar(0.50) spin(2000) time(time)
capture drop _all
input AC PAC
	0.4997   0.4997   
	0.2496  -0.0001   
	0.1243  -0.0005   
	0.0620   0.0001   
	0.0309   0.0001   
	0.0155   0.0000   
	0.0079   0.0002   
	0.0041   0.0001   
	0.0023   0.0002   
	0.0013   0.0001   
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar1_a.pdf", replace




*sim_arma X, nobs(`obs') ar(-0.50) spin(2000) time(time)
capture drop _all
input AC PAC
	-0.4998  -0.4998
	0.2497  -0.0002
	-0.1249  -0.0002
	0.0625   0.0001
	-0.0311   0.0002
	0.0152  -0.0004
	-0.0073   0.0002
	0.0033  -0.0002
	-0.0013   0.0003
	0.0005   0.0001
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar1_b.pdf", replace



*sim_arma X, nobs(`obs') ar(0.50 0.20) spin(2000) time(time)
capture drop _all
input AC PAC
	0.6250   0.6250
	0.5126   0.2002
	0.3814   0.0001
	0.2933  -0.0000
	0.2228  -0.0002
	0.1699  -0.0002
	0.1292  -0.0004
	0.0982  -0.0004
	0.0748   0.0001
	0.0569   0.0002
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar2_a.pdf", replace



* sim_arma X, nobs(`obs') ar(0.50 -0.20) spin(2000) time(time)
capture drop _all
input AC PAC
	0.4168   0.4168
	0.0085  -0.2000
	-0.0793  -0.0002
	-0.0414   0.0001
	-0.0053  -0.0005
	0.0056   0.0002
	0.0039   0.0000
	0.0007  -0.0002
	-0.0008  -0.0003
	-0.0005   0.0002
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar2_b.pdf", replace


                 
*sim_arma X, nobs(`obs') ar(-0.50 -0.20) spin(2000) time(time)
capture drop _all
input AC PAC
	-0.4164  -0.4164
	0.0084  -0.1996
	0.0788  -0.0001
	-0.0408   0.0003
	0.0046   0.0001
	0.0055  -0.0004
	-0.0036  -0.0002
	0.0006  -0.0002
	0.0008   0.0004
	-0.0006   0.0001
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar2_c.pdf", replace



* sim_arma X, nobs(`obs') ma(0.50) spin(2000) time(time)
capture drop _all
input AC PAC
	0.3998   0.3998
	-0.0002  -0.1904
	0.0000   0.0941
	0.0001  -0.0468
	0.0002   0.0235
	0.0002  -0.0115
	0.0002   0.0059
	0.0001  -0.0030
	-0.0000   0.0015
	-0.0004  -0.0013
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ma1_a.pdf", replace


                 
* sim_arma X, nobs(`obs') ma(-0.50) spin(2000) time(time) 
capture drop _all
input AC PAC
	-0.3999  -0.3999
	-0.0002  -0.1907
	0.0006  -0.0936
	-0.0006  -0.0471
	0.0003  -0.0234
	0.0003  -0.0113
	-0.0003  -0.0058
	-0.0000  -0.0029
	0.0004  -0.0010
	-0.0006  -0.0011
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ma1_b.pdf", replace




* sim_arma X, nobs(`obs') ma(0.50 0.20) spin(2000) time(time)
capture drop _all
input AC PAC
	0.4653   0.4653
	0.1555  -0.0780
	0.0004  -0.0531
	0.0001   0.0420
	-0.0002  -0.0107
	-0.0005  -0.0035
	-0.0004   0.0037
	0.0003  -0.0006
	0.0004  -0.0002
	0.0005   0.0007
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ma2_a.pdf", replace



                 
* sim_arma X, nobs(`obs') ma(0.50 -0.20) spin(2000) time(time)    
capture drop _all
input AC PAC
	0.3104   0.3104
	-0.1545  -0.2776
	0.0001   0.1769
	-0.0001  -0.1346
	0.0003   0.0992
	0.0003  -0.0747
	-0.0004   0.0559
	-0.0005  -0.0429
	0.0000   0.0326
	0.0004  -0.0245
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ma2_b.pdf", replace


 
* sim_arma X, nobs(`obs') ma(-0.50 -0.20) spin(2000) time(time)   
capture drop _all
input AC PAC
	-0.3103  -0.3103
	-0.1548  -0.2779
	-0.0002  -0.1772
	-0.0002  -0.1351
	0.0006  -0.0986
	-0.0004  -0.0750
	0.0004  -0.0560
	-0.0001  -0.0427
	-0.0006  -0.0331
	0.0005  -0.0247
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ma2_c.pdf", replace



* sim_arma X, nobs(`obs') ar(0.40) ma(0.40) spin(2000) time(time) 
capture drop _all
input AC PAC
	0.6270   0.6270
	0.2509  -0.2344
	0.1003   0.0926
	0.0398  -0.0375
	0.0157   0.0150
	0.0062  -0.0060
	0.0026   0.0026
	0.0012  -0.0009
	0.0005   0.0002
	0.0002  -0.0001
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar1ma1_a.pdf", replace



* sim_arma X, nobs(`obs') ar(-0.40) ma(-0.40) spin(2000) time(time)       
capture drop _all
input AC PAC
	-0.6274  -0.6274
	0.2511  -0.2350 
	-0.1004  -0.0929
	0.0400  -0.0372 
	-0.0158  -0.0146
	0.0062  -0.0059 
	-0.0026  -0.0027
	0.0014  -0.0006 
	-0.0010  -0.0006
	0.0004  -0.0006 
end
gen lag = _n
label var PAC "Partial Autocorrelations"
twoway dropline PAC lag, msymbol(O) yline(0)
graph export "theoretical_pacf_ar1ma1_b.pdf", replace





