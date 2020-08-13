*********************************************************
* Ch6: Seasonal ARMA(p,q) processes
*********************************************************

cd "C:\Users\johnlevendis\Google Drive\papers\book\ARIMA\"
set scheme s1mono

* Generating deterministically seasonal data
	local s = 4

	drop _all
	clear all
	set obs 100
	gen time = _n
	gen d = mod(t,`s')	
	gen D1 = 0
	gen D2 = 0
	gen D3 = 0
	gen D4 = 0
	replace D1 = 1 if d==1
	replace D2 = 1 if d==2
	replace D3 = 1 if d==3
	replace D4 = 1 if d==4
	drop d
	tsset time	
		
	set seed 1234
	gen error = rnormal()
	label var e "errors"

	local beta1 = .5
	local betaD1 = 5
	local betaD2 = 10
	local betaD3 = -3
	local betaD4 = 2

*Figure 6.2
		gen X= error in 1
		replace X= `betaD1'*D1 + `betaD2'*D2 + `betaD3'*D3 + `betaD4'*D4 + error in 2/L
		tsline X

**Seasonal (additive, quarterly) AR(1)	
*		drop X
*		gen X = error in 1
*		replace X = `beta1'*L.X + `betaD1'*D1 + `betaD2'*D2 + `betaD3'*D3 + `betaD4'*D4 + error in 2/L
*		graph twoway (connected X time)	
*		tsline X
	

	

* Downloading seasonal data 

	freduse LNU04000001 LNU04000002 HOUSTNSA RETAILNS LNU01000000, clear
	
	rename LNU04000001 MALE_UR
	rename LNU04000002 FEMALE_UR
	rename LNU01000000 CIVILIAN_LF
	rename HOUSTNSA HOUSING_STARTS
	rename RETAILNS RETAIL_SALES
	
	replace RETAIL = RETAIL/1000
	replace CIVILIAN_LF = CIVILIAN_LF/1000
	
	label var MALE_UR "Male Unemployment Rate"
	label var FEMALE_UR "Female Unemployment Rate"
	label var HOUSING_STARTS "Housing Starts, thousands"
	label var RETAIL_SALES "Retail Sales, nominal, millions USD"
	label var CIVILIAN_LF "Civilian labor force, thousands"
	
	drop date
	gen date=ym(year(daten), month(daten))  
	tsset date, monthly  
	format date %tm
	drop daten
	order date RETAIL HOUSING CIVILIAN MALE FEMALE
	
	save ".\data\seasonal_data.dta", replace
	

*Figure 6.1
	use ".\data\seasonal_data.dta", clear
	gen lnRETAIL = ln(RETAIL_SALES)
	
	tsline RETAIL_SALES if year(date)>1960, xtitle("")  tlabel(,alternate)
	tsline MALE_UR if year(date)>1960, xtitle("")  tlabel(,alternate)
	tsline HOUSING_STARTS if year(date)>1960, xtitle("")  tlabel(,alternate)
	tsline CIVILIAN_LF if year(date)>1960, xtitle("") tlabel(,alternate)

	
* Seasonal differencing
	use ".\data\seasonal_unemp_rate.dta", clear
	gen temp = unemp - L4.unemp
	list in 1/10


* Figure 6.3
	use ".\data\seasonal_data.dta", clear 
	* Seasonal differencing:
		gen RETAIL_sd = ln(RETAIL) - ln(L12.RETAIL)
		gen HOUSING_sd = HOUSING - L12.HOUSING
		gen CIVILIAN_sd = CIVILIAN - L12.CIVILIAN
		gen MALE_sd = MALE - L12.MALE
	* Label the variables
		label var RETAIL_sd "ln(RetailSales) - ln(L12.RetailSales)"
		label var HOUSING_sd "HousingStarts - L12.HousingStarts"
		label var CIVILIAN_sd "LaborForce - L12.LaborForce"
		label var MALE_sd "UnempRate - L12.UnempRate"
	* Graphing
		tsline RETAIL_sd, xtitle("")  tlabel(,alternate)
		tsline HOUSING_sd, xtitle("")  tlabel(,alternate)
		tsline CIVILIAN_sd, xtitle("") tlabel(,alternate)
		tsline MALE_sd, xtitle("")  tlabel(,alternate) 
 
 
 
* Figure 6.4: ARIMA(4,0,0)
	set seed 3543
	set more off
	drop _all
	set obs 10000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local beta4 = 0.50
	replace X = `beta4'*L4.X + e in 5/L 
	arima X, ar(12) nocons 
	ac X, lags(30)
	pac X, lags(30)
		


* Figure 6.5: ARIMA(0,0,4)
	set seed 9876
	set more off
	drop _all
	set obs 5000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local gamma12 = 0.50
	replace X = `gamma12'*L12.e + e in 13/L 
	arima X, ma(12) nocons 
	ac X, lags(60) 
	pac X, lags(60) 



	
************************************
* Exercise (1)
************************************


	* ARIMA(1,0,0)x(2,0,0)_{12}
	set seed 8765
	set more off
	drop _all
	set obs 1000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local beta1  = 0.10
	local beta12 = 0.40
	local beta24 = 0.40
	replace X = `beta1'*L.X + `beta12'*L12.X -`beta1'*`beta12'*L13.X -`beta1'*`beta24'*L25.X + e in 26/L
	tsline X
	ac X, lags(60)
	pac X, lags(60)


	
	
/*	

* Seasonal white noise data 
	drop _all
	clear all
	set more off
	local s = 12  	/* enter the seasonality, s, here */
	set obs 100		/* enter the number of observations here */
	gen time = _n
	gen d = mod(t,`s')
	replace d = `s' if d == 0
	tabulate d, gen(D)
	drop d
	tsset time	
	set seed 1234
	gen double error = rnormal()
	label var e "errors"

	gen Dtot = 0
	forvalues i = 1/`s' {
		local betaD`i' = `s'/(`i'+1)
		display "betaD`i'" 
		replace D`i' = D`i' * `betaD`i''
		replace Dtot = Dtot + D`i'
	}
	gen X= error + Dtot
	graph twoway (connected X time)



* ARIMA(1 12,0,0)
	set seed 9393
	set more off
	drop _all
	set obs 5000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local beta1  = 0.10
	local beta12 = 0.50
	replace X = `beta1'*L.X + `beta12'*L12.X + e in 14/L 
	arima X, ar(1 12) nocons 
	ac X, lags(60)
	pac X, lags(60)



* ARIMA(1,0,0)x(1,0,0)_{12}
	set seed 8765
	set more off
	drop _all
	set obs 5000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local beta1  = 0.10
	local beta12 = 0.50
	replace X = `beta1'*L.X + `beta12'*L12.X -`beta1'*`beta12'*L13.X + e in 14/L
	arima X, arima(1,0,0) sarima(1,0,0,12) nocons 
	* Alternatively, we could estimate: 
	arima X, arima(1,0,0) mar(1,12) nocons 
	ac X, lags(60)
	pac X, lags(60)



* ARIMA(0,0,1)x(0,0,1)_{12}
	set seed 7654
	set more off
	drop _all
	set obs 5000
	gen time = _n
	tsset time
	gen e = rnormal()
	gen X = e
	local gamma1  = 0.10
	local gamma12 = 0.50
	replace X = `gamma1'*L.e + `gamma12'*L12.e -`gamma1'*`gamma12'*L13.e + e in 14/L
	arima X, arima(0,0,1) sarima(0,0,1,12) nocons 
	ac X, lags(60) 
	pac X, lags(60)
	
*/	

