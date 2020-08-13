
***************************************
* Ch5: Non-stationarity and ARIMA processes
***************************************

set scheme s1mono


* Figure 5.1a
	* Nominal GDP per capita
	*wbopendata , country(usa) indicator(ny.gdp.pcap.cd) year(1960:2015) long nometadata clear
	*rename ny_gdp_pcap_cd NomPCGDP
	*save ".\data\gdppcap.dta", replace
	use ".\data\gdppcap.dta", clear
	graph twoway connected NomPCGDP year, ytitle("") xtitle("")

* Figure 5.1b
	* Consumer Price Index
	*wbopendata , country(usa) indicator(fp.cpi.totl) year(1960:2015) long nometadata clear
	*rename fp CPI
	*save ".\data\CPI.dta", replace
	use ".\data\CPI.dta", clear
	graph twoway connected CPI year, ytitle("") xtitle("")

* Figure 5.2a
	* Dow-Jones Industrial Average
	*fetchyahooquotes ^DJI, freq(d) start(1jan1995) end(31dec2015) 
	*rename adjclose__DJI DJIA
	*label var DJIA "Dow-Jones Industrial Average"
	*format date %tdCCYY
	*save ".\data\DJIA.dta", replace
	use ".\data\DJIA.dta", clear
	graph twoway line DJIA date, ytitle("") xtitle("")

* Figure 5.2b
	* Google Shares
	*fetchyahooquotes GOOG, freq(d) start(1jan1995) end(31dec2015)
	*rename adjclose_GOOG GOOG
	*label var GOOG "Google shares"
	*format date %tdCCYY
	*save ".\data\GOOG.dta", replace
	use ".\data\GOOG.dta", clear
	graph twoway line GOOG date, ytitle("") xtitle("") 
	graph twoway line D1.GOOG date, ytitle("") xtitle("")
	
* Figure 5.3
	use ".\data\DJIA.dta", clear
	graph twoway line D1.DJIA date if year(date)>=2000 & year(date)<=2010, ytitle("") xtitle("")


	
	

************************************	
* 5.1 Differencing 
* Exercise 2
************************************		

	set more off
	drop _all
	clear all
	set obs 100
	set seed 23456
	gen time = _n
	tsset time
	gen double error = rnormal()
	
	gen Z = error	// WN: Z~I(0)
	gen Y = sum(Z)  // RW: Y~I(1)
	gen X = sum(Y)	// X~I(2)

	gen A = error
	replace A = 0.50*L.A + error in 2/L
	gen B = sum(A)
	gen C = sum(B)
	
	*save ".\data\integrated012.dta", replace
	
	foreach var of varlist A B C X Y Z {
		tsline `var'
		tsline D.`var'
		tsline D2.`var'
	}
	
		
	
	
	
	


************************************************
* 5.6 Differencing and detrending appropriately
* and
* Exercise
************************************************

* Create the data
drop _all
clear all

local beta0 = 1
local beta1 = 1
local observations = 100 /* Exercise: Adjust to see how obs to see how this affects the variance. */
set obs  `observations'
set seed 1234
gen time = _n
tsset time
gen error = rnormal()

* Generate a Difference Stationary process (random walk with drift)
	gen y = error in 1
	replace y = `beta0' + L.y + error in 2/L

* Generate a Trend Stationary process
	gen x = error in 1
	replace x = `beta1'*time + error in 2/L


sum d.x d.y

* What if we detrend both datasets
	quietly reg y time
	predict dty, resid
	quietly reg x time
	predict dtx, resid

sum d.x d.y dtx dty	/* See how obs affects variance */

* Figure 5.4
	set scheme s1mono
	twoway connected x time, msymbol(i)
	
* Figure 5.5
	twoway connected y time, msymbol(i)

* Figure 5.6
	twoway connected d.x time, msymbol(i) 

* Figure 5.7
	twoway connected d.y time, msymbol(i)

* Figure 5.8
	twoway connected dtx time, msymbol(i)
	
* Figure 5.9
	twoway connected dty time, msymbol(i)
	
* Figure 5.10
	ac dty
	pac dty

	


************************************
* 5.6 Mistakenly differencing/detrending
************************************

* Exercise
	drop _all
	set seed 1776

	set obs 1000
	gen time = _n
	tsset time

	local variance = 4
	gen X = rnormal(0,sqrt(`variance'))

	qui sum X
	display "E(X)=" r(mean) " Var(X)=" r(Var)
	qui sum D1.X
	display "E(X)=" r(mean) " Var(X)=" r(Var)

	corrgram X, lags(5)
	corrgram D1.X, lags(5)

	arima D1.X, ma(1)
	estat aroots

	
	

************************************************
* 5.7 Replicating Granger and Newbold (1974)
************************************************

set scheme s1mono

drop _all
set seed 1234
set obs 100

gen double e1 = rnormal()
gen double e2 = rnormal()
gen time = _n
tsset time

* Creating the first random walk 
gen double X1 = e1 in 1
replace X1= L.X1 + e1 in 2/L

* Creating the second random walk 
gen double X2 = e2 in 1
replace X2= L.X2 + e2 in 2/L

* Figure 5.11a
graph twoway (line X1 time) (line X2 time, lpattern(dash))

reg X1 X2
predict X1hat

* Figure 5.11b
graph twoway (scatter X1 X2) (line X1hat X2)




*****************************************
* Replicating Granger and Newbold (1974)
*****************************************

capture program drop GN
program GN, rclass
	version 12.0
	drop _all
	*set obs $numobs
	set obs 50
	
	* Creating two pure random walk processes
	gen double e1 = rnormal()
	gen double e2 = rnormal()
	gen time = _n
	tsset time

	gen double X1 = e1 in 1
	replace X1= L.X1 + e1 in 2/L

	gen double X2 = e2 in 1
	replace X2= L.X2 + e2 in 2/L

	reg X1 X2
	estat dwatson

	return scalar DWstat = r(dw)
	return scalar R2 = e(r2)
	return scalar Pval = (2 * ttail(e(df_r), abs(_b[X2]/_se[X2])))

end

set seed 1234
simulate DWstat = r(DWstat) R2=r(R2) Pval=r(Pval), reps(200) saving(GNtemp2, replace) nodots nolegend: GN
summarize
count if Pval < 0.05











