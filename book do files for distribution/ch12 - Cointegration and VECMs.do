********************************************************************
* Cointegration and VECMs
********************************************************************

set scheme s1mono

*************************************
* Illustration of cointegrated variables
*************************************

clear all
graph drop _all
cd ".\graphs\"

set obs 50
set seed 4321
gen time = _n
tsset time
label var time "Time"

* Draw some random errors
	gen ex = rnormal(0,2)
	gen ey = rnormal(0,2)
	gen ez = rnormal(0,2)

* Generate the variables
	gen X = ex in 1
	replace X = 1 + L.X + ex in 2/L
	gen Y = 10 + X + ey
	gen Z = 10 + 3*X + ez


* Figure 12.1
	graph twoway  (connected Y time, msymbol(none) xtitle("")) (connected X time, msymbol(none) lpattern(dash)) , ///
		legend(ring(0) bplacement(seast) cols(1)) ///
		saving(coint01.gph, replace) 
	gen diff1 = Y - X
	label var diff1 "Y - X"
	graph twoway (connected diff1 time, msymbol(none) ytitle("")) , ///
		legend(on ring(0) bplacement(seast) label(1 "= Y - X")) ///
		saving(coint02.gph, replace)  
	graph combine coint01.gph coint02.gph, col(1)
		graph export coint0102.pdf, replace

		
		
* Figure 12.2
	graph twoway (connected Z time, msymbol(none) xtitle("")) (connected X time, msymbol(none) lpattern(dash)) , ///
		legend(ring(0) bplacement(seast)) ///
		saving(coint03.gph, replace) 
	gen diff2 = Z - X
	graph twoway (connected diff2 time, msymbol(none) ytitle("")) , ///
		legend(on ring(0) bplacement(seast) label(1 "= Z - X")) ///
		saving(coint04.gph, replace)
	graph combine coint03.gph coint04.gph, col(1)
		graph export coint0304.pdf, replace

* Figure 12.3
* We can tilt X up
	gen Xb = 3*X
	graph twoway (connected Xb time, msymbol(none) lpattern(dash)) (connected Z time, msymbol(none)) , ///
		legend(on ring(0) bplacement(seast) label(1 "X' = 3X")) xtitle("") ///
		saving(coint05.gph, replace)
	gen diff3 = Xb - Z
	graph twoway (connected diff3 time, msymbol(none)) , ///
		legend(on ring(0) bplacement(seast) label(1 "= X' - Z")) ytitle("") ///
		saving(coint06.gph, replace)
	graph combine coint05.gph coint06.gph, col(1)
		graph export coint0506.pdf, replace

* Figure 12.4
* Or we can tilt Z down
	gen Zb = Z/3
	graph twoway (connected Zb time, msymbol(none)) (connected X time, msymbol(none) lpattern(dash)) , ///
		legend(on ring(0) bplacement(seast) label(1 "Z'= Z/3")) xtitle("") ///
		saving(coint07.gph, replace)
	gen diff4 = X - Zb
	graph twoway (connected diff4 time, msymbol(none)) , ///
		legend(on ring(0) bplacement(seast) label(1 "= X - Z'")) ///
		saving(coint08.gph, replace)
	graph combine coint07.gph coint08.gph, col(1)
		graph export coint0708.pdf, replace

graph drop _all


* Figure 12.5
	*freduse PCEC GDP, clear
		*save "PCEC_GDP.dta", replace
		use  "..\PCEC_GDP.dta", clear
		* GDP  = "GDP, billions, quarterly, seaonally adjusted annual rate"
		* PCEC = "Personal Consumption Expenditures, quarterly, seasonally adjusted"
		gen lnGDP = ln(GDP)
		gen lnC = ln(PCEC)
		label var lnGDP "ln(GDP)"
		label var lnC "ln(personal consumption)"
	*Time
		gen time = mofd(daten)
		format time %tm
		gen year = year(daten)
		tsset time
	graph twoway (connected lnGDP time, msymbol(none)) ///
		(connected lnC time, symbol(none) lpattern(dash)) , ///
		legend(on ring(0) bplacement(seast) cols(1)) ///
		xtitle("") 
	graph export coint_ex_1.pdf, replace



* Figure 12.6	
	freduse AAA BAA GS10 GS5 , clear
		* AAA = Moody's Seasoned Aaa Corporate Bond Yield, monthly, not seasonally adjusted
		* Baa = Moody's Seasoned Baa Corporate Bond Yield, monthly, not seasonally adjusted
		* GS10 = 10-Year Treasury Constant Maturity Rate, monthly, pct, not seasonally adjusted
		* GS5 = 5-Year Treasury Constant Maturity Rate, pct, not seasonally adjusted 
	* Time
		gen time = mofd(daten)
		format time %tm
		gen year = year(daten)
		tsset time
		keep if year > 1980 & year < 2016

	graph twoway (connected BAA time, msymbol(none)) ///
		(connected AAA time, symbol(none) lpattern(dash)) , ///
		legend(on ring(0) bplacement(neast) cols(1) ///
		label(1 "Moody's Baa corporate bond yield") ///
		label(2 "Moody's Aaa corporate bond yield")) ///
		xtitle("") ytitle("Percent yield") 
	graph export coint_ex_2.pdf, replace

* Figure 12.7	
	graph twoway (connected GS10 time, msymbol(none)) ///
		(connected GS5 time, symbol(none) lwidth(medthick)) , ///
		legend(ring(0) bplacement(neast) cols(1)) ///
		ytitle("Percent yield") xtitle("") 
	graph export coint_ex_3.pdf, replace


	
	
	
******************************************************************
* 12.5 Engle and Granger's residual-based tests of cointegration
******************************************************************

****************************************
* 12.5.2 Engle-Granger approach
****************************************

clear all
set obs 50
set seed 4321
gen time = _n
tsset time

* Draw some random errors
	gen ex = rnormal(0,2)
	gen ey = rnormal(0,2)

* Generate the variables, X and Y
	gen X = ex in 1
	replace X = 1.0 + L.X + ex in 2/L
	gen Y = 10 + X + ey

* Are X and Y I(1)?
	dfuller X, drift
	dfuller Y, drift
	dfuller D.X
	dfuller D.Y

* Estimate the long-run relationship, and extract the residuals
	reg Y X
	predict ehat, residuals

* Verify that the residual is I(0), ie verify that Y and X are cointegrated
* We can get our test statistics from
	dfuller ehat, nocons
	*or
	reg d.ehat L.ehat, nocons
* which would give us the right test statistics, bu the wrong critical values. 
* We need to use MacKinnon critical values.
	
* MacKinnon (2010) critical values
* Using the tables, two variables, n=49, p=0.05, and a constant:
	di -3.33613 + (-6.1101)/49 + (-6.823)/(49^2)
	*or
	egranger Y X, regress 

* Estimate the ECM
	reg D.Y D.X L.ehat



	
********************************************************************
* 12.6.4 Johansen's tests and the rank of PI
********************************************************************

set scheme s1mono

cls
clear all
set more off
set obs 1000
set seed 4321
gen time = _n
tsset time
label var time "Time"

* Draw some random errors
	gen ex = rnormal(0,1)
	gen ey = rnormal(0,1)
	gen ez = rnormal(0,1)
	gen ev = rnormal(0,1)

* Generate the variables
	gen X = ex in 1
	replace X = 1 + L.X + ex in 2/L
	gen Y = 10 + 2*X + ey
	gen Z =  5 + 3*X + ez
	gen V = ev in 1
	replace V = 2 + L.V + ev in 2/L

* Vecrank
	vecrank X Y Z V, trend(constant) max	

* The next step is to estimate these long-run cointegrating equations,
* as well as the short-run adjustment terms.	
vec Z Y X , rank(2) trend(constant) 	

	

matrix list e(alpha)
matrix list e(beta) // which, in my notation is B' //
matrix list e(pi)
matrix list e(b)

* Figure 12.8
twoway (connected X time, msymbol(none) lpattern(solid) lwidth(thick)) ///
	   (connected Z time, msymbol(none) lpattern(dash)  lwidth(medthick)) ///
	   (connected Y time, msymbol(none) lpattern(solid) lwidth(medium)) ///
	   (connected V time, msymbol(none) lpattern(dash)  lwidth(thin)) ///
	   if time<=100 , legend(cols(4))

	
	
	


 






****************************************
* 12.11 Exercises
****************************************

* Exercises (4) and (5)
	* Set up
		clear all
		set more off
		set obs 500
		set seed 54321
		gen time = _n
		tsset time

	* Draw some random errors
		gen ex = rnormal(0,1)
		gen ey = rnormal(0,2)
		gen ez = rnormal(0,3)

	* Generate the variables, X and Y
		gen X = ex in 1
		replace X = 1.1 + L.X + ex in 2/L
		gen Y = X + ey
		gen Z = 2 + X + ez

	* Estimate the XY model.
		vec X Y, rank(1) trend(constant) alpha
		* The alphas look wrong, don't they?

	* Estimate the XZ model.
		vec X Z, rank(1) trend(constant)
		* The alphas look wrong, don't they?

	* (X) = (1.1672) + (-0.1227 0.0199)(DL.X) + (0.045)(L.X + 1.27 - 1.00 L.Y) + (e1)
	* (Y)   (-0.046)   (-0.4256 0.1263)(DL.Y)   (1.128)                          (e2)
		
	* (X) = (1.2449 ) + (-0.07975 -0.021304)(DL.X) + (-0.0164)(L.X + 3.2727 - 1.00 L.Z) + (e1)
	* (Z)   (0.01908)   (-0.22100 0.0566547)(DL.Z)   (1.07436)                            (e2)




