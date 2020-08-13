

*****************************************************
* Chapter 11: Vector Autoregressions II: Extensions
*****************************************************
	
	
**********************************
* 11.1 Orthogonalized IRFs
**********************************

set more off
clear all

* Generate data
	local obs = 10000
	set obs `obs'
	gen time = _n
	tsset time

	matrix sigma = (1 , .75 \ .75, 1)
	drawnorm eps_x eps_y, means(0,0) cov(sigma) double seed(1234)

	gen X = 0
	gen Y = 0

	forvalues i = 2/`obs'{
		qui replace X = 0.40*L.X + 0.10*L.Y + eps_x in `i'
		qui replace Y = 0.20*L.X + 0.10*L.Y + eps_y in `i'
	}
	*

qui var X Y, lags(1)
capture irf drop order1
capture irf drop order2
capture irf drop varYX 
capture iff drop varXY
irf create varXY, step(10) set(myirf1) replace
irf table irf, noci
irf table oirf, noci





matrix sigmahat = e(Sigma)
qui varstable, amat(betahat)

matrix L = cholesky(sigmahat)
matrix Lt = L'
*matrix LtInv = inv(L')
matrix LInv = inv(L)
*matrix list LInv

* Constructing the IRFs and OIRFs.
qui var X Y, lags(1)
	irf table irf, noci
	irf table oirf, noci
	
	matrix betahat0 = I(2)
	matrix betahat1 = betahat
	matrix betahat2 = betahat*betahat
	matrix betahat3 = betahat*betahat*betahat
	matrix betahat4 = betahat*betahat*betahat*betahat
	
	matrix phi_0 = L
	matrix phi_1 = betahat1*L
	matrix phi_2 = betahat2*L
	matrix phi_3 = betahat3*L
	matrix phi_4 = betahat4*L
	
	* IRFs
	matrix list betahat0
	matrix list betahat1
	matrix list betahat2
	matrix list betahat3
	matrix list betahat4
			
	* OIRFs
	matrix list phi_0
	matrix list phi_1 
	matrix list phi_2
	matrix list phi_3 
	matrix list phi_4 

	
*************************
* Using Cholesky to Orthogonalize Variables
*************************
	
* Showing that the transformed errors are now orthogonalized:
	drawnorm eps1 eps2, means(0,0) cov(sigma) double seed(1234)
	* Verify that they have the appropriate variance/covariance structure:
		correlate eps1 eps2, covariance
	* Now, multiply the error matrix by the inverse of the lower Cholesky factor:
		gen e_x = LInv[1,1]*eps_x + LInv[1,2]*eps_y
		gen e_y = LInv[2,1]*eps_x + LInv[2,2]*eps_y
	* Verify that the transformed errors have the appropriate variance/covariance structure:
		correlate e_x e_y, covariance
	* Graphs
		*twoway scatter eps_x eps_y, msymbol(point)
			*graph export choleskytrans_01.eps, replace
			*graph export choleskytrans_01.pdf, replace
		*twoway scatter e_x e_y, msymbol(point)
			*graph export choleskytrans_02.eps, replace
			*graph export choleskytrans_02.pdf, replace
	
	
	*****************************
	* Showing that order matters:
	*****************************
	
	qui var Y X, lags(1)
	capture irf drop order1
	capture irf drop order2
	capture irf drop varYX 
	capture iff drop varXY
	irf create varYX, step(10) set(myirf2) replace
	irf table irf, noci
	irf table oirf, noci

	matrix sigmahat = e(Sigma)
	qui varstable, amat(betahat)

	matrix L = cholesky(sigmahat)
	matrix Lt = L'
	*matrix LtInv = inv(L')
	matrix LInv = inv(L)

	* Constructing the IRFs and OIRFs.
	qui var Y X, lags(1)
		irf table irf, noci
		irf table oirf, noci

		matrix betahat0 = I(2)
		matrix betahat1 = betahat
		matrix betahat2 = betahat*betahat
		matrix betahat3 = betahat*betahat*betahat
		matrix betahat4 = betahat*betahat*betahat*betahat
		
		matrix phi_0 = L
		matrix phi_1 = betahat1*L
		matrix phi_2 = betahat2*L
		matrix phi_3 = betahat3*L
		matrix phi_4 = betahat4*L
		
		* IRFs
		matrix list betahat0
		matrix list betahat1
		matrix list betahat2
		matrix list betahat3
		matrix list betahat4
				
		* OIRFs
		matrix list phi_0
		matrix list phi_1 
		matrix list phi_2
		matrix list phi_3 
		matrix list phi_4 

		

		
		


**********************************
* 11.1 Orthogonalized IRFs
**********************************

* VAR(1)
clear all
set more off
local obs = 10000
set obs `obs'
gen time = _n-1
tsset time

matrix sigma = (1,.75 \.75, 1)
drawnorm e1 e2, means(0,0) cov(sigma) double seed(1234)

gen Y = 0
gen X = 0

quietly{
	forvalues i = 2/`obs'{
		replace Y =  0.30*L.Y - 0.50*L.X + e1 in `i'
		replace X = -0.40*L.Y + 0.30*L.X + e2 in `i'
	}
}
save "VAR_sim.dta", replace



**********************************
* Figure 11.1
**********************************
use "VAR_sim.dta", clear
var Y X , lags(1) nocons 
irf create irf , set(name1, replace)
irf graph irf oirf, noci plot2opts(lpattern(dash))
	*graph export irf_vs_oirf.eps, replace
	*graph export irf_vs_oirf.pdf, replace

*predict errorY, residual equation(Y)
*predict errorX, residual equation(X)
*correlate errorY errorX, covariance


**********************************
* Figure 11.1
**********************************
clear all
use "VAR_sim.dta", clear
var Y X , lags(1) nocons 

irf create order1, step(10) order(Y X) set(myirf1, replace) 
irf create order2, step(10) order(X Y)

irf ograph	(order1 X X oirf, noci) (order2 X X oirf, noci clpat(dash)) , saving(VAR_oirf_12_xx, replace)
irf ograph 	(order1 X Y oirf, noci) (order2 X Y oirf, noci clpat(dash)) , saving(VAR_oirf_12_xy, replace) 
irf ograph	(order1 Y X oirf, noci) (order2 Y X oirf, noci clpat(dash)) , saving(VAR_oirf_12_yx, replace)
irf ograph	(order1 Y Y oirf, noci) (order2 Y Y oirf, noci clpat(dash)) , saving(VAR_oirf_12_yy, replace)

graph combine VAR_oirf_12_yy.gph ///
	VAR_oirf_12_xy.gph ///
	VAR_oirf_12_yx.gph ///
	VAR_oirf_12_xx.gph, ///
	/*iscale(1)*/ ycommon saving(VAR_oirf_12, replace)
	*graph export VAR_oirf_12.eps, replace		
	*graph export VAR_oirf_12.pdf, replace		
			

			
			



	


*****************************
* 11.1.2 Cholesky decompositions
* Analytics
*****************************
			
matrix sigma = (1 , .75 \ .75, 1)
matrix list sigma
* A matrix is positive definite IFF all of its eigenvalues are positive
	matrix eigenvalues real immaginary = sigma 
	matrix list real
	matrix list immaginary
* Cholesky
	matrix L = cholesky(sigma)
* Verify
	matrix list L
	matrix Lt = L'
	matrix list Lt
	matrix temp = L*Lt
	matrix list temp

	
	
matrix beta = (1 , .30 , 0.20 \ .30, 1, .10 \ .20 , .10, 1)
matrix list beta
* A matrix is positive definite IFF all of its eigenvalues are positive
	matrix eigenvalues real immaginary = sigma 
	matrix list real
	matrix list immaginary
* Cholesky
	matrix L = cholesky(beta)
* Verify
	matrix list L
	matrix Lt = L'
	matrix list Lt
	matrix temp = L*Lt
	matrix list temp


	
	
	
	

*********************************************
* 11.1.2 Using Cholesky to Orthogonalize Variables
* Simulations
*********************************************

clear all

matrix sigma = (1 , .75 \ .75, 1)
matrix L = cholesky(sigma)  /* The lower triangular factor */
matrix Lt = L'              /* L-transpose, the upper factor */
matrix LInv = inv(L)        /* The inverse of L  */

* Generate the untransformed variables:
	set obs 10000
	gen time = _n-1
	tsset time
	drawnorm eps1 eps2, means(0,0) cov(sigma) double seed(1234)


* Verify that they have the appropriate variance/covariance structure:
	correlate eps1 eps2, covariance

* Now, multiply the error matrix by the inverse of the lower Cholesky factor:
	gen e1 = eps1*LInv[1,1] + eps2*LInv[1,2]
	gen e2 = eps1*LInv[2,1] + eps2*LInv[2,2]

* Verify that the transformed errors have the appropriate variance/covariance structure:
	correlate e1 e2, covariance

* Graphs
	twoway scatter eps1 eps2, msymbol(Oh)
		*graph export choleskytrans_01.eps, replace
		graph export choleskytrans_01.pdf, replace
	twoway scatter e1 e2, msymbol(Oh)
		*graph export choleskytrans_02.eps, replace
		graph export choleskytrans_02.pdf, replace

		




*****************************************
* 11.1.2 Expressing OIRFs as the 
* coefficients of the MA representation
*****************************************

set more off
clear all

* Generate data
	local obs = 10000
	set obs `obs'
	gen time = _n
	tsset time

	matrix sigma = (1 , .75 \ .75, 1)
	drawnorm eps_x eps_y, means(0,0) cov(sigma) double seed(1234)

	gen X = 0
	gen Y = 0

	forvalues i = 2/`obs'{
		qui replace X = 0.40*L.X + 0.10*L.Y + eps_x in `i'
		qui replace Y = 0.20*L.X + 0.10*L.Y + eps_y in `i'
	}
	*

* OIRF via Cholesky
	qui var X Y, lags(1)
	* OIRF
	capture irf drop order1
	capture irf drop order2
	capture irf drop orderYX 
	capture iff drop orderXY
	irf create orderXY, step(10) set(myirf, replace) replace
	irf table oirf, noci

* OIRF via SVAR
	matrix A = (1, 0 \ ., 1) 
	matrix B = (., 0 \ 0, .)
	qui svar X Y, aeq(A) beq(B) lags(1)
	* OIRF
	capture irf drop order1
	capture irf drop order2
	capture irf drop orderYX 
	capture iff drop orderXY
	irf create orderXY, step(10) set(myirf, replace) replace
	irf table oirf, noci


	

	
	
















	
	




************************************
* Rough replication of Sims (1980)
************************************

drop _all
freduse GNPC96 MANMM101USQ189S GNPDEF A576RC1Q027SBEA /// 
	B021RG3Q086SBEA USAURHARMQDSMEI
gen time=yq(year(daten), quarter(daten))  
tsset time, quarterly 
drop if year(time)<1960
drop date 

label var GNPC9 "Real GNP, base 2009"
label var MANMM "M1 for the US"
label var GNPDE "GNP Implicit Price Deflator"
label var A576R "Compensation of employees: Wages and salaries"
label var B021R "Imports of goods and services (chained price index)" 
label var USAUR "Harmonized Unemployment Rate: All Persons for the US"
label var USAUR "Unemployment Rate"

gen Ygr = ln(GNPC9) - ln(L.GNPC9)
gen Mgr = ln(MANMM) - ln(L.MANMM)
gen Pgr = ln(GNPDE) - ln(L.GNPDE)
gen Wgr = ln(A57) - ln(L.A57)
gen IPgr= ln(B02) - ln(L.B02)
gen Ugr = ln(USAUR) - ln(L.USAUR)

label var Ygr "Real GNP growth rate"
label var Mgr "M1 growth rate"
label var Pgr "Inflation rate"
label var IPgr "Import price growth rate"
label var Wgr "Wages and Salary growth rate"
label var Ugr "Unemployment growth rate"

* You may need to install the KPSS command
foreach var of varlist Ygr Mgr Pgr Wgr IPgr Ugr{
	kpss `var' if year(daten)< 1980
	tsline `var' if year(daten)< 1980
}

rename USAUR Unemp
rename Pgr Infl

twoway (connected Unemp time, msymbol(none) lpattern(longdash) lwidth(medthick) ) (connected Infl time, yaxis(2) msymbol(none)  )
*	graph export ts_phillips.pdf
*	graph export ts_phillips.eps

varsoc Infl Unemp
irf set, clear

var Infl Unemp, lag(1 2 3)
varstable
varstable , graph modlabel
varlmar

irf set temp, replace
irf create order1, step(10) order(Infl Unemp)
irf create order2, step(10) order(Unemp Infl)
irf graph oirf, irf(order1) noci
*	graph export irf_phillips_1.pdf
*	graph export irf_phillips_1.eps
irf graph oirf, irf(order2) noci
*	graph export irf_phillips_2.pdf
*	graph export irf_phillips_2.eps

