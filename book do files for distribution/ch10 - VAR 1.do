*****************************************************
* Chapter 10: Vector Autoregressions I: Basics
*****************************************************
	
set scheme s1mono



	

******************************
* 10.2 A Simple VAR(1) and how to estimate it
******************************

drop _all
freduse GNPC96 MANMM101USQ189S
gen time=yq(year(daten), quarter(daten))  
tsset time, quarterly 
gen GNPgr = ln(GNPC96) - ln(L.GNPC96)
gen M1gr  = ln(MANMM)  - ln(L.MANMM)
label var GNPgr "Growth rate of real GNP"
label var M1gr "Growth rate of M1"
drop GNPC MANMM date daten
drop if year(time)<1960
tsline GNPgr M1gr, lwidth(thick thin)
	*graph export VARdata01.eps, replace
	graph export VARdata01.pdf, replace
reg GNPgr L.GNPgr L.M1gr 		if time<yq(2017,3)
reg M1gr L.GNPgr L.M1gr			if time<yq(2017,3)
sureg (GNPgr L.GNPgr L.M1gr) (M1gr L.GNPgr L.M1gr) 
var GNPgr M1gr if time<yq(2017,3), lags(1)			
varbasic GNPgr M1gr if time<yq(2017,3), lags(1)
*irf graph irf, saving("IRF01")
irf graph irf
	*graph export IRF01.eps, replace
	*graph export IRF01.pdf, replace
varsoc GNPgr M1gr if time<yq(2017,3), maxlag(8)

	
	


********************************************************************
* 10.5 Stability
********************************************************************


cd "C:\Users\johnlevendis\Google Drive\papers\book\VAR\graphs\"
clear all
set more off
local obs = 1000000
set obs `obs'
gen time = _n-1
tsset time

matrix sigma = (1, 0 \ 0, 1)
drawnorm eX eY , means(0,0) cov(sigma) double seed(1234)

gen X = 0
gen Y = 0

quietly{
	forvalues i = 2/`obs'{
		replace X =  0.50*L.X + 0.60*L.Y + eX in `i'
		replace Y = -0.50*L.X + 0.50*L.Y +  eY in `i'
	}
}
var X Y, lags(1) nocons
estimates store Estimates
varstable , amat(A)
*varstable, graph modlabel
matrix list A

varstable, estimates(Estimates)


* Verify by hand
matrix beta = (0.50 , 0.60 \ -0.50 , 0.50)
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
display (real[1,1]^2 + immaginary[1,1]^2)^(1/2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)

varstable, estimates(Estimates) graph modlabel 


	

******************************
* 10.5 Exercises
******************************

* Exercise 1

matrix beta = (.1,-.2 \ -.3,.4)
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
*Since Stata orders eigenvalues from smallest to largest modulus, then
* we can just see whether the largest one is greater than or less than one.
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
local z1 =  ((beta[2,2] + beta[1,1]) + sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
local z2 =  ((beta[2,2] + beta[1,1]) - sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
di "z1 = " `z1'
di "z2 = " `z2'
di "1/z1 = " 1/`z1'
di "1/z2 = " 1/`z2'


matrix beta = (.2,.4 \ .4,.5)
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
local z1 =  ((beta[2,2] + beta[1,1]) + sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
local z2 =  ((beta[2,2] + beta[1,1]) - sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
di "z1 = " `z1'
di "z2 = " `z2'
di "1/z1 = " 1/`z1'
di "1/z2 = " 1/`z2'

matrix beta = (2,4 \ 4,5)
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
local z1 =  ((beta[2,2] + beta[1,1]) + sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
local z2 =  ((beta[2,2] + beta[1,1]) - sqrt( (beta[2,2] + beta[1,1])^2 - 4*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]) ))/(2*(beta[1,1]*beta[2,2] - beta[1,2]*beta[2,1]))
di "z1 = " `z1'
di "z2 = " `z2'
di "1/z1 = " 1/`z1'
di "1/z2 = " 1/`z2'

* Exercise 2

matrix beta = (-0.2,0.4,0.35 \ 0.3,-0.2,0.15\ 0.2,-0.3,0.4 )
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
display sqrt(real[1,3]^2 + immaginary[1,3]^2)

matrix beta = (-0.1,0.3,0.4 \ 0.9,-0.2,0.1\ 0.2,0.3,0.8 )
matrix list beta
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
display sqrt(real[1,3]^2 + immaginary[1,3]^2)


	
	
	
	



**********************************
* 10.8 Impulse response functions
**********************************

* Suppose:
* X_t = 0.40 X_{t-1} +  0.10 Y_{t-1}  + \epsilon_{x,t}
* Y_t = 0.20 X_{t-1} -  0.10 Y_{t-1}  + \epsilon_{y,t}

* Impulse from X
	drop _all
	set obs 11
	gen t = _n-1
	tsset t
	gen epsilonX = 0
	gen epsilonY = 0
	gen X = 0
	gen Y = 0
	label var X "X{subscript:t}"
	label var Y "Y{subscript:t}"
	replace epsilonX = 1 if t == 1
	forvalues i = 2/11{
	replace X = 0.40*L.X + 0.10*L.Y + epsilonX in `i'
	replace Y = 0.20*L.X - 0.10*L.Y + epsilonY in `i'
	}
	twoway connected X t, saving(IRFbyhandXtoX.gph, replace) note("Impulse=X; Response=X", pos(2) ring(0) size(small))
		*graph export IRFbyhandXtoX.pdf, replace
		*graph export IRFbyhandXtoX.eps, replace
	twoway connected Y t, saving(IRFbyhandXtoY, replace) note("Impulse=X; Response=Y", pos(2) ring(0) size(small))
		*graph export IRFbyhandXtoY.pdf, replace
		*graph export IRFbyhandXtoY.eps, replace

* Impulse from Y
	drop _all
	set obs 11
	gen t = _n-1
	tsset t
	gen epsilonX = 0
	gen epsilonY = 0
	gen X = 0
	gen Y = 0
	label var X "X{subscript:t}"
	label var Y "Y{subscript:t}"
	replace epsilonY = 1 if t == 1
	forvalues i = 2/11{
	replace Y = 0.20*L.X - 0.10*L.Y + epsilonY in `i'
	replace X = 0.40*L.X + 0.10*L.Y + epsilonX in `i'
	}
	twoway connected X t, saving(IRFbyhandYtoX, replace) note("Impulse=Y; Response=Y", pos(2) ring(0) size(small))
		*graph export IRFbyhandYtoX.pdf, replace
		*graph export IRFbyhandYtoX.eps, replace
	twoway connected Y t, saving(IRFbyhandYtoY, replace) note("Impulse=Y; Response=Y", pos(2) ring(0) size(small))
		*graph export IRFbyhandYtoY.pdf, replace
		*graph export IRFbyhandYtoY.eps, replace

* Combining the graphs
*	graph combine IRFbyhandXtoX.gph IRFbyhandXtoY.gph ///
*		IRFbyhandYtoX.gph IRFbyhandYtoY.gph, iscale(1) ycommon
*	graph export IRFbyhand.pdf, replace
*	graph export IRFbyhand.eps, replace

* MA representation of the IRF
* X_t = 0.40 X_{t-1} +  0.10 Y_{t-1}  + \epsilon_{x,t}
* Y_t = 0.20 X_{t-1} -  0.10 Y_{t-1}  + \epsilon_{y,t}
* So the companion matrix \beta is 
matrix beta = (0.40, 0.10 \ 0.20,-0.10)
matrix beta2 = beta*beta
matrix beta3 = beta*beta*beta
matrix list beta
matrix list beta2
matrix list beta3


	
	
	

	
******************************************************	
* 10.8.1 IRFs as the components of the MA coefficients
******************************************************

**********************************
* Example 2: Two-variable VAR(1)
**********************************
* Suppose:
* X_t = 0.50 X_{t-1} + 0.20 Y_{t-1}  + \epsilon_{x,t}
* Y_t = 0.30 X_{t-1} + 0.15 Y_{t-1}  + \epsilon_{y,t}

* So the companion matrix \beta is 
matrix beta = (0.50, 0.20 \ 0.30, 0.15)
matrix list beta
* Is this VAR stable?
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
* The lengths of the eigenvalues are 
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
* IRF
matrix beta2 = beta*beta
matrix beta3 = beta*beta*beta
matrix beta4 = beta*beta*beta*beta
matrix list beta
matrix list beta2
matrix list beta3
matrix list beta4

**********************************
* Example 3: Three-variable VAR(1)
**********************************
* Suppose:
* X_t = 0.25 X_{t-1} + 0.20 Y_{t-1}  + 0.15 Z_{t-1}  + \epsilon_{x,t}
* Y_t = 0.15 X_{t-1} + 0.30 Y_{t-1}  + 0.10 Z_{t-1} + \epsilon_{y,t}
* Z_t = 0.20 X_{t-1} + 0.25 Y_{t-1}  + 0.35 Z_{t-1} + \epsilon_{z,t}

* So the companion matrix \beta is 
matrix beta = (0.25, 0.20, 0.15 \ 0.15, 0.30, 0.10 \ 0.20, 0.25, 0.35)
matrix list beta
* Is this VAR stable?
matrix eigenvalues real immaginary = beta
matrix list real
matrix list immaginary
* The lengths of the eigenvalues are 
display sqrt(real[1,1]^2 + immaginary[1,1]^2)
display sqrt(real[1,2]^2 + immaginary[1,2]^2)
* IRF
matrix beta2 = beta*beta
matrix beta3 = beta*beta*beta
matrix beta4 = beta*beta*beta*beta
matrix list beta
matrix list beta2
matrix list beta3
matrix list beta4

**********************************
* Example 4: Two-variable VAR(2)
**********************************
* Suppose:
* X_t = 0.50 X_{t-1} + 0.20 Y_{t-1} + 0.10 X_{t-2} + 0.10 Y_{t-2} + \epsilon_{x,t}
* Y_t = 0.30 X_{t-1} + 0.15 Y_{t-1} + 0.20 X_{t-2} - 0.10 Y_{t-2} + \epsilon_{y,t}

clear all
set obs 10
gen t = _n-2
tsset t
gen X = 0
gen Y = 0
gen e_x = 0
gen e_y = 0
replace e_x = 1 if t==1
forvalues i=1/8{
replace X = 0.50*L.X + 0.20*L.Y + 0.10*L2.X + 0.10*L2.Y + e_x if t == `i'
replace Y = 0.30*L.X + 0.15*L.Y + 0.20*L2.X - 0.10*L2.Y + e_y if t == `i'
}
* Compare these values I just generated, with what I generate below using matrices.
* The companion matrix \beta is 
	matrix beta = (0.50, 0.20, 0.10, 0.10 \ 0.30, 0.15, 0.20, -0.10 \ 1, 0, 0, 0 \ 0, 1, 0, 0)
	matrix list beta
* IRF
	matrix beta2 = beta*beta
	matrix beta3 = beta*beta*beta
	matrix beta4 = beta*beta*beta*beta
	matrix list beta
	matrix list beta2
	matrix list beta3
	matrix list beta4

	




*******************************
* 10.9 Forecasting
*******************************

* generating fake time series for illustration purposes

	clear all
	local obs = 100
	set obs `obs'
	set seed 1234
	gen t = _n
	tsset t
	set more off
	
* Generate VAR(2) data

	gen double e1 = rnormal()
	gen double e2 = rnormal()

	gen double X = . 
	gen double Y = .

	qui replace X = e1 in 1/2
	qui replace Y = e2 in 1/2
	
	forvalues i = 3/`obs'{
		qui replace X = 0.30*L1.X + 0.20*L2.X + .10*L1.Y + .05*L2.Y + e1 in `i'
		qui replace Y = 0.35*L1.X + 0.25*L2.X + .15*L1.Y + .01*L2.Y + e2 in `i'
	}
	drop e1 e2
	var X Y, lags(1/2) nocons
	
* Estimate the VAR(2) so that it has exactly the coefficients I want it to have.
	constraint define 1 [X]L1.X = 0.30
	constraint define 2 [X]L2.X = 0.20
	constraint define 3 [X]L1.Y = 0.10
	constraint define 4 [X]L2.Y = 0.05
	constraint define 5 [Y]L1.X = 0.35
	constraint define 6 [Y]L2.X = 0.25
	constraint define 7 [Y]L1.Y = 0.15
	constraint define 8 [Y]L2.Y = 0.01
		
	var X Y, lags(1/2) nocons constraint(1 2 3 4 5 6 7 8)
	
* Add some empty observations
	tsappend, add(5)
	list t X Y

* Now forecast
	fcast compute E, nose step(5)
	list t X Y EX EY
	

	

	
	
	
**********************************
* 10.10 Granger causality
**********************************
clear all
freduse GNPDEF USAURHARMQDSMEI

gen time=yq(year(daten), quarter(daten))  
tsset time, quarterly 
drop if year(daten)<1960 | year(daten)>=1980

rename USAUR Unemp
label var Unemp "Unemployment Rate"

label var GNPDE "GNP Implicit Price Deflator"
gen Infl = ln(GNPDE) - ln(L.GNPDE)
label var Infl "Inflation rate"

var Infl Unemp if year(daten)< 1980, lag(1 2)
vargranger 
test [Infl] L1.Unemp [Infl] L2.Unemp
test [Unemp] L1.Infl [Unemp] L2.Infl	
		


******************************
* 10.10.1 Replicating Sims (1972)
******************************
drop _all

* Load the initial data
	set more off
	*freduse GNP AMBSL DISAMBSL CURRSL DISAMBNS , clear
	freduse GNP AMBSL, clear
	drop if year(daten)< 1947 | year(daten)> 1969  

	
* Create quarterly versions of these variables
	gen year = year(daten)
	gen quarter = quarter(daten)
	*collapse (mean) GNP AMBSL DISAMBSL CURRSL DISAMBNS, by(year quarter)
	collapse (mean) GNP AMBSL, by(year quarter)

	
* Create time variables
	gen date = yq(year,quarter)
	format date %tq
	tsset date
	sort date
	gen time = _n

	
* Gen Log versions of the variables
	gen Y = log(GNP)
	label var Y "ln(GNP)"
	gen MB = log(AMBSL)
	label var MB "ln(AMBSL) = log(StLouis adjusted monetary base, seasonally adjusted)

* Gen quarterly dummy variables
	gen q1 = 0
	gen q2 = 0
	gen q3 = 0
	gen q4 = 0
	replace q1 = 1 if quarter ==1
	replace q2 = 1 if quarter ==2
	replace q3 = 1 if quarter ==3
	replace q4 = 1 if quarter ==4
	
* VAR-GRANGER 
	*var Y MB, lags(1/8) exog(time quarter)
	*vargranger 
	
	quietly var Y MB, lags(1/8) exog(time q1-q4)
	vargranger 

	
* SIMS' UNIQUE APPROACH
	qui reg Y L(0/8).MB F(1/4).MB time quarter
		test F1.MB F2.MB F3.MB F4.MB
	qui reg MB  L(0/8).Y  F(1/4).Y time quarter
		test F1.Y F2.Y F3.Y F4.Y
	
* Exercise 1 (1970-2016)
	freduse GNP AMBSL, clear
	drop if year(daten)< 1970 | year(daten)> 2016
	gen year = year(daten)
	gen quarter = quarter(daten)
	collapse (mean) GNP AMBSL, by(year quarter)
	gen date = yq(year,quarter)
	format date %tq
	tsset date
	sort date
	gen time = _n
	gen Y = log(GNP)
	gen MB = log(AMBSL)
	qui var Y MB, lags(1/8) exog(time quarter)
	vargranger 





	



*************************
* 10.10.2 Indirect causality
*************************

clear all
set more off
local obs = 10000
set obs `obs'
gen time = _n-1
tsset time

matrix sigma = (1, 0, 0 \ 0, 1, 0 \ 0, 0, 1)
drawnorm e1 e2 e3, means(0,0,0) cov(sigma) double seed(1234)

gen X = 0
gen Z = 0
gen Y = 0

quietly{
	forvalues i = 2/`obs'{
		replace X = 0.00*L.Y + 0.00*L.X + 0.00*L.Z + e1 in `i'
		replace Z = 0.00*L.Y + 0.90*L.X + 0.00*L.Z + e2 in `i'
		replace Y = 0.00*L.Y + 0.00*L.X + 0.90*L.Z + e3 in `i'
	}
}
quietly varbasic X Z Y, lags(1)
vargranger


quietly varbasic X Z Y, lags(1)
	*graph export IRF_indirect_01.eps
	*graph export IRF_indirect_01.pdf, replace
		
	


********************************************************************
* 10.11 VAR Example: GNP and Unemployment
********************************************************************

clear all

freduse GNPC96 USAURHARMQDSMEI 
gen time=yq(year(daten), quarter(daten))  
tsset time, quarterly 

gen GNPgr = (ln(GNPC9) - ln(L.GNPC9))*100

label var GNPC9 "Real GNP"
label var USAUR "Unemployment Rate"
label var GNPgr "Real GNP growth rate (%)"

rename USAUR Unemp 

kpss GNPgr , notrend
kpss Unemp, notrend
kpss D.Unemp , notrend

* Graph the data
	twoway connected GNPC9 time, msymbol(none) xtitle("") name("VAR_example_GNPC9", replace) xlabel( ,alternate)
	twoway connected GNPgr time, msymbol(none) xtitle("") name("VAR_example_GNPgr", replace) xlabel( ,alternate)
	twoway connected Unemp time, msymbol(none) xtitle("") name("VAR_example_Unemp", replace) xlabel( ,alternate)
	twoway connected D.Une time, msymbol(none) xtitle("") name("VAR_example_DUnemp", replace) xlabel( ,alternate) ytitle("D.Unemployment") 
	graph combine VAR_example_GNPC9  VAR_example_Unemp  VAR_example_GNPgr  VAR_example_DUnemp
	graph export VAR_example_data.pdf, replace
	*graph export VAR_example_data.eps, repace

varsoc GNPgr D.Unemp, maxlag(8)

varbasic D.Unemp GNPgr , lag(1/2) irf
	*graph export VAR_example_IRF.pdf, replace
	*graph export VAR_example_IRF.eps, replace

varbasic D.Unemp GNPgr , lag(1/2) fevd
	*graph export VAR_example_FEVD.pdf, replace
	*graph export VAR_example_FEVD.eps, replace
	
varstable , graph modlabel
	*graph export VAR_example_varstable.pdf, replace
	*graph export VAR_example_varstable.eps, replace

varlmar
vargranger	





	
	

**********************************
* 10.12 Exercises
**********************************

* Exercise 2 
	* X_t = \beta_{xx1} L1.X + \beta_{xy1} L1.Y  + \beta_{xx2}L2.X  + \beta_{xy2}L2.Y  +  \epsilon_{x,t}
	* Y_t = \beta_{yx1} L1.X + \beta_{yy1} L1.Y  + \beta_{yx2}L2.X  + \beta_{yy2}L2.Y  +  \epsilon_{y,t}
*Write out the companion matrix:
	matrix beta = (0.30, -0.20, 0.15, 0.05 \ 0.30, -0.10, 0.05, 0.01 \ 1, 0, 0, 0 \ 0, 1, 0, 0 )
	matrix list beta
* Is this VAR stable?
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
* The lengths of the eigenvalues are 
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)
	display sqrt(real[1,3]^2 + immaginary[1,3]^2)
	display sqrt(real[1,4]^2 + immaginary[1,4]^2)
* IRF
	matrix beta2 = beta*beta
	matrix beta3 = beta*beta*beta
	matrix beta4 = beta*beta*beta*beta
	matrix list beta
	matrix list beta2
	matrix list beta3
	matrix list beta4
		
	
**********************************
* 10.06: Long-run levels: Including a constant
* and
* 10.12: Exercise 5 
**********************************

* VAR(1)
clear all
set more off
local obs = 10000
set obs `obs'
gen time = _n-1
tsset time

gen epsilonX = 0
gen epsilonY = 0

* What is the long-run level?
	gen X = 85.71429
	gen Y = 78.57143

	label var X "X{subscript:t}"
	label var Y "Y{subscript:t}"
	replace epsilonX = 0 if t == 1
	quietly{
		forvalues i = 2/`obs'{
		replace X = 100 + 0.20*L.X - 0.40*L.Y + epsilonX in `i'
		replace Y = 120 - 0.30*L.X - 0.20*L.Y + epsilonY in `i'
		}
	}
	*tsline Y X
	*sum 

* What is the IRF?	
drop epsilon* X Y

matrix sigma = (1, 0 \0, 1)
drawnorm epsilonX epsilonY, means(0,0) cov(sigma) double seed(1234)

gen Y = 0
gen X = 0

quietly{
	forvalues i = 2/`obs'{
		replace X = 100 + 0.20*L.X - 0.40*L.Y + epsilonX in `i'
		replace Y = 120 - 0.30*L.X - 0.20*L.Y + epsilonY in `i'
	}
}
qui var Y X , lags(1) nocons 
irf set temp, replace
irf create order1, step(10) 
irf table irf, irf(order1) noci step(4)

	
	
	
	
	
	