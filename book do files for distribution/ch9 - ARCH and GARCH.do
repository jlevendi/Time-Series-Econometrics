version 12.1
set scheme s1mono


*********************************
* Chapter 9: ARCH and GARCH
*********************************


****************************************
* Pictures of non-stationary variance
* A) Traditional heteroskedasticity
*    i) increasing variance
*   ii) decreasing variance
* B) Volatility clustering
*    i) artificial data
*   ii) real data
****************************************


* A) Traditional heteroskedasticity:
	drop _all
	clear all
	set more off
	set obs 250
	set seed 2345
	gen time = _n
	tsset time
	gen X = runiform()*100
	label var X " "
	
	gen e = rnormal()*30	// Normal errors
	gen u = e*6*X			// Errors (positively) proportional to X
	gen v = e*60/X			// Errors (inversely) proportional to X
	
	
	gen W = 10 + 2*X + e
	label var W " "
	scatter W X, title(No heteroskedasticity) ylabel(none) xlabel(none)
	graph save het_0, replace
	graph export het_0.pdf, replace


	gen Y = 10 + 2*X + u
	label var Y " "
	scatter Y X, title(Increasing variance) ylabel(none) xlabel(none)
	graph save het_1, replace
	graph export het_1.pdf, replace
	
	
	gen Z = 10 + 2*X + v
	label var Z " "
	scatter Z X, title(Decreasing variance) ylabel(none) xlabel(none)
	graph save het_2, replace
	graph export het_2.pdf, replace
	
	
	graph combine het_1.gph het_2.gph, rows(1) cols(2)  saving(het_1_2.gph, replace)
	graph export het_1_2.pdf, replace



* B) Volatility clustering

	* i) Artificial data
	
		* Procedure:
		* (1) generate sigma--a number that describes the variance of the error--as an AR process
		* (2) at each observation, draw an error term from a normal distribution with a stdev of sigma
	
		* Initialize
			drop _all
			local observations 300
			set obs `observations'
			set seed 2345
			gen time = _n
			tsset time
		
		* (1) Generate the stdev of X as an AR(1) process
			* gen e = runiform()
			* gen sigma = .
			* replace sigma = 12 in 1
			* replace sigma = 0.10 + .95*l.sigma + e in 2/L
		 	* graph twoway (line sigma time)
			* sum 
	
			* Generated as a simple regime-switching type model of high and low volatility
			gen sigma = 1
			replace sigma = 10 if (time < 20) | (time > 50 & time < 60) | (time > 90 & time < 130) | (time > 190 & time < 250)
	
		* (2) Generate the data, drawing errors from N(0,sigma)
			gen x = .
			forvalues t = 1/`observations'{
				quietly sum sigma in `t'
				quietly local sigma = r(mean)
				quietly replace x = rnormal(0,`sigma') in `t'
			}
			graph twoway (line x time, ylabel(none) title(Artificial Data))
			graph save   volatility_clustering    , replace
			graph export volatility_clustering.pdf, replace



*	* Real data (National Bank of Greece)
*		
*		drop _all	
*		fetchyahooquotes nbg, freq(d) chg(ln)
*		tsset date
*		drop adj
*		rename ln_nbg NBG_returns
*		label var NBG "Daily returns of NBG"
*		
*		gen year=year(date)
*		gen month = month(date)
*
*		graph twoway (line NBG_returns date if (year > 2005 & year < 2014), title(National Bank of Greece returns))
*		graph save   NBG    , replace
*		graph export NBG.pdf, replace
*
*
	graph combine volatility_clustering.gph NBG.gph, rows(1) cols(2)  saving(vol_clust_combined.gph, replace)
	graph export vol_clust_combined.pdf, replace


****************************************
* Simulating various GARCH data
****************************************

**************************
* 9.3.1 ARCH(1)
* A simple ARCH(1) model
**************************

* Y_t = beta0 + e_t
* e_t = u_t * (alpha0 + alpha1*(e_{t-1})^2)^(1/2)
* u_t = rnormal(0,1)

drop _all
set obs 1000
set seed 12345
gen time = _n
tsset time

* Enter the parameter values
local beta0 = 10
local alpha0 = 0.4
local alpha1 = 0.5

* Generate the data
gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1
replace e = u*(`alpha0' + `alpha1'*(L.e^2))^(1/2) in 2/L
gen Y = `beta0' + e

*save "arch-1.dta", replace

arch Y, arch(1)
predict varhat1, variance

graph twoway line Y time
	graph save ARCH-1.gph, replace
	graph export ARCH-1.pdf, replace

qnorm e
	graph save qqARCH-1.gph, replace
	graph export qqARCH-1.pdf, replace

histogram e, normal
	graph save histARCH-1.gph, replace
	graph export histARCH-1.pdf, replace
	
tabstat Y, statistics(kurtosis)

* Predicting the conditional variance:
	tsline varhat1
		graph save ARCH-1variance.gph, replace
		graph export ARCH-1variance.pdf, replace

* Predicting the unconditional variance
	local new = _N + 100
	set obs `new'
	replace time = _n
	predict varhat2, variance
	*predict varhat2, variance dynamic(.)
	*predict varhat3, variance dynamic(1000)
	*list in 995/1010
	display "Unconditional variance is: " [ARCH]_b[_cons] / (1 - [ARCH]_b[l1.arch])
	* .87119537
	tsline varhat2
		graph save ARCH-1uncvariance.gph, replace
		graph export ARCH-1uncvariance.pdf, replace





**************************
* 9.3.2 AR(1)-ARCH(1)
**************************


* Y_t = beta0 + beta1*Y_{t-1} + e_t
* e_t = u_t * (alpha0 + alpha1*(e_{t-1})^2)^(1/2)
* u_t = rnormal(0,1)

drop _all
set obs 1000
set seed 12345
gen time = _n
tsset time

local beta0 = 10
local beta1 = 0.10
local alpha0 = 0.40
local alpha1 = 0.50

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1
replace e = u*(`alpha0' + `alpha1'*(L.e^2))^(1/2) in 2/L
gen Y = .
replace Y = 11 in 1
replace Y = `beta0' + `beta1'*L.Y + e in 2/L

arch Y L.Y, arch(1)

graph twoway line Y time
	graph save AR-1-ARCH-1.gph, replace
	graph export AR-1-ARCH-1.pdf, replace

qnorm Y
	graph save qqAR-1-ARCH-1.gph, replace
	graph export qqAR-1-ARCH-1.pdf, replace
	
histogram Y, normal
	graph save histAR-1-ARCH-1.gph, replace
	graph export histAR-1-ARCH-1.pdf, replace


tabstat Y, statistics(kurtosis)



**************************
* 9.3.3 ARCH(2)
**************************

* Y_t = beta0 + e_t
* e_t = u_t * (alpha0 + alpha1*(e_{t-1})^2 + alpha2*(e_{t-2})^2 )^(1/2)
* u_t = rnormal(0,1)

drop _all
set obs 1000
set seed 12345
gen time = _n
tsset time

local beta0 = 10
local alpha0 = 0.20
local alpha1 = 0.30
local alpha2 = 0.40

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1/2
replace e = u*(`alpha0' + `alpha1'*(L.e^2) + `alpha2'*(L2.e^2) )^(1/2) in 3/L
gen Y = `beta0' + e

graph twoway line Y time
	graph save AR-0-ARCH-2.gph, replace
	graph export AR-0-ARCH-2.pdf, replace

qnorm Y
	graph save qqAR-0-ARCH-2.gph, replace
	graph export qqAR-0-ARCH-2.pdf, replace
	
histogram Y, normal
	graph save histAR-0-ARCH-2.gph, replace
	graph export histAR-0-ARCH-2.pdf, replace

tabstat Y, statistics(kurtosis)

arch Y, arch(1/2)



**************************
* 9.3.4: Testing for ARCH (example)
**************************

drop _all
*fetchyahooquotes GM, freq(d) start(01jan2000) end(11nov2013) 
*gen time = _n
*tsset time
*gen Y = ln(adjclose_) - ln(L.adjclose_)
* save "ARCH-GM.dta", replace
use "ARCH-GM.dta", clear

reg Y
estat archlm, lags(1/10)

*Table 9.1: AIC and BIC lag-selection
	forvalues lags=1/10{
	display "Lags = " `lags'
	quietly arch Y, arch(1/`lags')
	estat ic
	di ""
	di ""
	di ""
	}
	*






**************************
* 9.3.5 Example 1: Toyota Motor Company
**************************

* Load data
	*fetchyahooquotes TM, freq(d) start(01jan2000) end(31dec2010) 
	*gen time = _n
	*tsset time
	*gen TM = ln(adjclose_TM) - ln(L.adjclose_TM)
	*save "ARCH-TM.dta", replace
	use "ARCH-TM.dta", replace

* Graph data
	tsline TM
	*graph save ARCH-TM.gph, replace
	*graph export ARCH-TM.pdf, replace

* LM test quickly:
	quietly reg TM
	estat archlm, lags(1/10)

* Or by hand:
	quietly reg TM
	predict e, resid
	gen e2 = e^2
	drop e
	reg e2 L(1/10).e2

test L1.e2 L2.e2 L3.e2 L4.e2 L5.e2 L6.e2 L7.e2 L8.e2 L9.e2 L10.e2 

* Even more precisely, the "by hand" test should be:
	scalar teststat = e(r2)*e(N)
	display e(r2)*e(N)
	display chi2tail(10,teststat)

*Ljung-Box pre-estimation test
	wntestq e2

* There are ARCH effects of at least lag length = 1. How many lags should we include?
	forvalues lags=1/10{
	display "Lags = " `lags'
	qui arch TM, arch(1/`lags')
	estat ic
	di ""
	di ""
	}
	*

* Given that there are arch effects of length at least equal to 10, let's estimate the ARCH model:
	arch TM, arch(1/10)

* Do the coefficients indicate stationarity?
	display [ARCH]L1.arch +  [ARCH]L2.arch + [ARCH]L3.arch + [ARCH]L4.arch + ///
		[ARCH]L5.arch + [ARCH]L6.arch + [ARCH]L7.arch + [ARCH]L8.arch + ///
		[ARCH]L9.arch + [ARCH]L10.arch
	test [ARCH]L1.arch +  [ARCH]L2.arch + [ARCH]L3.arch + [ARCH]L4.arch + ///
		[ARCH]L5.arch + [ARCH]L6.arch + [ARCH]L7.arch + [ARCH]L8.arch + ///
		[ARCH]L9.arch + [ARCH]L10.arch  = 1
	* yes, stationary







**************************
* 9.3.6 Example 2: Ford Motor Company 
**************************

* Load the data
	*fetchyahooquotes F, freq(d) start(01jan1990) end(31dec1999) 
	*gen time = _n
	*tsset time
	*gen F = ln(adjclose_F) - ln(L.adjclose_F)
	*save "ARCH-F.dta", replace
	use "ARCH-F.dta", replace

* Graph the data
	tsline F
	*graph save ARCH-Ford.gph, replace
	*graph export ARCH-Ford.pdf, replace

* Engle's LM test 
	quietly reg F
	estat archlm, lags(1/5)

* Ljung-Box
	quietly reg F
	predict e, resid
	gen e2 = e^2
	wntestq e2, lags(5)

* What lag length ARCH model?
	forvalues lags=1/20{
	display "Lags = " `lags'
	qui arch F, arch(1/`lags')
	estat ic
	di ""
	di ""
	}
	*

* Given that there are ARCH, let's estimate the ARCH model:
	arch F, arch(1/15) nolog

* Do the coefficients indicate stationarity?
	display [ARCH]L1.arch +  [ARCH]L2.arch + [ARCH]L3.arch + ///
		[ARCH]L4.arch + [ARCH]L5.arch + [ARCH]L6.arch + ///
		[ARCH]L7.arch + [ARCH]L8.arch + [ARCH]L9.arch + ///
		[ARCH]L10.arch + [ARCH]L11.arch + [ARCH]L12.arch + ///
		[ARCH]L13.arch +  [ARCH]L14.arch + [ARCH]L15.arch

	test [ARCH]L1.arch +  [ARCH]L2.arch + [ARCH]L3.arch + ///
		[ARCH]L4.arch + [ARCH]L5.arch + [ARCH]L6.arch + ///
		[ARCH]L7.arch + [ARCH]L8.arch + [ARCH]L9.arch + ///
		[ARCH]L10.arch + [ARCH]L11.arch + [ARCH]L12.arch + ///
		[ARCH]L13.arch +  [ARCH]L14.arch + [ARCH]L15.arch = 1






****************************************
* 9.4 GARCH models
****************************************





**************************
* 9.4.1 GARCH(1,1)
**************************

* Y_t = beta0 + e_t
* e_t = sigma_t * u_t 
* sigma^2_t =  alpha0 + alpha1*e^2_{t-1} + gamma1*sigma^2_{t-1}
* u_t = rnormal(0,1)

drop _all
set obs 5000
set seed 12345
gen time = _n
tsset time

local beta0 = 10
local alpha0 = 0.2
local alpha1 = 0.4
local gamma1 = 0.6

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1
gen e2 = .
replace e2 = e^2 
gen sigma2 = .
replace sigma2 = 1 in 1

forvalues i=2/`=_N'{
	display "i is  `i'"
	replace sigma2 =  `alpha0' + `alpha1'*L.e2 +  `gamma1'*L.sigma2 in `i'
	replace e = sqrt(sigma2)*u in `i'
	replace e2 = e^2 in `i'
}


gen Y = `beta0' + e

arch Y , arch(1) garch(1)


graph twoway line Y time
	*graph save AR-0-GARCH-1-1.gph, replace
	*graph export AR-0-GARCH-1-1.pdf, replace


	
	
	


*********************
* 9.4.2 GARCH(2,1) example on simulated data
*********************

* Simulate the data
drop _all
set obs 100000
set seed 345
gen time = _n
tsset time

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1/3
gen e2 = .
replace e2 = e^2 
gen sigma2 = .
replace sigma2 = 1 in 1/3

quietly{
	forvalues i=4/`=_N'{
		replace sigma2 =  0.10 + 0.20*L.e2 + 0.30*L2.e2  + 0.20*L1.sigma2  in `i'
		replace e = sqrt(sigma2)*u in `i'
		replace e2 = e^2 in `i'
	}
}
gen Y = 0.10 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , arch(1/2) garch(1/1)

* It is always best to begin by graphing the data. 
tsline Y

* Is there volatility clustering? 
* Ljung-Box test
quietly reg Y
predict v, resid
gen v2 = v^2
wntestq v2
drop v v2

* Given the volatility clustering, which type of model best fits our data? AICs and BICs:
arch Y , arch(1) 
	estat ic
arch Y , arch(1) garch(1)
	estat ic
arch Y , arch(1) garch(1/2)
	estat ic
arch Y , arch(1/2) 
	estat ic
arch Y , arch(1/2) garch(1)
	estat ic
arch Y , arch(1/2) garch(1/2)
	estat ic

* Estimate the model:
arch Y, arch(1/2) garch(1/1)

* Post-estimation tests for specification: 
* Ljung-Box test
drop sigma2
predict w, resid
predict sigma2, variance
gen w2 = (w^2)/sigma2
wntestq w2
drop sigma2

* Are the coefficients significant?
test  [ARCH]L1.arch [ARCH]L2.arch [ARCH]L1.garch  
* yes

* Is the model stationary?
test  [ARCH]L1.arch + [ARCH]L2.arch + [ARCH]L1.garch = 1 
* yes

	
	


*********************
* 9.4.2 GARCH(p,q) example on IBM data
*********************

* Load the dataset
	*drop _all
	*fetchyahooquotes IBM, freq(d) start(01jan2001) end(31dec2010) //Download the data
	*gen time = _n
	*tsset time
	*gen IBM = ln(adjclose) - ln(L.adjclose) //Calculate the percentage daily returns. 
	*save ".\datasets\GARCH-IBM.dta", replace
	use ".\datasets\GARCH-IBM", replace

* Use a Ljung-Box test to determine whether these daily returns exhibit volatility clustering. 
	quietly reg IBM
	predict e, resid
	gen e2 = e^2
	wntestq e2

* Table 9.3
	* Given that our results indicate possible ARCH/GARCH effects, 
	* estimate an appropriate GARCH(p,q) model (using a q<= 2). 
		local maxp 2
		local maxq 2

	forvalues p = 1/`maxp'{
		display "ARCH(`p')"
		display "GARCH(0)"
		quietly arch IBM, arch(1/`p') 
		estat ic
		matrix rS = r(S)
		display "rS_AIC = " rS[1,5]
		display "rS_BIC = " rS[1,6]

		forvalues q = 1/`maxq'{
			display "ARCH(`p')"
			display "GARCH(`q')"
			quietly arch IBM, arch(1/`p') garch(1/`q')
			estat ic
			matrix rS = r(S)
			display "rS_AIC = " rS[1,5]
			display "rS_BIC = " rS[1,6]

		}
	}


* Estimate the model
	arch IBM, arch(1) garch(1) nolog

* Post-estimation Ljung-Box
	predict w, resid
	predict sigma2, variance
	gen w2 = (w^2)/sigma2
	wntestq w2

* Joint significance and stationarity
	test [ARCH]L1.arch [ARCH]L1.garch 
	test [ARCH]L1.arch + [ARCH]L1.garch =1







*****************************************
* 9.5.1 GARCH-t
* Replicating Bollerslev's (1987) GARCH(1,1)-t model
*****************************************
use ".\Bollerslev 1987\Bollerslev1987.dta", clear

foreach var of varlist USUK USDE{
	di "Variable is `var'. One of the exchange rate variables."
	gen Y = ln(`var') - ln(L.`var')

*Pre-tests
	qui reg Y
	qui predict e, resid
	qui gen e2 = e^2
	wntestq e, lags(10) 
	wntestq e2, lags(10) 
	tabstat Y, statistics(kurtosis) 
	qui drop e*

* Garch(1,1)-t estimate
	arch Y, arch(1) garch(1) distribution(t) nolog	
	qui predict e, resid
	qui predict v, variance

* Post-tests
	qui gen eSqrtV = e/sqrt(v)
	qui gen e2v = (e/sqrt(v))^2
	wntestq eSqrtV, lags(10)
	wntestq e2v, lags(10)
	tabstat eSqrtV, statistics(kurtosis)
	drop e* v Y

	di ""
	di "***************************************"
	di ""
}
*





****************************
* 9.5.2 GARCH-M or GARCH-IN-MEAN
* Example: GARCH-M in stock index returns
****************************

*cd "C:\Users\johnlevendis.BUSINESS\Google Drive\papers\book\ARCH\datasets\"
cd "C:\Users\johnlevendis\Google Drive\papers\book\ARCH\datasets\"

* Get the data
	drop _all
	set more off
	*fetchyahooquotes DJIA ^GSPC COMP ^IRX, freq(d) start(01jan1960) end(31dec2012) ff3
	*save "temp_indexes.dta", replace

* tsset
	*gen year = year(date)
	*sort date
	*gen time = _n
	*tsset time

* Rename variables
	*rename adjclose_DJIA DJIA
	*rename adjclose__GSPC SP
	*rename adjclose_COMP NASDAQ
	*rename adjclose__IRX TBILL
	
*drop ff3_M ff3_S ff3_H
*save "GARCH-M example.dta", replace

use "GARCH-M example.dta", clear

* Calculate excess returns
	gen retDJIA = log(DJIA) - log(L.DJIA) - ff3_RF			
	gen retSP = log(SP) - log(L.SP)  - ff3_RF
	gen retNASDAQ = log(NASDAQ) - log(L.NASDAQ) - ff3_RF 

* GARCH-M estimation and results
	qui arch retDJIA, arch(1) garch(1) archm 
		estimates store DJIA
	qui arch retSP, arch(1) garch(1) archm 
		estimates store SP
	qui arch retNASDAQ, arch(1) garch(1) archm 
		estimates store NASDAQ
	esttab DJIA SP NASDAQ,  star(* 0.10 ** 0.05 *** 0.01)






*********************
* GJR-GARCH(1,1) 
*********************

* Simulate the data
drop _all
set obs 100000
set seed 345
gen time = _n
tsset time

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1/3
gen e2 = .
replace e2 = e^2 
gen sigma2 = .
replace sigma2 = 1 in 1/3
gen D = 0

quietly{
	forvalues i=4/`=_N'{
		replace sigma2 =  0.05 + (0.20*L.e2 + 0.30*L.D*L.e2) + 0.30*L1.sigma2  in `i'
		replace e = sqrt(sigma2)*u in `i'
		replace e2 = e^2 in `i'
		replace D=1 if e >= 0 in `i'
	}
}
gen Y = 0.10 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , arch(1) garch(1) tarch(1)







*********************
* GJR-GARCH(2,1) 
*********************

* Simulate the data
drop _all
set obs 100000
set seed 345
gen time = _n
tsset time

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1/3
gen e2 = .
replace e2 = e^2 
gen sigma2 = .
replace sigma2 = 1 in 1/3
gen D = 0

quietly{
	forvalues i=4/`=_N'{
		replace sigma2 =  0.01 + (0.10*L.e2 + 0.20*L.D*L.e2) +   (0.15*L2.e2 + 0.25*L2.D*L2.e2) + 0.30*L1.sigma2  in `i'
		replace e = sqrt(sigma2)*u in `i'
		replace e2 = e^2 in `i'
		replace D=1 if e >= 0 in `i'
	}
}
gen Y = 0.05 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , arch(1/2) garch(1) tarch(1/2)







*********************
* GJR-GARCH(2,1) where only the first lag is asymmetric
*********************

* Simulate the data
drop _all
set obs 500000
set seed 345
gen time = _n
tsset time

gen u = rnormal(0,1)
gen e = .
replace e = 0 in 1/3
gen e2 = .
replace e2 = e^2 
gen sigma2 = .
replace sigma2 = 1 in 1/3
gen D = 0

quietly{
	forvalues i=4/`=_N'{
		replace sigma2 =  0.01 + (0.10*L.e2 + 0.20*L.D*L.e2) +  (0.15*L2.e2) + 0.30*L1.sigma2  in `i'
		replace e = sqrt(sigma2)*u in `i'
		replace e2 = e^2 in `i'
		replace D=1 if e >= 0 in `i'
	}
}
gen Y = 0.05 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , arch(1/2) garch(1) tarch(1)



******************************************
* 9.5.3 Asymmetric responses in GARCH
******************************************

*********************
* 9.5.3 E-GARCH(1,1) simulation
*********************

* Simulate the data
drop _all
set obs 500000
set seed 345
gen time = _n
tsset time

gen z = rnormal(0,1)
gen e = .
replace e = 0 in 1/3
gen sigma2 = .
replace sigma2 = 1 in 1/3
gen lnsigma2 = .
replace lnsigma2 = log(sigma2) in 1/3

quietly{
	forvalues i=4/`=_N'{
		replace lnsigma2 =  0.05 + 0.10*L.z +  0.20*(abs(L.z) - sqrt(2/_pi)) + 0.30*L1.lnsigma2  in `i'
		*replace lnsigma2  =  0.05 - 0.20*sqrt(2/_pi) + 0.10*L.z +  0.20*abs(L.z)  + 0.30*L1.lnsigma2  in `i'
		replace sigma2 = exp(lnsigma2) in `i'
		replace sigma = sqrt(sigma2) in `i'
		replace e = sigma*z in `i'
	}
}
gen Y = 0.10 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , earch(1) egarch(1)





*********************
* 9.5.3 E-GARCH(2,3) simulation
*********************

* Simulate the data
drop _all
set obs 100000
set seed 345
gen time = _n
tsset time

gen z = rnormal(0,1)
gen e = .
replace e = 0 in 1/4
gen sigma2 = .
replace sigma2 = 1 in 1/4
gen lnsigma2 = .
replace lnsigma2 = log(sigma2) in 1/4

local alpha0  0.05
local alpha11 0.20
local alpha12 0.15
local alpha21 0.10
local alpha22 0.05
local gamma1  0.15
local gamma2  0.10
local gamma3  0.05

quietly{
	forvalues i=4/`=_N'{
		replace lnsigma2 =  `alpha0' + (`alpha11'*L1.z +  `alpha12'*(abs(L1.z) - sqrt(2/_pi))) ///
							+ (`alpha21'*L2.z +  `alpha22'*(abs(L2.z) - sqrt(2/_pi))) ///
							+ `gamma1'*L1.lnsigma2  + `gamma2'*L2.lnsigma2  + `gamma3'*L3.lnsigma2  in `i'
		replace sigma2 = exp(lnsigma2) in `i'
		replace sigma = sqrt(sigma2) in `i'
		replace e = sigma*z in `i'
	}
}
gen Y = 0.10 + e

* To verify that our simulation was correct, we estimate the model:
arch Y , earch(1/2) egarch(1/3)







*********************
* 9.5.3 Example: Comparing Asymmetric GARCH models
*********************

* Get the data
drop _all
fetchyahooquotes DJIA, freq(d) start(01jan2005) end(31dec2012)
gen DJIA = log(adj) - log(L.adj)
*tsline DJIA

* (1) Estimate a normal GARCH(1,1) model
quietly arch DJIA, arch(1) garch(1) nolog
estimates store GARCH
predict GARCH_hat, variance

* (2) Estimate a GJR-GARCH(1,1) model
quietly arch DJIA, arch(1) garch(1) tarch(1) nolog
estimates store GJR_GARCH
predict GJR_GARCH_hat, variance

* (3) Estiate an E-GARCH(1,1) model
quietly arch DJIA, earch(1) egarch(1) nolog
estimates store E_GARCH
predict E_GARCH_hat, variance

* The output
* estimates table GARCH GJR_GARCH E_GARCH
esttab GARCH GJR_GARCH E_GARCH, not
estimates stats GARCH GJR_GARCH E_GARCH
correlate GARCH_hat E_GARCH_hat GJR_GARCH_hat






*********************
* 9.5.4 I-GARCH or Integrated GARCH
*********************

	*fetchyahooquotes DJIA, freq(d) start(01jan1970) end(31dec2015)
	*gen DJIA_returns = log(adj)-log(L.adj)
	*drop adj
	*tsset date
	*save "IGARCH-DJIA.dta", replace

	use "IGARCH-DJIA.dta", clear
	constraint 1 [ARCH]_b[L1.arch] + [ARCH]_b[L1.garch] = 1
	arch DJIA_returns, arch(1) garch(1) nolog constraints(1) 
	display [ARCH]_b[L1.arch] + [ARCH]_b[L1.garch] 


	
	
	
	
	
	
******************************************
* 9.6 EXERCISES 
******************************************

********************
*Exercise 4) LM tests of arch-1 data
********************
	use "arch-1.dta", clear
	reg Y
	estat archlm, lags(1/10)


********************
*Exercise 6) AICs and BICs on MMM
********************
	local var MMM
	local maxlag 20

	fetchyahooquotes `var', freq(d) start(01jan2000) end(31dec2010) 
	gen time = _n
	tsset time
	gen `var' = ln(adjclose_`var') - ln(L.adjclose_`var')



*********************
* Exercise 7) ARCH(10) vs GARCH(1,1)
*********************	
	
* We use the old Ford data
	use "ARCH-F.dta", clear"

* The two models
	arch F, arch(1/15) nolog
	estat ic
	predict varhat_arch, variance
	arch F, arch(1) garch(1) nolog
	estat ic
	predict varhat_garch, variance

* Show that the predicted variances are very similar
	correlate varhat_arch varhat_garch
	
	
********************
* Exercise 12) Bollerslev (1987)
********************
use "BollerslevSP.dta", clear

* Generate the necessary auxiliary data
	gen Y = .
	gen temp = log(adj/L.adj)
	sum temp
	replace Y = r(mean) in 1/5
	drop temp
	sum time
	local last = r(max)
	replace Y = log(adj/L.adj) - 0.2680*L.Y - 0.0718*L2.Y + 0.0192*L3.Y - 0.0052*L4.Y in 6/`last'

*Pre-tests
	qui reg Y
	qui predict e, resid
	qui gen e2 = e^2
	wntestq e, lags(10) /* Bollerslev gets Q = 15*/
	wntestq e2, lags(10) /* Bollerslev gets Q2 = 26.36*/
	tabstat e, statistics(kurtosis) /* Bollerslev gets K = 4.06 */
	qui drop e e2

* Garch(1,1)-t etimate
	arch Y, arch(1) garch(1) distribution(t) nolog	
	qui predict e, residual
	qui predict v, variance

* Post-tests
	gen eSqrtV = e/sqrt(v)
	gen e2v = (e/sqrt(v))^2
	wntestq eSqrtV, lags(10)	/* Bollerslev gets 14.08*/
	wntestq e2v, lags(10)	/* Bollerslev gets 8.08*/

* Conditional kurtosis? (Bollerslev gets 3.89)
	tabstat eSqrtV, statistics(kurtosis)	 






	
	
	
	
	
