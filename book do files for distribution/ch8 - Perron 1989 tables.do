

set scheme s1mono


************************************************************************
* PERRON (1989) replication
* "The Great Crash, the Oil Price Shock, and the Unit Root Hypothesis"
************************************************************************


********************************
* Getting the data in order
********************************

use NelsonPlosserData.dta, clear

set more off
keep l* year bnd


* Dummy variable indicating all periods after the Great Depression
	gen DL = 0
	replace DL = 1 if year > 1929

* Dummy variable indicating the one period after the shock
	gen DP = 0
	replace DP = 1 if year == 1930

gen NewSlope = year*DL



************************************************
* Table 8.1
* sample Autocorrelations of the Detrended Data
************************************************


foreach var of varlist lrgnp lgnp lpcrgnp lip lemp lprgnp lcp lwg lm lvel bnd {
	quietly reg `var'  year DL
	quietly predict dt_`var', resid
	label var dt_`var' "Detrended `var'"
	display "`var'"
	corrgram dt_`var', lags(6) noplot
	display ""
	display "*******************************************"
	display ""
}
*

foreach var of varlist lrwg lsp500 {
	quietly reg `var'  year DL NewSlope
	quietly predict dt_`var', resid
	label var dt_`var' "Detrended `var'"
	display "`var'"
	corrgram dt_`var', lags(6) noplot
	display ""
	display "*******************************************"
	display ""
}
*


**************************************************
* Table 8.2
* Tests For A Unit Root, w/ known structural break
**************************************************


reg lrgnp DL year DP L.lrgnp    L(1/8)D.lrgnp

*******************************************

local var = "lrgnp"
local lags = "8"
	reg `var' DL year DP L.`var' L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lgnp"
local lags = "8"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lpcrgnp"
local lags = "7"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lip"
local lags = "8"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lemp"
local lags = "7"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************


local var = "lprgnp"
local lags = "5"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lcpi"
local lags = "2"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lwg"
local lags = "7"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************




local var = "lm"
local lags = "6"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "lvel"
local lags = "0"
	*quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	quietly reg `var' DL year DP L.`var'    
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



local var = "bnd"
local lags = "2"
	quietly reg `var' DL year DP L.`var'    L(1/`lags')D.`var'
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************



* Real Wages and the S&P-500 have a different model specification.
* Same spec as above, but we add DT

local var = "lsp500"
local lags = "1"
	quietly reg `var' DL year DP L.`var'  NewSlope  L(1/`lags')D.`var' 
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "beta(NewSlope) = " %4.3f _b[NewSlope]
	display "   t(NewSlope) = " %4.2f _b[NewSlope]/_se[NewSlope]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************


local var = "lrwg"
local lags = "8"
	quietly reg `var' DL year DP L.`var'  NewSlope  L(1/`lags')D.`var' 
	display "Variable = `var'"
	display "      beta(DL) = " %4.3f _b[DL]
	display "         t(DL) = " %4.2f _b[DL]/_se[DL]
	display "    beta(time) = " %4.3f _b[year]    
	display "       t(time) = " %4.2f _b[year]/_se[year]
	display "      beta(DP) = " %4.3f _b[DP]               
	display "         t(DP) = " %4.2f _b[DP]/_se[DP]
	display "beta(NewSlope) = " %4.3f _b[NewSlope]
	display "   t(NewSlope) = " %4.2f _b[NewSlope]/_se[NewSlope]
	display "         alpha = " %4.3f _b[L.`var']       
	display "      t(alpha) = " %4.2f (_b[L.`var']-1)/_se[L.`var']
	display "      S(e hat) = " %4.3f e(rmse)
*******************************************












