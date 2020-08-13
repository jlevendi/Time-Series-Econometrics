version 12.1

*********************************
* Zivot and Andrews structural break test
*********************************

use "NelsonPlosserData.dta" , clear
	
keep l* year bnd
order year lrgnp lgnp lpcrgnp lip lemp lprgnp lcpi lwg lm lvel bnd lrwg lsp500
drop lun


foreach var in lrgnp lgnp lpcrgnp lip lemp lprgnp lcpi lwg lm lvel bnd {
quietly{
	local optimal_t = 100
	local optimal_year = 0

	forvalues bp = 1860(1)1970{
	capture drop DL
	gen DL = 0
	replace DL = 1 if year > `bp'
		
		local lags = 8
		reg `var' year DL    Ld(0/`lags').`var' 
		local tstat = _b[Ld`lags'.`var']/_se[Ld`lags'.`var']
		local alpha_tstat =   (_b[L1.`var']-1)/_se[L1.`var']
		if `alpha_tstat' < `optimal_t' {
			local optimal_t = `alpha_tstat'
			local optimal_year = `bp'
			}

		while abs(`tstat')<1.60 & `lags'>=0 {
			local lags = `lags'-1
			reg `var' year DL         Ld(0/`lags').`var' 
			local tstat = _b[Ld`lags'.`var']/_se[Ld`lags'.`var']
			if `alpha_tstat' < `optimal_t' {
				local optimal_t = `alpha_tstat'
				local optimal_year = `bp'
			}
		}
	}

		capture drop DL
		gen DL = 0
		replace DL = 1 if year > `optimal_year'
		* Given the optimal year, what is the optimal K?
			local lags = 8
			reg `var' year DL    Ld(0/`lags').`var' 
			local tstat = _b[Ld`lags'.`var']/_se[Ld`lags'.`var']
			while abs(`tstat')<1.60 & `lags'>=0 {
				local lags = `lags'-1
				reg `var' year DL         Ld(0/`lags').`var' 
				local tstat = _b[Ld`lags'.`var']/_se[Ld`lags'.`var']
			}
	}


quietly reg `var' year DL    Ld(0/`lags').`var' 

display ""

display "`var'"
display "   Breakpoint: " `optimal_year'
display "   K (lags):   " `lags'
display "   Beta(DL): " %4.3f  _b[DL]       ///
		"   t(DL): (" %4.2f _b[DL]/_se[DL] ")" 
display "   Beta(t):   " %4.3f _b[year]       ///
		"   t(t):  (" %4.2f _b[year]/_se[year] ")"
display "   Alpha:  " %4.3f _b[l1.`var']  ///
		"   t(alpha): (" %4.2f (_b[L1.`var']-1)/_se[L1.`var'] ")"
display "   S(e):   " %4.2f e(rmse)
display "------------------------------------------"
display ""
}
*











*******************************************
* Table 8.3: zandrews minimum t-statistucs
*******************************************

zandrews lrgnp, break(intercept) lagmethod(input) maxlags(8)
zandrews lgnp , break(intercept) lagmethod(input) maxlags(8)
zandrews lpcrgnp  , break(intercept) lagmethod(input) maxlags(7)
zandrews lip  , break(intercept) lagmethod(input) maxlags(8)
zandrews lemp , break(intercept) lagmethod(input) maxlags(7)
zandrews lprgnp  , break(intercept) lagmethod(input) maxlags(5)
zandrews lcpi , break(intercept) lagmethod(input) maxlags(2)
zandrews lwg  , break(intercept) lagmethod(input) maxlags(7)
zandrews lm   , break(intercept) lagmethod(input) maxlags(6)
zandrews lvel , break(intercept) lagmethod(input) maxlags(2)
zandrews bnd  , break(intercept) lagmethod(input) maxlags(2)

zandrews lsp500, break(both)  lagmethod(input) maxlags(1)
zandrews lrwg, break(both)  lagmethod(input) maxlags(8)

* For all of these, the date of structural change reported in zandrews is always one 
* year later than the date reported in the paper by Zivot and Andrews (1992).
* This is just a labelling issue.




