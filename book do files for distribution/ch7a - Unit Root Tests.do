****************************************************
* 7.3.1 A Random Walk vs a zero-mean AR(1) process
****************************************************

use ".\data\RWdata.dta", clear
dfuller X, noconstant regress
tsline X /* Figure 7.1 */

use ".\data\ARexamples.dta", clear
keep if _n <= 100
dfuller X, noconstant regress
tsline X /* Figure 7.2 */


****************************************************
* 7.3.2 A Random Walk vs an AR(1) model with a constant
****************************************************

use ".\data\RWdata.dta", clear
dfuller X, regress

use ".\data\AR1nonzero_b.dta", clear
dfuller X, regress


****************************************************
* 7.3.3 A RandomWalk with Drift vs a Deterministic Trend
****************************************************

use ".\data\RWdata.dta", clear
dfuller Y, regress trend
tsline Y /* Figure 7.3 */

use ".\data\DTdata.dta", clear
dfuller X, regress trend
tsline X /* Figure 7.4 */


****************************************************
* 7.3.6 Choosing the lag length in DF-type tests
****************************************************
drop _all
freduse UNRATE
gen time = _n
tsset time
drop if time > 838
dfgls UNRATE


****************************************************
* 7.4 Phillips-Perron tests
****************************************************
use ".\data\RWdata.dta", clear
pperron Y, regress trend
pperron D.Y, trend


****************************************************
* 7.5 KPSS tests, Exercises
****************************************************
*net install sts15_2.pkg

use ".\data\NelsonPlosserData.dta", clear
keep l* year bnd
order lrgnp lgnp lpcrgnp lip lemp lun lprgnp lcpi lwg lrwg lm lvel bnd lsp500	
set more off

* Exercise 1
* 1a)
	foreach v of var lrgnp-lsp500 {
	display ""
	display ""
	display ""
	display ""
	display ""
	display "-----------------------------------------"
	display "`v'"
	kpss  `v' , maxlag(8)  notrend
	}

* 1b)
	foreach v of var lrgnp-lsp500 {
	display ""
	display ""
	display ""
	display ""
	display ""
	display "-----------------------------------------"
	display "`v'"
	kpss  `v' , maxlag(8) 
	}


* Exercise 2
* 1a)
	foreach v of var lrgnp-lsp500 {
	display ""
	display ""
	display ""
	display ""
	display ""
	display "-----------------------------------------"
	display "`v'"
	kpss  `v' , maxlag(8)  notrend auto qs
	}

* 1b)
	foreach v of var lrgnp-lsp500 {
	display ""
	display ""
	display ""
	display ""
	display ""
	display "-----------------------------------------"
	display "`v'"
	kpss  `v' , maxlag(8) auto qs
	}










