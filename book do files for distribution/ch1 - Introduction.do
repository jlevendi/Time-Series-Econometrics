************************
* Introduction
************************

set scheme s1mono

* Taking lags and first-differences
	use ".\DandL.dta", clear
	list
	gen LagX = L.X
	gen SecLagX = L2.X
	gen Y = D.X		//Y is the first difference of X
	gen Z = D.Y		//Z is the first difference of Y, and 2nd diff of X.
	list

* Difference between cross-sections and time-series.
	* Cross-sections
		drop _all
		set seed 1234
		set obs 20
		gen errors = rnormal()
		gen X = rnormal()
		gen Y = 10 + 2*X + errors
		twoway scatter Y X

	* Time-series
		drop _all
		set seed 567
		set obs 20
		sim_arma Y, ar(1) time(time)
		replace Y = 2+Y
		twoway connected Y time

