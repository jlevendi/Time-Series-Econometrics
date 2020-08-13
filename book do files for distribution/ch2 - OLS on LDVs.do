*****************************************
* OLS on a LDV is biased
*****************************************

* PLEASE NOTE: THIS FILE TAKES A LONG TIME TO RUN

drop _all

* generating fake time series for illustration purposes

	set seed 123
	local replications = 100000
	
	forvalues obs = 20(20)100{
		capture postclose buffer
		postfile buffer beta1_`obs' using ".\data\ols_on_ldv_a_`obs'.dta", replace
		
		forvalues i = 1/`replications'{
			qui drop _all
			qui set obs `obs'
			qui gen time = _n
			qui tsset time
			qui gen double error = rnormal()
			qui label var e "errors"

			* AR(1)
			qui gen double Y = .
			qui replace Y = error in 1
			qui replace Y = 0.50*L1.Y + error in 2/L
			
			qui reg Y L1.Y
			
			post buffer (_b[L.Y])
		}
		di "obs = `obs'"
		postclose buffer
	}
	*


* MERGING THE DATA
	drop _all
	use ".\data\ols_on_ldv_a_20.dta", clear
	merge 1:1 _n using ".\data\ols_on_ldv_a_40.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_a_60.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_a_80.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_a_100.dta"
	drop _merge	


* Kdensity
forvalues obs = 20(20)100{
	kdensity beta1_`obs', n(500) generate(x_`obs' f_`obs') kernel(gaussian) nograph
	label variable f_`obs' "N=`obs'"
}
save ".\data\ols_on_ldv_a.dta", replace
*

* Figure 2.3
use ".\data\ols_on_ldv_a.dta", clear
graph twoway (line f_20 x_20) (line f_40 x_40) (line f_60 x_60) (line f_80 x_80) (line f_100 x_100) , ///
	note("100,000 replications. Y{sub:t} = 0.50Y{sub:t-1}+ e{sub:t}") ///
	legend(cols(3))
	graph save ".\graphs\ols_on_ldv_is_biased.gph" , replace
	graph export ".\graphs\ols_on_ldv_is_biased.pdf", replace

* How bad is the bias?
	sum beta1*	
	
* Erase not-needed datafiles
	erase ".\data\ols_on_ldv_a.dta"
	erase ".\data\ols_on_ldv_a_20.dta"
	erase ".\data\ols_on_ldv_a_40.dta"
	erase ".\data\ols_on_ldv_a_60.dta"
	erase ".\data\ols_on_ldv_a_80.dta"
	erase ".\data\ols_on_ldv_a_100.dta"

	
	
	
	
	
	
	
	
	
	
*****************************************
* OLS on a LDV w/ AR errors is inconsistent
*****************************************

drop _all

* generating fake time series for illustration purposes

	set seed 123
	local replications = 100000
	
	forvalues obs = 20(20)100{
		capture postclose buffer
		postfile buffer beta1_`obs' using ".\data\ols_on_ldv_b_`obs'", replace
		
		forvalues i = 1/`replications'{
			qui drop _all
			qui set obs `obs'
			qui gen time = _n
			qui tsset time
			qui gen double u = rnormal()
			qui label var u "errors"

			* AR(1) errors
			qui gen e = u in 1
			qui replace e = 0.20*L1.e + u in 2/L
			
			* AR(1)
			qui gen double Y = .
			qui replace Y = e in 1
			qui replace Y = 0.50*L1.Y + e in 2/L
			
			qui reg Y L1.Y
			
			post buffer (_b[L.Y])
		}
		di "obs = `obs'"
		postclose buffer
	}
	*


* MERGING THE DATA
	drop _all
	use ".\data\ols_on_ldv_b_20.dta", clear
	merge 1:1 _n using ".\data\ols_on_ldv_b_40.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_b_60.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_b_80.dta"
	drop _merge
	merge 1:1 _n using ".\data\ols_on_ldv_b_100.dta"
	drop _merge	
	

* Kdensity
forvalues obs = 20(20)100{
	kdensity beta1_`obs', n(500) generate(x_`obs' f_`obs') kernel(gaussian) nograph
	label variable f_`obs' "N=`obs'"
}
save ".\data\ols_on_ldv_b.dta", replace

* Figure 2.4
use ".\data\ols_on_ldv_b.dta", clear
graph twoway (line f_20 x_20) (line f_40 x_40) (line f_60 x_60) (line f_80 x_80) (line f_100 x_100) , ///
	note("100,000 replications. Y{sub:t} = 0.50Y{sub:t-1}+ e{sub:t} with e{sub:t} = 0.20e{sub:t-1} + u{sub:t}") ///
	legend(cols(3))
	graph save ".\graphs\ols_on_ldv_is_inconsistent.gph" , replace
	graph export ".\graphs\ols_on_ldv_is_inconsistent.pdf", replace

* How bad is the bias?
	sum beta1*
	
* Erase not-needed datafiles
	erase ".\data\ols_on_ldv_b.dta"
	erase ".\data\ols_on_ldv_b_20.dta"
	erase ".\data\ols_on_ldv_b_40.dta"
	erase ".\data\ols_on_ldv_b_60.dta"
	erase ".\data\ols_on_ldv_b_80.dta"
	erase ".\data\ols_on_ldv_b_100.dta"
	


	
	


*****************************************
* OLS on LDVs snippet
*****************************************

drop _all
set seed 123
postfile buffer beta1 using filename, replace

forvalues i = 1/1000 	/*The number of replications*/ {
	qui drop _all
	qui set obs 20 		/*The number of observations*/
	qui gen time = _n
	qui tsset time
	qui gen double error = rnormal()
		
	* Generate Y
	qui gen double Y = .
	qui replace Y = error in 1
	qui replace Y = 0.50*L1.Y + error in 2/L
	
	qui reg Y L1.Y		/*Estimates beta1*/
	
	post buffer (_b[L.Y]) /*Saves the estimate*/
	}

postclose buffer		/*Closes the file with all the estimates*/


* Kdensity
	use filename.dta , clear
*	kdensity beta1, kernel(gaussian) 
	kdensity beta1 


	

	
