***************************************
* Ch4: Stationarity and Invertibility 
***************************************

set scheme s1mono



****************************************
* Generate some data, some stationary, some not
****************************************
	drop _all
	set obs 50
	set seed 1234
	gen time = _n
	tsset time
	gen double error = rnormal()
	*W: ar(-0.80 0.30) 	
		gen W = error 
		replace W = -0.80*L.W + 0.30*L2.W + error in 3/L
	*Z: ar(0.80 0.30) 	 
		gen Z = error 
		replace Z = 0.80*L.Z + 0.30*L2.Z + error in 3/L
	*A: ar(-0.15 1.10)
		gen A = error 
		replace A = -0.15*L.A + 1.10*L2.A + error in 3/L
	*X: ar(0 1.10)
		sim_arma X, nobs(100) ar(1.10)          time(time) 
		arima X, ar(1/2) nolog nocons
	*Y: ar(0.70 0.10)
		sim_arma Y, nobs(100) ar(0.70 0.10)     time(time) 
	*C: ar(1.30 -0.50) 	
		sim_arma C, nobs(100) ar(1.30 -0.50) time(time)
	*B: ar(0.50 0.10)
		sim_arma B, nobs(100) ar(0.50 0.10)     time(time) 
	*D: linear + error
		gen D = time + error
	drop error
	order time X Y Z W A B C D
	save ".\data\stationarity_ex.dta", replace

	
			
	
	
*******************************************
* 4.3.2 Restrictions on AR(2) coeffients
*******************************************

* Figure 4.2
	drop _all
	set seed 1234
	sim_arma X, nobs(100) ar(1.10)          spin(2000) time(time) 
	sim_arma Y, nobs(100) ar(0.70 0.10)     spin(2000) time(time) 
	sim_arma Z, nobs(100) ar(0.80 0.30) 	spin(2000) time(time) 
	sim_arma W, nobs(100) ar(-0.80 0.30) 	spin(2000) time(time) 

	twoway line X time if time< 50, title(X = 1.10*L.X + e) ylabel(none)
	twoway line Y time, title(Y = 0.70*L.Y + 0.10*L2.Y + e) ylabel(none)
	twoway line Z time if time<50, title(Z = 0.80*L.Z + 0.30*L2.Z + e) ylabel(none)
	twoway line W time, title(W =-0.80*L.W + 0.30*L2.W + e) ylabel(none)

		
* Example, and Figure 4.3a
	use ".\data\stationarity_ex.dta", clear
	arima X, ar(1/2) nolog nocons
	estat aroots
	
* Example, and Figure 4.3b
	arima Y, ar(1/2) nolog nocons
	estat aroots
		
* Exercises
	qui arima Z, ar(1/2) diffuse nolog nocons
		estat aroots, nograph
	qui arima W, ar(1/2) diffuse nolog nocons
		estat aroots, nograph
	qui arima A, ar(1/2) diffuse nolog nocons
		estat aroots, nograph
	qui arima C, ar(1/2) diffuse nolog nocons
		estat aroots, nograph
	qui arima B, ar(1/2) diffuse nolog nocons
		estat aroots, nograph
	qui arima D, ar(1/2) nocons nolog
		estat aroots, nograph


		
************************************
* 4.3.4: Exercises: AR(2) stability
************************************

* Exercises
drop _all
set obs 50
gen t = _n
tsset t

foreach i in a b c d e f g{
	gen Y`i' = 0
	replace Y`i' = 1 in 1/2
}
/*a*/ replace Ya =-0.10*L.Ya + 0.15*L2.Ya in 3/L
/*b*/ replace Yb = 0.20*L.Yb + 0.80*L2.Yb in 3/L
/*c*/ replace Yc = 0.15*L.Yc + 0.80*L2.Yc in 3/L
/*d*/ replace Yd = 0.15*L.Yd - 0.80*L2.Yd in 3/L
/*e*/ replace Ye =-0.40*L.Ye + 0.60*L2.Ye in 3/L
/*f*/ replace Yf = 0.30*L.Yf + 0.10*L2.Yf in 3/L
/*g*/ replace Yg = 0.80*L.Yg + 0.10*L2.Yg in 3/L

/*Solution: a*/
	matrix beta = [-0.10, 0, 0.15 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: b*/
	matrix beta = [0.20, 0, 0.80 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: c*/
	matrix beta = [0.15, 0, 0.80 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: d*/
	matrix beta = [0.15, 0, -0.80 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: e*/
	matrix beta = [-0.40, 0, 0.60 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: f*/
	matrix beta = [0.30, 0, 0.10 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)

/*Solution: g*/
	matrix beta = [0.80, 0, 0.10 , 0 \ 0,0,0,0 \ 1,0,0,0\0,1,0,0]
	matrix list beta
	matrix eigenvalues real immaginary = beta
	matrix list real
	matrix list immaginary
	display sqrt(real[1,1]^2 + immaginary[1,1]^2)
	display sqrt(real[1,2]^2 + immaginary[1,2]^2)


		
*********************************
* Stralkowski triangle
*********************************

	twoway (function beta2 = 1+ x, range(-3 2) lcolor(black) ) ///
		(function beta2 = 1-x, range(-2 3) lcolor(black) ) ///
		(function beta2 = -1, range(-4 4) lcolor(black) ) ///
		(function beta2 = -(x^2)/4 , range(-3 3) lcolor(black) ) ///
		(function beta2 = 1+x, range(-2 0) lcolor(black) lwidth(thick) ) ///
		(function beta2 = 1-x, range(0 2) lcolor(black) lwidth(thick) ) ///
		(function beta2 = -1, range(-2 2) lcolor(black) lwidth(thick) ) ///
		, yline(0) xline(0) legend(off) plotregion(style(none)) ///
		xtitle("{&beta}{sub:1}") ytitle("{&beta}{sub:2}") 
	
	



	
	
	
	