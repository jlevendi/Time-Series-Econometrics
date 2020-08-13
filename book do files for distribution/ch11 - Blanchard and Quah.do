*****************************
* Blanchard and Quah (via Schenck)
*****************************
cd "C:\Users\johnlevendis\Google Drive\papers\book\VAR\graphs\"
set scheme s1mono


*freduse GNPC96 UNRATE, clear
*save BQdata.dta, replace
use BQdata.dta, clear

gen year = year(daten)
gen quarter = quarter(daten)

collapse (mean) GNP UNRATE (first) daten, by(year quarter)
sort year quarter
gen time = _n
tsset time

gen GNPgr = 400*(log(GNPC96) - log(L.GNPC96))

*This follows Blanchard and Quah
	gen quarterly = yq(year,quarter)
	format quarterly %tq
	keep if quarterly >= yq(1952,2) & quarterly <= yq(1987,4)

*This follows Schenck
	*keep if year>=1952 & year <= 1987	

quietly regress UNRATE time
predict unemp, residual

gen Pre1974 = 0
gen Post1974 = 0
replace Pre1974 = 1 if year < 1974
replace Post1974 = 1 if year >=1974

quietly reg GNPgr Pre1974 Post1974, noconstant
predict Ygr, residual

matrix C = (., 0 \ . , .)
svar Ygr unemp , lags(1/8) lreq(C)

irf create Blanchard_Quah, set(lrirf, replace) step(40) replace
irf graph sirf, yline(0) byopts(yrescale)
graph export BQ_IRF_1.pdf, replace



