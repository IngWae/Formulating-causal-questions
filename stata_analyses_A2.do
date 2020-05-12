*cd "C:\Users\slecessie\Dropbox\Stratos-causality\May June 2018\Analyses with stata"
cd "C:\Users\Bianca DeStavola\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS"
*cd "C:\Users\sejjbld\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS"

*analyses of version 2 of the data (created in April 2018)
cap log close
set more off
set scheme s2mono
clear all
discard
set seed 13745

*log using "Analyses with stata\Appendix2_A2_stata.log", replace


*****
capture : file close myfile
file open myfile using "Analyses with stata\results_for_Appendix 2_0305202.txt", write replace
	file write myfile "Today is  $S_DATE " 
	
	file write myfile "Results for Appendix 2  ATE: " _n

*read new version of the data
import delimited "PROBITsim_april2018.csv"
des 
su
rename wgt3observed wgt3
rename a2observed a2
rename a3observed a3

/******************************
FOR THE TEXT
Estimate the causal effect of intervention (ATE of $A_1$) on infant's weight at 3 months: same as ITT
*/
reg wgt3 i.a1, base
ta a2 a1,col
ta a3 a1,col
ta a3 a2 if a1==1,col


/******************************
TABLE 5 OF THE PAPER:   Estimate the ATE of $A_2$ using: true value 165g
(a)	regression adjustment
(b)	PS stratification
(c)	PS adjustment
(d)	PS matching
(e)	IPW
(f)	AIPW
*/

*First: some data manipulation: centre MATERNAL age and BW and then create quadratic terms for them
su age
gen cage=age-r(mean)
gen cage2=cage^2
su wgt0
gen cwgt0=wgt0-r(mean)
gen cwgt02=cwgt0^2


*(a) 	using regression adjustment
reg wgt3  a2, base 
file write myfile "&Crude regression        & " %7.1f (_b[a2]) "& (" %7.1f (_se[a2]) ") \\" _n
reg wgt3 a2 i.educ i.smoke i.allergy cage cage2 i.location , base
file write myfile "&Regression adjustment {\tiny(simple) }    & " %7.1f (_b[a2]) "& (" %7.1f (_se[a2]) ") \\" _n  
teffects ra (wgt3 i.educ i.smoke i.allergy cage cage2 i.location) (a2)
matrix b = e(b)
matrix V =e(V)
file write myfile "&Regression adjustment {\tiny(with interactions) }    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

/*PROPENSITY SCORE MODEL*/
cap drop ps
logit a2 i.location i.educ i.smoke cage cage2 i.allergy,  asis	
predict ps
*examine overlap
cap drop x_* y_*
kdensity ps if a2==0, generate(x_ps0 y_ps0) kernel(triangle) 
kdensity ps if a2==1, generate(x_ps1 y_ps1) kernel(triangle) 
label var y_ps0 "a2=0"
label var y_ps1 "a2=1"
line y_ps0 x_ps0 || line y_ps1 x_ps1, title("Density distribution of the predicted propensity score") ///
      subtitle("By exposure status") ytitle(Density) xtitle(Propensity score) ylabel(,angle(h))  name(plot1,replace)
	  */// 	  xlabel(0(.2)1)  xscale(range(0 1))
graph export Figures\ps_A2.pdf,replace

*examine balance checks
teffects ipw (wgt3) (a2 i.location i.educ i.smoke cage cage2 i.allergy, logit) 
teffects overlap,
tebalance summarize, baseline
tebalance summarize, 
ex

*tebalance box, 
*tebalance density hba1c, name(bal_hb,replace)
*tebalance overid,bconly nolog

cap drop weight
gen weight=(a2/ps)+((1-a2)/(1-ps))
su ps weight
/*
foreach v of varlist age {	
tabstat `v',by(a2)
cap drop `v'0 `v'1
gen `v'1=a2*weight*`v'
gen `v'0=(1-a2)*weight*`v'
su `v' `v'0 `v'1
}
foreach v of varlist smoke allergy   {	
ta `v' a2,chi2 col nokey
cap drop `v'0 `v'1
gen `v'1=a2*weight*`v'
gen `v'0=(1-a2)*weight*`v'
su `v' `v'0 `v'1
}
foreach v of varlist location {
ta `v' a2,chi2 col nokey
  foreach j of numlist 1/4{
	cap drop `v'0`j' `v'1`j'
	gen `v'1`j'=a2*weight*(`v'==`j')
	gen `v'0`j'=(1-a2)*weight*(`v'==`j')
	}
	su `v' `v'0* `v'1*
}
foreach v of varlist  educ {
ta `v' a2,chi2 nokey col
  foreach j of numlist 1/3{
	cap drop `v'0`j' `v'1`j'
	gen `v'1`j'=a2*weight*(`v'==`j')
	gen `v'0`j'=(1-a2)*weight*(`v'==`j')
	}
	su `v' `v'0* `v'1*
}
*/
*(b)	using PS stratification*by stratification (with 6 strata)
*with bootstrap

*check cutpoints:
capture drop ps 
capture drop psblock
logit a2 i.location i.educ i.smoke cage cage2 i.allergy, asis	
qui predict ps
qui xtile psblock=ps, nq(6)
ta psblock

capture program drop ACE_psblock
program define ACE_psblock, rclass
preserve
capture drop ps 
capture drop psblock
logit a2 i.location i.educ i.smoke cage cage2 i.allergy, asis	
qui predict ps
qui xtile psblock=ps, nq(6)
qui regress wgt3  a2  if psblock==1
local ACE_1=_b[a2 ]
forvalues i=2(1)6 {
	qui  regress wgt3  a2  if psblock==`i'
	local ACE_`i'=_b[a2 ]
}
local ACE=(`ACE_1'+`ACE_2'+`ACE_3'+`ACE_4'+`ACE_5'+`ACE_6' )/6
drop ps psblock
scalar ACE=`ACE'
restore
return scalar ACE=`ACE'
end
bootstrap r(ACE), reps(1000): ACE_psblock
matrix b = e(b)
matrix se =e(se)
file write myfile "& PS stratification$^*$ {\tiny(6 strata)}&     " %7.1f (b[1,1]) "& (" %7.1f (se[1,1]) ") \\" _n  




*(c)	using PS adjustment
*with bootstrap
capture program drop ACE_regress
program define ACE_regress, rclass
preserve
capture drop ps 
logit a2 i.location i.educ i.smoke cage cage2 i.allergy, asis	
qui predict ps
qui regress wgt3  a2  ps
local ACE=_b[a2 ]
return scalar ACE=`ACE'
restore
end
bootstrap r(ACE), reps(1000): ACE_regress
matrix b = e(b)
matrix se =e(se)

file write myfile "& Regression with PS $^*$ &     " %7.1f (b[1,1]) "& (" %7.1f (se[1,1]) ") \\" _n  

*(d)	using PS matching
*by matching on 1 nearest neighbour (caliper=0.1)
*robust SE
teffects psmatch (wgt3) (a2  i.location i.educ cage cage2 i.smoke i.allergy, logit), caliper(.1) nneighbor(1) vce(robust, nn(2))
teffects overlap,
tebalance summarize, baseline
tebalance summarize, 
tebalance box, 

*assuming homoskedastic conditional outcome variance
teffects psmatch (wgt3) (a2  i.location i.educ cage cage2 i.smoke i.allergy, logit), caliper(.1) nneighbor(1) vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(1 match) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

	
*by matching on 3 nearest neighbours (caliper=0.1)
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), caliper(.1) nneighbor(3)vce(robust, nn(4))
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), caliper(.1) nneighbor(3)vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(3 matches) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

	
*by matching on 1 nearest neighbour (caliper=0.2)
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02, logit), caliper(.08) nneighbor(1) vce(iid)
*teffects overlap, ptlevel(1)

*(e)	using IPW
teffects ipw (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS IPW   & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

*(f)	using AIPW
teffects aipw (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS DR IPW    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

*IV
 ivregress 2sls wgt3 (a2=a1)
file write myfile "IV & " %7.1f (_b[a2]) "& (" %7.1f (_se[a2]) ") \\" _n  



 *********************  ATT   ********************************************************************************************
 
************************************************************************************************************************************
*----------------------------*
*For a2: the ATT
*****
file write myfile "Results for Appendix 2  ATT: " _n

reg wgt3 a2 i.educ i.smoke i.allergy cage cage2 i.location  if (a2==1), base
file write myfile "&Regression adjustment {\tiny(simple) }    & " %7.1f (_b[a2]) "& (" %7.1f (_se[a2]) ") \\" _n  

*regression adjustment 
teffects ra (wgt3 i.educ i.smoke i.allergy cage cage2 i.location ) (a2), atet
matrix b = e(b)
matrix V =e(V)
file write myfile "&Regression adjustment {\tiny(with interactions) }  & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  


*stratified by PS cats (with 6 strata)
capture program drop ATT_psblock
program define ATT_psblock, rclass
preserve
capture drop ps 
capture drop psblock
logit a2 i.location i.educ cage cage2 i.smoke i.allergy, asis	
qui predict ps
qui xtile psblock=ps, nq(6)
qui regress wgt3 a2  if psblock==1
su if psblock==1 & a2==1
local n1=r(N)
local ATT_1=_b[a2]
forvalues i=2(1)6 {
	qui regress wgt3 a2  if psblock==`i'
	su if psblock==`i' & a2==1
	local n`i'=r(N)
	local ATT_`i'=_b[a2]
}
local ATT=(`n1'*`ATT_1'+`n2'*`ATT_2'+`n3'*`ATT_3'+`n4'*`ATT_4'+`n5'*`ATT_5'+`n6'*`ATT_6' )/(`n1'+`n2'+`n3'+`n4'+`n5'+`n6')
drop ps psblock
return scalar ATT=`ATT'
restore
end
bootstrap r(ATT), reps(1000): ATT_psblock
matrix b = e(b)
matrix se =e(se)
file write myfile "& PS stratification$^*$ {\tiny(6 strata)}&     " %7.1f (b[1,1]) "& (" %7.1f (se[1,1]) ") \\" _n  


*by matching
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), ///
          nneighbor(1) caliper(.1) atet vce(robust, nn(2))
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), ///
          nneighbor(1) caliper(.1) atet vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(1 match) }   & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), ///
          nneighbor(3) caliper(.1) atet vce(robust, nn(4))
teffects psmatch (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), ///
          nneighbor(3) caliper(.1) atet vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(3 matches) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

		  
*by IPW
teffects ipw (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy, logit), ///
           atet vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS IPW    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

*teffects aipw (wgt3) (a2 i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02, logit), ///
 *          atet
*matrix b = e(b)
*matrix V =e(V)



file close myfile


ex





