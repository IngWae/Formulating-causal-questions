cap log close
set more off
set scheme s2mono
clear all
discard
set seed 50159
*cd "C:\Users\slecessie\Dropbox\Stratos-causality\May June 2018\Analyses with stata"
cd "C:\Users\Bianca DeStavola\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS"
*cd "C:\Users\sejjbld\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS"


*log using "Analyses with stata\Appendix2_A3_A1_0_stata.log", replace

*****
capture : file close myfile
file open myfile using "Analyses with stata\Results_for_A3_A1_0_03052020.txt", write replace

	file write myfile "Today is  $S_DATE " 

	file write myfile "Results for A3 in the subgroup randomized for control  ATE: " _n

*analyses of version 2 of the data (created in April 2018)
*use PROBITsim2018_v12.dta,clear
import delimited "PROBITsim_april2018.csv"
des 
rename wgt3observed wgt3
rename a2observed a2
rename a3observed a3


*Examine the uptake of the program and of breast feeding by allocation
tab a1 a3, row nokey

* S E L E C T I O N OF C O N T R O L    G R O U P 
keep if a1==0

*centre MATERNAL age and bweight and create quadratic term
su age
gen cage=age-r(mean)
gen cage2=cage^2
su wgt0
gen cwgt0=wgt0-r(mean)
gen cwgt02=cwgt0^2

*crude
reg wgt3  a3, base 

file write myfile "&Crude regression        & " %7.1f (_b[a3]) "& (" %7.1f (_se[a3]) ") \\" _n

*regression adjustment 

reg wgt3 a3  i.educ i.smoke i.allergy cage cage2 i.location i.cesarean i.sex cwgt0 cwgt02 i.a2, base

file write myfile "&Regression adjustment {\tiny(simple) }    & " %7.1f (_b[a3]) "& (" %7.1f (_se[a3]) ") \\" _n  

*regression adjustment that allows for all X*C interactions
teffects ra (wgt3  i.educ i.smoke i.allergy cage cage2 i.location i.cesarean i.sex cwgt0  cwgt02 i.a2) (a3)
matrix b = e(b)
matrix V =e(V)

file write myfile "&Regression adjustment {\tiny(with interactions) }    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

/*PROPENSITY SCORE MODEL*/
cap drop ps
logit a3  i.location i.educ i.smoke cage cage2 i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2 , asis	
predict ps
*examine overlap
cap drop x_* y_*
kdensity ps if a3==0, generate(x_ps0 y_ps0) kernel(triangle) 
kdensity ps if a3==1, generate(x_ps1 y_ps1) kernel(triangle) 
label var y_ps0 "a3=0"
label var y_ps1 "a3=1"
line y_ps0 x_ps0 || line y_ps1 x_ps1, title("Density distribution of the predicted propensity score") ///
      subtitle("By exposure status") ytitle(Density) xtitle(Propensity score) ylabel(,angle(h))  
*graph export Figures\ps_A30.pdf,replace

*examine balance checks
teffects ipw (wgt3) (a3  i.location i.educ i.smoke cage cage2 i.allergy i.cesarean i.sex cwgt0 cwgt02, logit) 
teffects overlap,
tebalance summarize, baseline
tebalance summarize,
ex


cap drop weight
gen weight=(a3/ps)+((1-a3)/(1-ps))
su ps weight

foreach v of varlist age wgt0 {	
tabstat `v',by(a3)
cap drop `v'0 `v'1
gen `v'1=a3*weight*`v'
gen `v'0=(1-a3)*weight*`v'
su `v' `v'0 `v'1
}
foreach v of varlist smoke allergy  cesarean sex  a2{	
ta `v' a3,chi2 col nokey
cap drop `v'0 `v'1
gen `v'1=a3*weight*`v'
gen `v'0=(1-a3)*weight*`v'
su `v' `v'0 `v'1
}
foreach v of varlist location {
ta `v' a3,chi2 col nokey
  foreach j of numlist 1/4{
	cap drop `v'0`j' `v'1`j'
	gen `v'1`j'=a3*weight*(`v'==`j')
	gen `v'0`j'=(1-a3)*weight*(`v'==`j')
	}
	su `v' `v'0* `v'1*
}
foreach v of varlist  educ {
ta `v' a3,chi2 nokey col
  foreach j of numlist 1/3{
	cap drop `v'0`j' `v'1`j'
	gen `v'1`j'=a3*weight*(`v'==`j')
	gen `v'0`j'=(1-a3)*weight*(`v'==`j')
	}
	su `v' `v'0* `v'1*
}
capture drop ps 
logit a3  i.location i.educ i.smoke cage cage2 i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, asis	
qui predict ps
qui regress wgt3  a3  ps


*by adjustment
*with bootstrap
capture program drop ACE_regress
program define ACE_regress, rclass
preserve
capture drop ps 
logit a3  i.location i.educ i.smoke cage cage2 i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, asis	
qui predict ps
qui regress wgt3  a3  ps
local ACE=_b[a3 ]
return scalar ACE=`ACE'
restore
end
bootstrap r(ACE), reps(1000): ACE_regress
matrix b = e(b)
matrix se =e(se)

file write myfile "& Regression with PS $^*$ &     " %7.1f (b[1,1]) "& (" %7.1f (se[1,1]) ") \\" _n  


*by stratification (with 6 strata)
*with bootstrap
capture program drop ACE_psblock
program define ACE_psblock, rclass
preserve
capture drop ps 
capture drop psblock
logit a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, asis	iter(100)
qui predict ps
qui xtile psblock=ps, nq(6)
qui regress wgt3  a3  if psblock==1
local ACE_1=_b[a3 ]
forvalues i=2(1)6 {
	qui  regress wgt3  a3  if psblock==`i'
	local ACE_`i'=_b[a3 ]
}
local ACE=(`ACE_1'+`ACE_2'+`ACE_3'+`ACE_4'+`ACE_5'+`ACE_6' )/6
drop ps psblock
scalar ACE=`ACE'
return scalar ACE=`ACE'
restore
end
bootstrap r(ACE), reps(1000): ACE_psblock
matrix b = e(b)
matrix se =e(se)

file write myfile "& PS stratification$^*$ {\tiny(6 strata)}&     " %7.1f (b[1,1]) "& (" %7.1f (se[1,1]) ") \\" _n  


*by matching on 1 nearest neighbour (caliper=0.1)
*robust SE
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.1) nneighbor(1) vce(robust, nn(2))
*assuming homoskedastic conditional outcome variance
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.1) nneighbor(1) vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(1 match) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

	
*by matching on 3 nearest neighbours (caliper=0.1)
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.1) nneighbor(3) vce(robust, nn(4))
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.1) nneighbor(3) vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(3 matches) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

	
*by matching on 1 nearest neighbour (caliper=0.2)
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.08) nneighbor(1) vce(robust, nn(2))
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit), caliper(.08) nneighbor(1) vce(iid)

teffects overlap, ptlevel(1)


*IPW
teffects ipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit),  
teffects ipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit),   vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS IPW   & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

* double robust
teffects aipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit),  
teffects aipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean ///
    i.sex cwgt0 cwgt02 i.a2, logit),  	vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS DR IPW    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

*IV
 *ivregress 2sls wgt3 (a3=a1)
*file write myfile "IV & " %7.1f (_b[a3]) "& (" %7.1f (_se[a3]) ") \\" _n  





 *********************  ATT   ********************************************************************************************
 
************************************************************************************************************************************
*----------------------------*
*For a3: the ATT

file write myfile _n _n "Results for a3  ATT: " _n

*regression adjustment 
teffects ra (wgt3  i.educ i.smoke i.allergy cage cage2 i.location i.cesarean i.sex cwgt0 cwgt02 i.a2) (a3), atet
matrix b = e(b)
matrix V =e(V)
file write myfile "&Regression adjustment {\tiny(with interactions) }  & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  


*-stratified by PS cats (with 6 strata)
capture program drop ATT_psblock
program define ATT_psblock, rclass
preserve
capture drop ps 
capture drop psblock
logit a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, asis iter(100)
qui predict ps
qui xtile psblock=ps, nq(6)
qui regress wgt3 a3  if psblock==1
su if psblock==1 & a3==1
local n1=r(N)
local ATT_1=_b[a3]
forvalues i=2(1)6 {
	qui regress wgt3 a3  if psblock==`i'
	su if psblock==`i' & a3==1
	local n`i'=r(N)
	local ATT_`i'=_b[a3]
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
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, logit), ///
          nneighbor(1) caliper(.1) atet vce(robust, nn(2))
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, logit), ///
          nneighbor(1) caliper(.1) atet vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(1 match) }   & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, logit), ///
          nneighbor(3) caliper(.1) atet vce(robust, nn(4))
teffects psmatch (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, logit), ///
          nneighbor(3) caliper(.1) atet vce(iid)
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS matching  {\tiny(3 matches) }     & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

		  
*by IPW
teffects ipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02 i.a2, logit), ///
           atet 	vce(boot,reps(1000))
matrix b = e(b)
matrix V =e(V)
file write myfile "&PS IPW    & " %7.1f (b[1,1]) "& (" %7.1f (sqrt(V[1,1])) ") \\" _n  

*teffects aipw (wgt3) (a3  i.location i.educ cage cage2 i.smoke i.allergy i.cesarean i.sex cwgt0 cwgt02, logit), ///
 *          atet
*matrix b = e(b)
*matrix V =e(V)



file close myfile

ex
