*load the data; 
proc import out= probit datafile = "C:\Users\slecessie\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS\Data\probitsim2018_v12.dta" replace;
run;

* 2. DATA DESCRIPTION;
proc contents data=probit;
run;
proc means data =probit;
run;
proc freq data=probit;
tables a1  a3 a4 location educ smoke cesarean sex allergy ;
run;

*Uptake of the program and number of breast feeding by allocation;
proc freq data=probit;
tables a1* a1*a3 /nocol nopercent;
run;

* 3 THE EFFECT OF RANDOMISATION (A1);
* Estimate the average mean dierence in infant's weight at 3 months between those assigned (A1=1)and
those not assigned to the program (A1=0);
proc reg data=probit  PLOTS=NONE;
  model wgt3 = a1;
run;

* Effect of A3 separately for A1=1 and A1=0;

* A1=0;
data probit0; set probit;
if a1=0;
run;

* * 4 THE EFFECT OF Breastfeeding (A3);
*

*centre maternal age and birthweight, to create quadratic term;
proc means data = probit0 mean;
  var age wgt0;
run;
data probit0; set probit0;
   cage = age - 24.2666628;
   cage2=cage**2;
   cwgt0 = (wgt0 - 3061.01)/100;
   cwgt02=cwgt0**2;
run;


*crude estimate (which is confounded);
 proc reg data=probit0 PLOTS=NONE;
  model wgt3 = a3/clb;
   ods output ParameterEstimates = A3_ATE_crude;
 run;

* Simple adjustment in linear regression (GLM because of categorical variables);
proc glm data=probit0 PLOTS=NONE;
class educ location;
model wgt3 = a3  educ smoke allergy cage cage2 location cesarean sex cwgt0 cwgt02  /solution clparm; 
 ods output ParameterEstimates = A3_ATE_adjsimple;
run;
* OUTCOME REGRESSION, fit two potential outcomes models for A3=1 and A3=0;
* Then average across the sample;
* This can be done in SAS/STAT 14.2 by proc causaltrt;
proc causaltrt data=probit0 method=regadj ;
class educ location a3;
psmodel a3(ref=FIRST);
model wgt3 =  educ smoke allergy cage cage2 location cesarean sex cwgt0 cwgt02  ; 
ods output CausalEffects = A3_ATE_adjinter;
run;

* Propensity based methods;
* Build a propensity model using logistic regression with a3 as dependent and 
  location educ smoke cage cage2 allergy  cesarean sex cwgt0 cwgt02 as independent variables;
PROC LOGISTIC data=probit0 descending;
  class educ location ;
  model a3 =  location educ smoke cage cage2 allergy cesarean sex cwgt0 cwgt02 ; 
  OUTPUT OUT=ps1 PREDICTED=ps1;
  run;
*examine overlap;
data ps1; set ps1;
  if a3 = 1 then ps1_1 = ps1;
  if a3 = 0 then ps1_0 = ps1;
 run;
proc sgplot data=ps1;
  density ps1_1 / type=kernel  LEGENDLABEL= "a3=1";
  density ps1_0 / type=kernel LEGENDLABEL= "a3=0";
run;

* b.assessing balance and propensity stratification in sas/stat 14.2;

proc psmatch data =ps1 region=allobs;
  class educ location  a3;
  psmodel a3 (Treated='yes')=  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ;  
  strata nstrata=6 KEY=NONE; 
  assess ps var=(smoke cage cage2 allergy ) / varinfo plots=(boxplot barchart);
  output out(obs=all)=ps_strata;
run;
proc print data=ps_strata(obs=10);
run;
* calculate effect per stratum of propensity score and marginalize the effects;
* se's not corrected for estimation of propensity score;
proc glm data = ps_strata plots=none;
      class a3 _strata_;
      model wgt3=a3 _strata_ a3*_strata_;
 	  estimate 'diff treated vs not treated' a3 -1 1;
  	  ODS output Estimates= A3_ATE_prop_strat;
run;
   
* propensity score adjustment;
proc reg data=ps1 PLOTS=NONE;
  model wgt3 = a3 ps1/clb;
 run;


*Propensity weighting;
* Note that the standard errors calculated here also take into account that weights are estimated;
* Bootstrap;
proc causaltrt data=probit0 ;
bootstrap seed = 462;
class educ location a3;
psmodel a3(ref=FIRST) =  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ; 
model wgt3;
ods output CausalEffects = A3_ATE_IPW;
run;


* Double robust;
*bootstrap;
proc causaltrt data=probit0 method=AIPW;
bootstrap seed = 462;
class educ location a3;
psmodel a3(ref=FIRST) =  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ; 
model wgt3=educ smoke allergy  cesarean sex cwgt0 cwgt02 cage cage2 location ; 
ods output CausalEffects = A3_ATE_DR; 
run;
;

run;

* IV for complete data;
proc syslin data= probit 2sls first ;
   endogenous  a3;
   instruments a1;
   model wgt3= a3;  
   ODS output ParameterEstimates= _IV;
run;
*************************************************************************************************
*
*  OUTPUT ATE RESULTS                                                                        *******;
*
************************************************************************************************;

proc print data= A3_ATE_crude;
 run;
data A3_ATE_crude; set A3_ATE_crude;
if Variable = "a3";
Estimand = "A3: Crude regression                  "; 
keep Estimand Estimate StdErr ;
run;
proc print data= A3_ATE_adjsimple;
run;
data A3_ATE_adjsimple; set A3_ATE_adjsimple;
if Parameter = "a3";
Estimand = "Regression adjustment (simple)"; 
keep Estimand Estimate StdErr ;
run;

proc print data= A3_ATE_adjinter;
run;
data A3_ATE_adjinter; set A3_ATE_adjinter;
if Parameter = "ATE";
Estimand ="Regression adjustment (with interactions)";
keep Estimand Estimate StdErr ;
run;

proc print data=A3_ATE_prop_strat;
run;
data A3_ATE_prop_strat; set A3_ATE_prop_strat;
Estimand ="PS stratification (6 strata)";
keep Estimand Estimate StdErr;
run;

proc print data= A3_ATE_IPW;
run;
data A3_ATE_IPW; set A3_ATE_IPW;
if Parameter = "ATE";
Estimand ="PS IPW";
keep Estimand Estimate StdErr BTStdErr;
run;


proc print data= A3_ATE_DR; 
run;
data A3_ATE_DR; set A3_ATE_DR;
if Parameter = "ATE";
Estimand ="PS DR";
keep Estimand Estimate StdErr BTStdErr;
run;
proc print data= _IV;
run;
data _IV; set _IV;
if Variable = "a3";
Estimand = "IV"; 
keep Estimand Estimate StdErr ;
run;



data ate; set A3_ATE_crude A3_ATE_adjsimple A3_ATE_adjinter A3_ATE_prop_strat A3_ATE_IPW A3_ATE_DR _IV;
run;
proc print; var  Estimand Estimate StdErr BTStdErr; run;


**************************  ATT  ***************************************;

* OUTCOME REGRESSION, fit two potential outcomes models for a3=1 and a3=0;
* Then average across the sample;
* This can be done in SAS/STAT 14.2 by proc causaltrt;
proc causaltrt data=probit0 method=regadj att;
class educ location a3;
psmodel a3(ref=FIRST);
model wgt3 =  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ; 
ods output CausalEffects = A3_ATT_adjinter;
run;


* propensity matching;
proc psmatch data =ps1 region=allobs;
  class educ location a3;
  psmodel a3 (Treated='yes')=  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ;  
  match method=greedy(k=1) caliper=0.1;
  output out(obs=match)=ps_match matchid= MatchID;
run;
* ps_matching matches treated to untreated. The treated with a control match are in dataset ps_match ;
proc print data=ps_match(obs=10);
run;

* effect in matched subset;
proc reg data=ps_match PLOTS=NONE;
  model wgt3 = a3 /clb;
  ods output ParameterEstimates = A3_ATT_match1;
 run;

* propensity matching, 3 controls to one treated;
proc psmatch data =ps1 region=allobs;
  class educ location a3;
  psmodel a3 (Treated='yes')=  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ;  
  match method=greedy(k=3) caliper=0.1;
  output out(obs=match)=ps_match matchid= MatchID;
run;
* ps_matching matches treated to untreated. The treated with a control match are in dataset ps_match ;

* effect in matched subset;
proc reg data=ps_match PLOTS=NONE;
 weight _matchwgt_;
  model wgt3 = a3 /clb;
   ods output ParameterEstimates = A3_ATT_match3;
run;


*Propensity weighting;
* Note that the standard errors calculated here also take into account that weights are estimated;
proc causaltrt data=probit0 att;
class educ location a3;
bootstrap seed = 462;
psmodel a3(ref=FIRST) =  educ smoke allergy cesarean sex cwgt0 cwgt02 cage cage2 location ; 
model wgt3;
ods output CausalEffects = A3_ATT_IPW;
run;


*************************************************************************************************
*
*  OUTPUT ATT RESULTS                                                                        *******;
*
************************************************************************************************;

proc print data= A3_ATT_adjinter;
run;
data A3_ATT_adjinter; set A3_ATT_adjinter;
if Parameter = "ATT";
Estimand ="Regression adjustment (with interactions)";
keep Estimand Estimate StdErr ;
run;

*proc print data=_ATT_prop_strat;
*run;
*data A3_ATT_prop_strat; set A3_ATT_prop_strat;
*Estimand ="PS stratification (6 strata)";
*keep Estimand Estimate StdErr;
*run;
data A3_ATT_match1; set A3_ATT_match1;
if Variable = "a3";
Estimand = "PS matching (1 match)          "; 
keep Estimand Estimate StdErr ;
run;
data A3_ATT_match3; set A3_ATT_match3;
if Variable = "a3";
Estimand = "PS matching (3 matches)          "; 
keep Estimand Estimate StdErr ;
run;

proc print data= A3_ATT_IPW;
run;
data A3_ATT_IPW; set A3_ATT_IPW;
if Parameter = "ATT";
Estimand ="PS IPW";
keep Estimand Estimate StdErr BTStdErr;
run;


data ATT; set  A3_ATT_adjinter A3_ATT_match1 A3_ATT_match3 A3_ATT_IPW;
run;
proc print; var  Estimand Estimate StdErr BTStdErr; run;

* All results;
proc print data=ATE noobs; title"results A3 ATE in subgroup A1=0"; 
var  Estimand Estimate StdErr BTStdErr; run;
proc print data=ATT noobs; title"results A3 ATT in subgroup A1=0";  
var  Estimand Estimate StdErr BTStdErr; run;
