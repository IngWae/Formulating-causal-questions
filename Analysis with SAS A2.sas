*load the data; 
proc import out= probit datafile = "C:\Users\slecessie\Dropbox\Stratos-causality\simulations\CURRENT VERSION SIMULATION + DATA ANALYSIS\Data\PROBITsim2018_v12.dta" replace;
run;

* 2. DATA DESCRIPTION;
proc contents data=probit;
run;
proc means data =probit;
run;
proc freq data=probit;
tables a1 a2 a3 a4 location educ smoke cesarean sex allergy ;
run;

*Uptake of the program and number of breast feeding by allocation;
proc freq data=probit;
tables a1*a2 a1*a3 /nocol nopercent;
run;

* 3 THE EFFECT OF RANDOMISATION (A1);
* Estimate the average mean dierence in infant's weight at 3 months between those assigned (A1=1)and
those not assigned to the program (A1=0);
proc reg data=probit  PLOTS=NONE;
  model wgt3 = a1;
  ods output ParameterEstimates = A1effect;
run;

* 4 THE EFFECT OF UPTAKE OF THE BFE PROGRAMME (A2);
*centre maternal age and birthweight, to create quadratic term;
proc means data = probit mean;
  var age wgt0;
run;
data probit; set probit;
   cage = age - 24.2666628;
   cage2=cage**2;
   cwgt0 = wgt0 - 3061.01;
   cwgt02=cwgt0**2;
run;

*crude estimate (which is confounded);
 proc reg data=probit PLOTS=NONE;
  model wgt3 = a2/clb;
    ods output ParameterEstimates = A2_ATE_crude;
 run;

* Simple adjustment in linear regression (GLM because of categorical variables);
proc glm data=probit PLOTS=NONE;
class educ location ;
model wgt3 = a2 educ smoke allergy cage cage2 location /solution clparm; 
 ods output ParameterEstimates = A2_ATE_adjsimple;
run;

* OUTCOME REGRESSION, fit two potential outcomes models for A2=1 and A2=0;
* Then average across the sample;
* This can be done in SAS/STAT 14.2 by proc causaltrt;
proc causaltrt data=probit method=regadj ;
class educ location a2;
psmodel a2(ref=FIRST);
model wgt3 = educ smoke allergy cage cage2 location ; 
 ods output CausalEffects = A2_ATE_adjinter;
run;

* Propensity based methods;
* Build a propensity model using logistic regression with a2 as dependent and 
  location educ smoke cage cage2 allergy as independent variables;
PROC LOGISTIC data=probit descending;
  class educ location ;
  model a2 = location educ smoke cage cage2 allergy ; 
  OUTPUT OUT=ps1 PREDICTED=ps1;
  run;

*examine overlap;
data ps1; set ps1;
  if a2 = 1 then ps1_1 = ps1;
  if a3 = 0 then ps1_0 = ps1;
 run;
proc sgplot data=ps1;
  density ps1_1 / type=kernel  LEGENDLABEL= "a2=1";
  density ps1_0 / type=kernel LEGENDLABEL= "a2=0";
run;


* b.assessing balance and propensity stratification in sas/stat 14.2;
proc psmatch data =probit region=allobs;
  class educ location a2;
  psmodel a2 (Treated='yes')= location educ smoke cage cage2 allergy ; 
  strata nstrata=6 KEY=NONE; 
  assess ps var=(smoke cage cage2 allergy ) / varinfo plots=(boxplot barchart);
  output out(obs=all)=ps_strata;
run;
proc print data=ps_strata(obs=10);
run;
  
* calculate effect per stratum of propensity score and marginalize the effects;
* se's not corrected for estimation of propensity score;
proc glm data = ps_strata plots=none;
      class a2 _strata_;
      model wgt3=a2 _strata_ a2*_strata_;
 	  estimate 'diff treated vs not treated' a2 -1 1;
	  ODS output Estimates= A2_ATE_prop_strat;
run;
   
* propensity score adjustment;
proc reg data=ps1 PLOTS=NONE;
  model wgt3 = a2 ps1/clb;
 run;



*Propensity weighting;
* Note that the standard errors calculated here also take into account that weights are estimated;
proc causaltrt data=probit ;
class educ location a2;
psmodel a2(ref=FIRST) = educ smoke allergy cage cage2 location ; 
model wgt3;
run;
* Bootstrap;
proc causaltrt data=probit ;
bootstrap seed = 462;
class educ location a2;
psmodel a2(ref=FIRST) = educ smoke allergy cage cage2 location ; 
model wgt3;
 ods output CausalEffects = A2_ATE_IPW;
run;



* Double robust;
proc causaltrt data=probit method=AIPW;
class educ location a2;
psmodel a2(ref=FIRST) = educ smoke allergy cage cage2 location ; 
model wgt3=educ smoke allergy cage cage2 location ; ;
run;
*bootstrap;
proc causaltrt data=probit method=AIPW;
bootstrap seed = 462;
class educ location a2;
psmodel a2(ref=FIRST) = educ smoke allergy cage cage2 location ; 
model wgt3=educ smoke allergy cage cage2 location ;
 ods output CausalEffects = A2_ATE_DR; 
run;

* IV;
proc syslin data= probit 2sls first ;
   endogenous  a2;
   instruments a1;
   model wgt3= a2;  
   ODS output ParameterEstimates= A2_IV;
run;

*************************************************************************************************
*
*  OUTPUT ATE RESULTS                                                                        *******;
*
************************************************************************************************;

proc print data=A1effect;
run;
data A1effect; set A1effect;
if Variable = "a1";
Estimand = "ATE A1"; 
keep Estimand Estimate StdErr ;
run;
proc print data= A2_ATE_crude;
 run;
data A2_ATE_crude; set A2_ATE_crude;
if Variable = "a2";
Estimand = "A2: Crude regression                  "; 
keep Estimand Estimate StdErr ;
run;
proc print data= A2_ATE_adjsimple;
run;
data A2_ATE_adjsimple; set A2_ATE_adjsimple;
if Parameter = "a2";
Estimand = "Regression adjustment (simple)"; 
keep Estimand Estimate StdErr ;
run;

proc print data= A2_ATE_adjinter;
run;
data A2_ATE_adjinter; set A2_ATE_adjinter;
if Parameter = "ATE";
Estimand ="Regression adjustment (with interactions)";
keep Estimand Estimate StdErr ;
run;

proc print data=A2_ATE_prop_strat;
run;
data A2_ATE_prop_strat; set A2_ATE_prop_strat;
Estimand ="PS stratification (6 strata)";
keep Estimand Estimate StdErr;
run;

proc print data= A2_ATE_IPW;
run;
data A2_ATE_IPW; set A2_ATE_IPW;
if Parameter = "ATE";
Estimand ="PS IPW";
keep Estimand Estimate StdErr BTStdErr;
run;


proc print data= A2_ATE_DR; 
run;
data A2_ATE_DR; set A2_ATE_DR;
if Parameter = "ATE";
Estimand ="PS DR";
keep Estimand Estimate StdErr BTStdErr;
run;
proc print data= A2_IV;
run;
data A2_IV; set A2_IV;
if Variable = "a2";
Estimand = "IV"; 
keep Estimand Estimate StdErr ;
run;



data ate; set A2_ATE_crude A2_ATE_adjsimple A2_ATE_adjinter A2_ATE_prop_strat A2_ATE_IPW A2_ATE_DR A2_IV;
run;
proc print; var  Estimand Estimate StdErr BTStdErr; run;


**************************  ATT  ***************************************;

* OUTCOME REGRESSION, fit two potential outcomes models for A2=1 and A2=0;
* Then average across the sample;
* This can be done in SAS/STAT 14.2 by proc causaltrt;
proc causaltrt data=probit method=regadj att;
class educ location a2;
psmodel a2(ref=FIRST);
model wgt3 = educ smoke allergy cage cage2 location ; 
ods output CausalEffects = A2_ATT_adjinter;
run;


* propensity matching;
proc psmatch data =probit region=allobs;
  class educ location a2;
  psmodel a2 (Treated='yes')= location educ smoke cage cage2 allergy ; 
  match method=greedy(k=1) caliper=0.1;
  assess ps var=(smoke cage cage2 allergy ) / varinfo plots=(boxplot barchart) weight=none;
  output out(obs=match)=ps_match matchid= MatchID;
run;
* ps_matching matches treated to untreated. The treated with a control match are in dataset ps_match ;
proc print data=ps_match(obs=10);
run;

* effect in matched subset;
proc reg data=ps_match PLOTS=NONE;
  model wgt3 = a2 /clb;
   ods output ParameterEstimates = A2_ATT_match1;
 run;

* propensity matching, 3 controls to one treated;
proc psmatch data =probit region=allobs;
  class educ location a2;
  psmodel a2 (Treated='yes')= location educ smoke cage cage2 allergy ; 
  match method=greedy(k=3) caliper=0.1;
  assess ps var=(smoke cage cage2 allergy ) / varinfo plots=(boxplot barchart) weight=none;
  output out(obs=match)=ps_match matchid= MatchID;
run;
* ps_matching matches treated to untreated. The treated with a control match are in dataset ps_match ;
proc print data=ps_match(obs=10);
run;

* effect in matched subset;
proc reg data=ps_match PLOTS=NONE;
 weight _matchwgt_;
  model wgt3 = a2 /clb;
   ods output ParameterEstimates = A2_ATT_match3;
run;


*Propensity weighting;
* Note that the standard errors calculated here also take into account that weights are estimated;
proc causaltrt data=probit att;
class educ location a2;
bootstrap seed = 462;
psmodel a2(ref=FIRST) = educ smoke allergy cage cage2 location ; 
model wgt3;
ods output CausalEffects = A2_ATT_IPW;
run;


*************************************************************************************************
*
*  OUTPUT ATT RESULTS                                                                        *******;
*
************************************************************************************************;

proc print data= A2_ATT_adjinter;
run;
data A2_ATT_adjinter; set A2_ATT_adjinter;
if Parameter = "ATT";
Estimand ="Regression adjustment (with interactions)";
keep Estimand Estimate StdErr ;
run;

*proc print data=A2_ATT_prop_strat;
*run;
*data A2_ATT_prop_strat; set A2_ATT_prop_strat;
*Estimand ="PS stratification (6 strata)";
*keep Estimand Estimate StdErr;
*run;
data A2_ATT_match1; set A2_ATT_match1;
if Variable = "a2";
Estimand = "PS matching (1 match)          "; 
keep Estimand Estimate StdErr ;
run;
data A2_ATT_match3; set A2_ATT_match3;
if Variable = "a2";
Estimand = "PS matching (3 matches)          "; 
keep Estimand Estimate StdErr ;
run;

proc print data= A2_ATT_IPW;
run;
data A2_ATT_IPW; set A2_ATT_IPW;
if Parameter = "ATT";
Estimand ="PS IPW";
keep Estimand Estimate StdErr BTStdErr;
run;


data ATT; set  A2_ATT_adjinter A2_ATT_match1 A2_ATT_match3 A2_ATT_IPW;
run;
proc print; var  Estimand Estimate StdErr BTStdErr; run;

* All results;
proc print data=ATE noobs; title"results A2 ATE"; 
var  Estimand Estimate StdErr BTStdErr; run;
proc print data=ATT noobs; title"results A2 ATT";  
var  Estimand Estimate StdErr BTStdErr; run;



