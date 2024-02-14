*********************************************************************************
* AUTHOR		: MAGHFIRA RAMADHANI											*
* PROJECT		: PS1															*
* COURSE		: Health Economics II											*
* DESCRIPTION	: Main Code														*
* INPUT			: NA	    													*
* OUTPUT		: .\output\table, .\output\figure, .\output\log					*
* STATA VERSION	: Stata/MP 18.0													*
*********************************************************************************

clear
version 18.0
macro drop _all
set linesize 255
set more off, permanently
capture log close
capture graph drop _all
matrix drop _all

*********************************************************************************
* Setup the profile of your machine

	* Select option to install Stata packages (list package in profile.do)
	
	global install_stata_packages 0 // Set to 1 for first time running, 0 o/w
	
	* Select option to export log
	
	global export_log 0 // Set to 1 if you want to export log, 0 o/w

	* Set the location of project directory location
	
	global path "C:\Users\mramadhani3\OneDrive - Georgia Institute of Technology\Documents\Spring-24\health-econ-ii\health-ii-ps\ps-1"
	global data_path "$path\data"
	global temp_path "$path\temp"
	global code_path "$path\code" 
	global table_path "$path\output\table" 
	global figure_path "$path\output\figure"

	* ON IAC VLAB server, you will need to uncomment this line and run this:
	*sysdir set PERSONAL \\iac.nas.gatech.edu\mramadhani3

	* Set the location of Python and R executable
	
	*global RSCRIPT_PATH "C:\Program Files\R\R-4.2.2\bin\x64\Rscript.exe"
	*Change the following line to your Conda env, and uncomment the following line for first time run
	*python set exec C:\Users\mramadhani3\AppData\Local\anaconda3\envs\economics\python.exe
	
	*python set userpath "C:\Users\mramadhani3\AppData\Local\anaconda3\Lib\site-packages" "C:\Users\mramadhani3\OneDrive - Georgia Institute of Technology\Documents\Spring-24\environmental-econ-ii\phdee-24-MR\homework-2\code"


	* Set machine profile
	
	do "$code_path\0_profile.do"

*********************************************************************************
* Q1 Generate summary statistics of brand-size pair
*********************************************************************************

	*Load datasets
	import delimited "$data_path\otc_headache.csv", clear
	
	* Generate auxiliary variable
	egen total_sales=sum(sales), by(store week)
	gen share=sales/total_sales
	gen unitprice=price/size
	gen per50price=unitprice*50
	
	est clear
	foreach b in Advil Bayer "Store brand" Tylenol {
	foreach s in 25 50 100{
	count if brand == "`b'" & size==`s'
	scalar obs=r(N)
	if obs!=0 {
	eststo: quietly estpost summarize share unitprice per50price cost if brand == "`b'" & size==`s' //Assume cost is the wholesale price
	}
	}
	}
	esttab
	
	esttab using "$table_path\summary.tex", replace label ///
		cell( mean(pattern(1 1 1 1 1 ///
							1 1 1 1 1 1) fmt(3)) ) ///
		mgroups ("Advil" "Bayer" "Store Brand" "Tylenol", pattern(1 0 0 1 0 0 1 0 1 0 0 ) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("25" "50" "100" "25" "50" "100" "50" "100" "25" "50" "100") collabels(none) ///
		coeflabels(share "Market share of sales" unitprice "Unit price" per50price "Price per 50 tab" cost "Wholesale price") nonum nonote noobs
	
*********************************************************************************
* Q2-4 Demand model estimates by OLS
*********************************************************************************

	gen ln_share = ln(share)
	encode brand, gen(brand_fac)
	
	est clear
	* Estimate the demand model by OLS using price and promotion as product characteristics.
	eststo ols1: ivreghdfe ln_share price prom 
	estadd local FE "-"
	estadd scalar R2 = e(r2_a)
	* Estimate the demand model by OLS using price and promotion as product characteristics, also including brand-size dummies.
	eststo ols2: ivreghdfe ln_share price prom , absorb(brand_fac#size)
	estadd local FE "Product"
	estadd scalar R2 = e(r2_a)
	* Estimate the demand model by OLS using price and promotion as product characteristics, also including store-brand-size (the interaction of store and product) dummies.
	eststo ols3: ivreghdfe ln_share price prom , absorb(store#brand_fac#size)
	estadd local FE "Store-Product"
	estadd scalar R2 = e(r2_a)
	
	*** Using .tex
	esttab ols1 ols2 ols3 using "$table_path\OLS.tex", label replace ///
		keep(prom price) ///
		b(3) se(3) ////
		mgroups ("OLS", pattern(1 0 0) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("(1)" "(2)" "(3)") collabels(none) nostar nonote nonum ///
		coeflabels(prom "Promotions" price "Price") ///
		scalars("FE Fixed Effect" "R2 Adjusted-R$^2$") obslast

*********************************************************************************
* Q5-6 Demand model estimates by IV
*********************************************************************************		
	*Generate Hausman instrument
	bysort brand size week: gen otherstorecount=_N
	bysort brand size week: egen hausman=sum(price)
	replace hausman=(hausman-price)/(otherstorecount-1)
	est clear
	
	eststo iv1: ivreghdfe ln_share (price=cost) prom 
	estadd local FE "-"
	estadd local IV "Wholesale Price"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a)
	scalar alpha_1=e(b)[1,1]
	
	eststo iv2: ivreghdfe ln_share  (price=cost) prom, absorb(brand_fac#size)
	estadd local FE "Product"
	estadd local IV "Wholesale Price"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a)
	scalar alpha_2=e(b)[1,1]
	
	eststo iv3: ivreghdfe ln_share  (price=cost) prom, absorb(store#brand_fac#size)
	estadd local FE "Store-Product"
	estadd local IV "Wholesale Price"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a)
	scalar alpha_3=e(b)[1,1]
	
	eststo iv4: ivreghdfe ln_share  (price=hausman) prom
	estadd local FE "-"
	estadd local IV "Hausman"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a)
	scalar alpha_4=e(b)[1,1]
	
	eststo iv5: ivreghdfe ln_share  (price=hausman) prom, absorb(brand_fac#size)
	estadd local FE "Product"
	estadd local IV "Hausman"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a)
	scalar alpha_5=e(b)[1,1]
	
	eststo iv6: ivreghdfe ln_share  (price=hausman) prom, absorb(store#brand_fac#size)
	estadd local FE "Store-Product"
	estadd local IV "Hausman"
	estadd local Fstat = string(e(widstat),"%1.0f")
	estadd scalar R2 = e(r2_a) 
	scalar alpha_6=e(b)[1,1]
	
	*** Using .tex
	esttab iv1 iv2 iv3 iv4 iv5 iv6 using "$table_path\IV.tex", label replace ///
		keep(prom price) ///
		b(3) se(3) ////
		mgroups ("IV: Wholesale Price" "IV: Hausman Instrument", pattern(1 0 0 1 0 0) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("(1)" "(2)" "(3)" "(4)" "(5)" "(6)") collabels(none) nostar nonote nonum ///
		coeflabels(prom "Promotions" price "Price") ///
		scalars("FE Fixed Effect" "Fstat Montiel-Pflueger F-statistics" "R2 Adjusted-R$^2$") obslast

*********************************************************************************
* Q8 Mean own price elasticities
*********************************************************************************	

	gen price_x_1_share=price*(1-share)
	
	* Calling estimates saved from previous model
	est clear
	foreach e of num 1/6 {
	gen own_price_elasticity_`e'= alpha_`e' * price_x_1_share
	}
	
	foreach b in Advil Bayer "Store brand" Tylenol {
	foreach s in 25 50 100 {
	count if brand == "`b'" & size==`s'
	scalar obs=r(N)
	if obs!=0 {
	eststo: quietly estpost summarize own_price_elasticity_* if brand == "`b'" & size==`s' //Assume cost is the wholesale price
	}
	}
	}
	
	esttab using "$table_path\elasticity.tex", replace label ///
		cell( mean(pattern(1 1 1 1 1 ///
							1 1 1 1 1 1) fmt(3)) ) ///
		mgroups ("Advil" "Bayer" "Store Brand" "Tylenol", pattern(1 0 0 1 0 0 1 0 1 0 0 ) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("25" "50" "100" "25" "50" "100" "50" "100" "25" "50" "100") collabels(none) ///
		coeflabels(own_price_elasticity_1 "Specification (1)" ///
		own_price_elasticity_2 "Specification (2)" own_price_elasticity_3 "Specification (3)" ///
		own_price_elasticity_4 "Specification (4)" own_price_elasticity_5 "Specification (5)" ///
		own_price_elasticity_6 "Specification (6)") nonum nonote noobs
	
	
	* Alternative format but not used
	* Generate dummy for groups
	gen own_price_elasticity= alpha_1 * price_x_1_share
	levelsof brand_fac
	foreach b in `r(levels)' {
		gen i_`b'=0
		replace i_`b'=1 if brand_fac==`b'
	}
	
	est clear
	foreach s in 25{
	eststo: xi: reg own_price_elasticity i_1 i_2 i_4 if size==`s', noconstant //Assume cost is the wholesale price
	}
	foreach s in 50 100{
	eststo: xi: reg own_price_elasticity i_* if size==`s', noconstant //Assume cost is the wholesale price
	}
	esttab
	
	esttab using "$table_path\elasticity1.tex", replace label ///
		cell( b(pattern(1 1 1) fmt(3)) ) order(i_1 i_2 i_3 i_4) /// ///
		mgroups ("Own Price Elasticity", pattern(1 0 0) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("25" "50" "100") collabels(none) ///
		coeflabels(i_1 "Advil" i_2 "Bayer" i_3 "Store Brand" i_4 "Tylenol") nonum nonote noobs nostar

*********************************************************************************
* Q9 Mean change in consumer surplus
*********************************************************************************	
	
	gen ln_1_share=ln(1-share)
	* Calling estimates saved from previous model
	est clear
	foreach e of num 1/6 {
	gen surplus_`e'= (-1/alpha_`e') * ln_1_share
	}
	
	foreach b in Advil Bayer "Store brand" Tylenol {
	foreach s in 25 50 100 {
	count if brand == "`b'" & size==`s'
	scalar obs=r(N)
	if obs!=0 {
	eststo: quietly estpost summarize surplus_* if brand == "`b'" & size==`s' //Assume cost is the wholesale price
	}
	}
	}
	
	esttab using "$table_path\surplus.tex", replace label ///
		cell( mean(pattern(1 1 1 1 1 ///
							1 1 1 1 1 1) fmt(3)) ) ///
		mgroups ("Advil" "Bayer" "Store Brand" "Tylenol", pattern(1 0 0 1 0 0 1 0 1 0 0 ) ///
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
		mtitle("25" "50" "100" "25" "50" "100" "50" "100" "25" "50" "100") collabels(none) ///
		coeflabels(surplus_1 "Specification (1)" ///
		surplus_2 "Specification (2)" surplus_3 "Specification (3)" ///
		surplus_4 "Specification (4)" surplus_5 "Specification (5)" ///
		surplus_6 "Specification (6)") nonum nonote noobs

*********************************************************************************
* End of code
if $export_log == 1{
	log close
	}
