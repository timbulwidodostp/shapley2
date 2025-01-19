capture program drop shapley2 
*!version 1.7 2024-03-02
program define shapley2 , eclass 
version 9.2
syntax [anything] , stat(str) [Command(str) Indepvars(str) Depvars(str) GRoup(str) MEMory FORCE Noisily]

// CHECK FIRST THAT THE USER DOESN'T USE TWICE THE shapley2 command in a row
capture confirm matrix e(shapley)
if(_rc==0){
	di as error "You must not use twice shapley2 in a row. Restore first your regression (See {help estimates store:estimates store} for help)"
	exit
}
qui{
	// create local macros for temporary files
	tempfile orgdb temp
	// store a copy of the original data
	save `orgdb'

	est store myreg

	local full=e(`stat')


	if("`command'"==""){
		local command=e(cmd)
	}
	if("`depvar'"==""){
		local depvar=e(depvar)
	}
	if("`indepvars'"==""){
		local indepvarstemp:colnames e(b)		// load indepvars
		
		// control if all is ok (this eliminates some additional columns like _cons
		local indepvars=""
		foreach var of local indepvarstemp{
			capture confirm variable `var'
				if(_rc==0){
					local indepvars="`indepvars' `var'"	
				}
			}
	}

	// if the original command had an if condition, use this to recover it
	local if_condition: regexmatch "`e(cmdline)'" " (if [^,]+),?"
	if `if_condition'!=0{
		// there was an if condition used, so recover it
		local if_condition: regexs 1
	}
	else{
		// store an empty string which will not impact the estimation commands
		local if_condition ""
	}

	
	// CHECK THAT NO FACTOR VARIABLES ARE USED
	local cmdline=e(cmdline)
	local cmdline=subinstr(lower("`cmdline'"),lower("`command'"),"",1)
	local cmdline=subinstr(lower("`cmdline'"),lower("`depvar'"),"",1)
	
	if(strpos("`cmdline'","i.")>0 | strpos("`cmdline'","l.")>0 | strpos("`cmdline'","f.")>0 | strpos("`cmdline'","#")>0){
		noisily display as error "Factor variables are not supported" _n "Please create the variables manually"
		use `orgdb', clear
		exit
	}
		
	if("`group'"!=""){ // this is the algorithm for the group specific shapley value
		gl stop=0
		local g=1
		
		tokenize "`group'", parse(",")
		
		while($stop==0){
			local group`g' `1'
			
			macro shift
			if("`1'"==","){ // ANOTHER "," more groups are expected
				local g=`g'+1
				macro shift
			} //end if
			else{
				gl stop=1
			} // end else
		} // end while

		//special treatment for the reghdfe command's absorb option
		// will treat the absorbed FEs as a separate group
		if "`command'"=="reghdfe"{
			local absorbed_FEs: regexmatch "`e(cmdline)'" ",.*a[b]*[s]*[o]*[r]*[b]*\(([^\)]*)\)"
			if `absorbed_FEs'!=0{
				// there was an absorbed option provided, so recover it
				local absorbed_FEs: regexs 1
				local group0 `absorbed_FEs'
			}
			else{
				// no absorbed option provided, store empty string for absorbed FEs
				local absorbed_FEs ""
			}
		}
		else{
			// store an empty sting to simplify conditional logic below (if this local
			// macro is an empty string, then can execute simple command, without
			// any absorp options)
			local absorbed_FEs ""	
		}
		
		if(`g'<=12 & c(version)<12){ // for stata version <12, adapt matsize 
			local newmatsize=max(2^`g'+200,c(matsize))
			capture set matsize `newmatsize'
			}
		if(`g'>12){
			local runs=2^`K'
			noisily di as error "Too many groups defined (`runs' needed)"
			noisily di as error "A maximum of 12 groups is allowed"
			exit
		}
		
		preserve

			drop _all
			set obs 2

			if "`absorbed_FEs'"!=""{
				local starting_group_count 0
			}
			else{
				local starting_group_count 1
			}

			forvalues j=`starting_group_count'/`g'{
				gen _group`j'=1 in 1/1
				replace _group`j'=0 in 2/2
			}
			fillin _group*
			
			local allvars ""
			// this loop will start from 1 since group 0 involves absorbed FEs
			// which will be handled separately
			forvalues j=1/`g'{
				foreach var of local group`j'{
					gen `var'=_group`j'
					local allvars "`allvars' `var'"
				}
			}

			// create the special column to indicate inclusion of absorbed variables
			if "`absorbed_FEs'"!=""{
				gen absorbed = _group0
			}

			drop _fillin
			gen result=.
			mkmat * , matrix(combinations) 
			matrix list combinations

		restore
		

		local numcomb=rowsof(combinations)
		local numcols=colsof(combinations)

		// special case of no variables - 0 for the statistic of interest
		matrix combinations[1,`numcols']=0

		forvalues i=2/`numcomb'{
			local thisvars=""
			foreach var of local allvars{
				
				matrix mymat=combinations[`i',"`var'"]
				local test=mymat[1,1]
				
				if(`test'==1){
					
					local thisvars "`thisvars' `var'"
					
				}
			}

			if "`absorbed_FEs'"==""{
				// can assume that no absorption treatment is needed, so use the
				// approach used in shapley version 1.5
				`noisily' `command' `depvar' `thisvars' `if_condition'
			}
			else{
				// we have some absorption logic to incorporate, so need to
				// check if the command should include the absorbed FEs or not
				// NOTE that right now all options apart from absorb are ignored
				local include_absorbed = combinations[`i', "absorbed"]
				if `include_absorbed'==1{
					
					`noisily' `command' `depvar' `thisvars' `if_condition', absorb("`absorbed_FEs'")
				}
				else{
					`noisily' regress `depvar' `thisvars' `if_condition'
				}

			}

			matrix combinations[`i',`numcols']=e(`stat')
				
		}

		matrix list combinations
		preserve
		drop _all
		matrix list combinations
		svmat combinations,names(col) 

	}
	else{ 
		// no group variable, hence use all indepvars individually
		local K=wordcount("`indepvars'")

		//special treatment for the reghdfe command's absorb option
		// will treat the absorbed FEs as a separate variable
		if "`command'"=="reghdfe"{
			local absorbed_FEs: regexmatch "`e(cmdline)'" ",.*a[b]*[s]*[o]*[r]*[b]*\(([^\)]*)\)"
			if `absorbed_FEs'!=0{
				// there was an absorbed option provided, so recover it
				local absorbed_FEs: regexs 1
				local K=`K'+1
			}
			else{
				// no absorbed option provided, store empty string for absorbed FEs
				local absorbed_FEs ""
			}
		}
		else{
			// store an empty sting to simplify conditional logic below (if this local
			// macro is an empty string, then can execute simple command, without
			// any absorp options)
			local absorbed_FEs ""	
		}


		if(`K'>20 & "`force'"==""){
			local runs=2^`K'
			noisily di as error "Too many independent variables (`runs' needed)"
			noisily di as error "If you really want to proceed, use the option 'force'"
			exit
		}
		if(`K'<=12 & c(version)<12){
			local newmatsize=max(2^`K'+200,c(matsize))
			capture set matsize `newmatsize'
			
		}

		if(2^`K'<=c(matsize)){

			preserve

				drop _all
				set obs 2

				foreach var of local indepvars {
					gen `var'=1 in 1/1
					replace `var'=0 in 2/2
				}

				if "`absorbed_FEs'"!=""{
					gen absorbed=1 in 1/1
					replace absorbed=0 in 2/2
					fillin `indepvars' absorbed
				}
				else{
					// fill in using only indepvars
					fillin `indepvars'
				}
					
				drop _fillin
				gen result=.

				mkmat *, matrix(combinations) 
				matrix list combinations

			restore

			local numcomb=rowsof(combinations)
			local numcols=colsof(combinations)

			//di as error "I have to perform `numcomb' regressions"
			matrix combinations[1,`numcols']=0
			forvalues i=2/`numcomb'{
				local thisvars=""
				foreach var of local indepvars{
					matrix mymat=combinations[`i',"`var'"]
					local test=mymat[1,1]
					
					if(`test'==1){
						local thisvars "`thisvars' `var'"
						}
				}
		
				if "`absorbed_FEs'"==""{
					// can assume that no absorption treatment is needed, so use the
					// approach used in shapley version 1.5
					`noisily' `command' `depvar' `thisvars' `if_condition'
				}
				else{
					// we have some absorption logic to incorporate, so need to
					// check if the command should include the absorbed FEs or not
					// NOTE that right now all options apart from absorb are ignored
					local include_absorbed = combinations[`i', "absorbed"]
					if `include_absorbed'==1{
						
						`noisily' `command' `depvar' `thisvars' `if_condition', absorb("`absorbed_FEs'")
					}
					else{
						`noisily' regress `depvar' `thisvars' `if_condition'
					}

				}
				matrix combinations[`i',`numcols']=e(`stat')
			}

			preserve
			drop _all
			matrix list combinations
			svmat combinations,names(col)

		}
		else{
			// if the matsize is to big

			if "`absorbed_FEs'"!=""{
				di as error "Slow algorithm with absorbption is not supported."
				exit
			}

			if("`mem'"=="mem"){
				clear
				capture set mem 5000m
				while(_rc!=0){
				capture set mem `i'm
				local i=round((`i')*0.9)
				}
				use `orgdb'
			}


			di as error "Slow algorithm chosen. Try to increase matsize to enable the faster algorithm"
			drop _all
			set obs 2

			foreach var of local indepvars {
				gen `var'=1 in 1/1
				replace `var'=0 in 2/2
				}
			compress	
			fillin `indepvars'
			drop _fillin
			gen result=.
			
			local numcomb=_N
			
			di "`numcomb' combinations!"
			qui:replace result=0 in 1/1
			forvalues i=2/`numcomb'{
				local thisvars=""
				foreach var of local indepvars{
					local test=`var' in `i'/`i'
					if(`test'==1){
						local thisvars "`thisvars' `var'"
					}
				}
			
				//di "`thisvars'"
				preserve

					use `orgdb', clear
					di "`command' `depvar' `thisvars' `if_condition'"
					qui: `command' `depvar' `thisvars' `if_condition'

				restore
				
				qui: replace result=e(`stat') in `i'/`i'
			
			}
		}
	}

	/* Start computing the shapley value*/

	if("`group'"!=""){
		foreach var of varlist _group*{
			local IVlist "`IVlist' `var'"
		}
		egen t=rowtotal(_group*)
		sum t
		local Kgroup=r(max)
		replace t=t-1
		
		gen _weight = round(exp(lnfactorial(abs(t))),1) * round(exp(lnfactorial(`Kgroup'-abs(t)-1)),1)
		drop t
		
		keep _group* _weight result
		save `temp', replace
		matrix newshapley=[.]
		foreach var of local IVlist{
			local i=subinstr("`IVlist'","`var'","",1)
			reshape wide result _weight, i(`i') j(`var')
			gen _diff = result1-result0
			sum _diff [iweight = _weight1]
			use `temp',clear
			
			matrix newshapley = (newshapley \ r(mean))
		}
		matrix newshapley = newshapley[2...,1]
		matrix shapley=newshapley
		matrix shapley_rel=shapley/`full'
		
		// end if group
	} 
	else{

		if "`absorbed_FEs'"==""{
			// no absorption needed, so use the original logic from v1.5
			egen t=rowtotal(`indepvars')
		}
		else{
			egen t=rowtotal(`indepvars' absorbed)
		}
		replace t=t-1
		gen _weight = round(exp(lnfactorial(abs(t))),1) * round(exp(lnfactorial(`K'-abs(t)-1)),1)
		drop t

		save `temp', replace

		if "`absorbed_FEs'"==""{
			// no absorption, so use the original logic from v1.5
			matrix newshapley=[.]
			foreach var of local indepvars{
				local i=subinstr("`indepvars'","`var'","",1)
				reshape wide result _weight, i(`i') j(`var')
				gen _diff = result1-result0
				sum _diff [iweight = _weight1]
				use `temp',clear
				
				matrix newshapley = (newshapley \ r(mean))
			}
		}
		else{

			matrix newshapley=[.]
			foreach var of local indepvars{
				local i=subinstr("`indepvars'","`var'","",1)
				reshape wide result _weight, i(`i' absorbed) j(`var')
				gen _diff = result1-result0
				sum _diff [iweight = _weight1]
				use `temp',clear
				
				matrix newshapley = (newshapley \ r(mean))
			}

			// add the last row for absorbed FEs
			reshape wide result _weight, i(`indepvars') j(absorbed)
			gen _diff = result1-result0
			sum _diff [iweight = _weight1]
			use `temp',clear
			matrix newshapley = (newshapley \ r(mean))
		}

		matrix newshapley = newshapley[2...,1]

		matrix shapley=newshapley
		matrix shapley_rel=shapley/`full'

	}

	// GENERATE THE NORMALIZED VERSION
	matrix result=(shapley'\shapley_rel')
	restore

} // end quietly


// START OUTPUT
//di as text  "---------1---------2---------3---------4---------5---------6---------7---------8---------9---------10--------+"
di as text "Factor" _col(12) "{c |}" " Shapley value " _col(23) "{c |}  Per cent " //_col(40) "{c |} Shapley value" _col(45) "{c |}   Per cent  "
di as text  _col(12) "{c |}" "  (estimate)   " _col(23) "{c |} (estimate)" //_col(40) "{c |} (normalized) " _col(45) "{c |} (normalized)"
di as text "{hline 11}{c +}{hline 15}{c +}{hline 11}{c +}" //{hline 14}{c +}{hline 13}"
local i=0
if("`group'"!=""){
	// if absorbed FEs are provided, report their contribution
	if "`absorbed_FEs'"!=""{
		// note that because group0 contains the absorbed FEs, the other groups have to be offset by 1
		local i=1
		forvalues j=1/`g'{
			local i=`i'+1
			local varname="Group `j'" 
			di as text "`varname'" _col(12) "{c |}" as result %6.5f _col(15) el(result,1,`i') as text _col(28) "{c |}" _col(31) as result %6.2f 100*el(result,2,`i') as text " %" ///
			_col(40) "{c |}" //as result %6.5f _col(42) el(result,3,`i') as text _col(55) "{c |}" _col(57) as result %6.2f 100*el(result,4,`i') as text " %"
		}		
		//now report the contribution of the absorbed FEs, they should be in the first col since they correspond to group0
		local varname="Absorb. FE" 
		di as text "`varname'" _col(12) "{c |}" as result %6.5f _col(15) el(result, 1, 1) as text _col(28) "{c |}" _col(31) as result %6.2f 100*el(result, 2, 1) as text " %" ///
		_col(40) "{c |}" //as result %6.5f _col(42) el(result,3,`i') as text _col(55) "{c |}" _col(57) as result %6.2f 100*el(result,4,`i') as text " %"

	}
	else{
		// this is the same reporting as in v1.5
		forvalues j=1/`g'{
			local i=`i'+1
			local varname="Group `j'" 
			di as text "`varname'" _col(12) "{c |}" as result %6.5f _col(15) el(result,1,`i') as text _col(28) "{c |}" _col(31) as result %6.2f 100*el(result,2,`i') as text " %" ///
			_col(40) "{c |}" //as result %6.5f _col(42) el(result,3,`i') as text _col(55) "{c |}" _col(57) as result %6.2f 100*el(result,4,`i') as text " %"
		}
	}

}
else{
	foreach var of local indepvars{
		local i=`i'+1
		local varname=abbrev("`var'",10)
		di as text "`varname'" _col(12) "{c |}" as result %6.5f _col(15) el(result,1,`i') as text _col(28) "{c |}" _col(31) as result %6.2f 100*el(result,2,`i') as text " %" ///
		_col(40) "{c |}" //as result %6.5f _col(42) el(result,3,`i') as text _col(55) "{c |}" _col(57) as result %6.2f 100*el(result,4,`i') as text " %"
	}

	if "`absorbed_FEs'"!=""{
		// add the contribution of the absorbed FEs
		local i=`i'+1
		local varname=abbrev("Absorb. FE",10)
		di as text "`varname'" _col(12) "{c |}" as result %6.5f _col(15) el(result,1,`i') as text _col(28) "{c |}" _col(31) as result %6.2f 100*el(result,2,`i') as text " %" ///
		_col(40) "{c |}" //as result %6.5f _col(42) el(result,3,`i') as text _col(55) "{c |}" _col(57) as result %6.2f 100*el(result,4,`i') as text " %"

	}
}
di as text "{hline 11}{c +}{hline 15}{c +}{hline 11}{c +}"
//di as text "Residual" _col(12) "{c |}" as result %6.5f _col(15) `full'-`sum' as text _col(28) "{c |}" _col(31) as result %6.2f 100*(1-`sum'/`full') as text " %" _col(40) "{c |}" _col(55) "{c |}"
//di as text "{hline 11}{c +}{hline 15}{c +}{hline 11}{c +}{hline 14}{c +}{hline 13}"
di as text "TOTAL" _col(12) "{c |}" as result %6.5f _col(15) `full' as text _col(28) "{c |}" _col(31) as result %6.2f 100 as text " %" ///
				   _col(40) "{c |}" //as result %6.5f _col(42) `full' as text _col(55) "{c |}" _col(57) as result %6.2f 100 as text " %"
di as text "{hline 11}{c +}{hline 15}{c +}{hline 11}{c +}"
if("`group'"!=""){ //display the groups
	di as text "Groups are:"
	forvalues j=1/`g'{
		di as text "Group `j':" _col(10) as result "`group`j''"
	}
	if "`absorbed_FEs'"!=""{
		di as text "Group absorbed FEs:" _col(10) as result "`group0'"
	}
}

quietly{
	use `orgdb', clear
	est restore myreg
	est drop myreg
	ereturn matrix shapley shapley
	ereturn matrix shapley_rel shapley_rel
	ereturn local estat_cmd="shapley2"
}
end
