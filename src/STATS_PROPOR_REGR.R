#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.2

# history
# 25-Sep-2013 - original version
# 16-Oct=2013 - remove non-ascii hyphen from help text


helptext="This procedure estimates the mean and variance of
a beta-distributed proportion as a generalized linear mode as a
function of a set of predictors for the mean and variance parameters.

STATS PROPOR REGR DEPENDENT=varname
	INDEP=variable list INDEPPRE=variable list
	ID=varname
	LINK = LOGIT or PROBIT or CLOGLOG or CAUCHIT or LOG or LOGLOG
	PLINK = LOG or SQRT or IDENTIFY
	OFFSET = variable list
/OPTIONS ESTIMATOR=ML or BC or BR
	MISSING=OMIT or STOP
	MAXITER = number FSMAXITER = number FSTOL = value
	RESIDS = YES or NO COOKS=YES or NO LEVERAGE = YES or NO
	RESIDSLIN = YES or NO HALFNORMAL = YES or NO PREDOBS = YES or NO
	PRESIDTYPE = PEARSON or SWEIGHTED2 or DEVIANCE or 
		RESPONSE or WEIGHTED or SWEIGHTED
/SAVE DATASET = new dataset name
	RESIDUALS = PEARSON or SWEIGHTED2 or DEVIANCE or RESPONSE or
		WEIGHTED or SWEIGHTED
	FITTED = RESPONSEMEAN or LINK or RESPONSEVAR
	RETAIN = YES or NO
	RFILE = filespec
/HELP


The DEPENDENT and INDEP keywords are required.

Example:
STATS PROPOR REGR DEPENDENT = yield INDEP = batch temperature pressure.

DEPENDENT and INDEP specify the dependent and independent variables.
Categorical variables are converted to factors.

INDEPPRE can specify variables for the precision model.  A variable
can appear in both the INDEP and INDEPPRE list.  Categorical variables
are converted to factors.

ID specifies and ID variable that will be used if an output dataset is
created.

LINK specifies the generalized linear model link function.

PLINK specifyies the link function for the precision model.

OFFSET specifies an offset variable for the mean model.

ESTIMATOR specifies the estimator to be used.  The choices are
maximum likelihood, maximum likelihood with bias correction,
or maximum likelihood with bias reduction.

MISSING specifies whether to omit cases with missing values or
stop if any are found.

MAXITER specifies the maximum number of main iterations.  The
default is 5000.
FSMAXITER specifies the maximum number of Fischer scoring iterations.
The default is 200.

The remaining keywords on OPTIONS specify what plots to create
and the type of residuals to be used in these plots.  SWEIGHTED2
residuals are recommended but may be too time consuming with
large datasets.
See Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). 
On Beta Regression Residuals. Journal of Applied Statistics, 35(4), 407-419.
for detailed definitions.

The SAVE subcommand specifies an output dataset for residuals and
predicted values from the model and whether to keep the R workspace.
DATASET must specify a name for a new dataset that will hold these
variables.
RESIDUALS specifies the type of residuals to save.
FITTED specifies the type of fitted values.

RFILE is the name of an R workspace file to hold the model results
for future use in prediction.

RETAIN specifies whether to retain the model results in the active
R workspace for use for prediction in the current session.

STATS PROPOR REGR /HELP prints this help and does nothing else.
"

proext <- function(dep, indep, indeppre=NULL, id=NULL, link="logit", plink="log",
	offsetv=NULL, etype="ml", missingv="omit", maxiter=5000, fsmaxiter=200,
	fstol=1e-08, dataset=NULL, resids=NULL, fittedv=NULL,
	rfile=NULL, retain=FALSE,
	residso=FALSE, cooks=FALSE, leverage=FALSE, residslin=FALSE, 
	halfnormal=FALSE, predobs=FALSE, presidtype="pearson"
	) {
	
	nindep = length(indep)
	nindeppre = length(indeppre)
	etype = toupper(etype)
	allargs = as.list(environment())  # for passing external spec to other functions

    setuplocalization("STATS_PROPOR_REGR")

    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Proportional Variable Regression")
    warningsprocname = gtxt("Proportional Variable Regression: Warnings")
    omsid="STATSPROREG"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(betareg, quietly=TRUE), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","betareg"),
            dostop=TRUE)
    }
    )

	if (missingv == "omit") {  #ifelse does not work here
		naaction = na.omit
	} else {
		naaction = na.fail
	}
	if (is.null(dataset) && (!is.null(resids) || !is.null(fittedv))) {
		warns$warn(gtxt("A dataset name must be specified if residuals or fitted values are requested"),
			dostop=TRUE)
	}
	dscheck(dataset, warns)
	frml = buildformula(dep, indep, indeppre, offsetv)
	control = betareg.control(maxit=maxiter, fsmaxit=fsmaxiter, fstol=fstol)	
	allvars = c(dep, indep, indeppre, offsetv)

	dta = spssdata.GetDataFromSPSS(allvars, row.label=id, missingValueToNA=TRUE,
		factorMode="labels")
	allargs["ncases"] = nrow(dta)

	if (!is.null(offsetv) && is.factor(dta[[offsetv]])) {
	  warns$warn(gtxt("The offset variable must have a scale (continuous) measurement level"),
		dostop=TRUE)
	}
	# Can't just pass the name of the option variable, offsetv, because
	# betareg interprets that as the name of a variable in the data frame
	# A null offset just gets replicated as a column of zeros by betareg anyway.
	if (any(residso,cooks,leverage,residslin,halfnormal,predobs)) {
		keep = TRUE
	} else {
		keep = FALSE
	}
	res2 = tryCatch(
		{
		betareg(formula=frml, data=dta, na.action=naaction,
			link=link, link.phi=plink, type=etype,
			control=control, model=TRUE, y=keep, x=keep)
		},
		error=function(e) {warns$warn(e$message, dostop=TRUE)
		}
	)
	res = summary(res2)

    # print results and, if requested, save variables
    # 
	allargs["modeldate"] = date()
    displayresults(allargs, allvars, res, res2, warns)
	createdataset(res, res2, dta, allargs)

    # clean up workspace
	if (retain || !is.null(rfile)) {
        assign("stats_propor_res", res2, envir=.GlobalEnv)
		assign("stats_propor_ressum", res, envir=.GlobalEnv)
        assign("stats_propor_allargs", allargs, envir=.GlobalEnv)
        rm(res, allargs, dta)
        if (!is.null(rfile)) {
            save(stats_propor_res, stats_propor_allargs, file=rfile)
        }
    }
    if (!retain) {
        res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
    }
}


displayresults = function(info, allvars, res, res2, warns) {

    StartProcedure(gtxt("Proportional Variable Regression"), "STATSPROREG")
    
    # summary results
	# input specifications
	lbls = c(gtxt("Dependent Variable"),
		gtxt("Link Function"),
		gtxt("Precision Model Link Function"),
		gtxt("Offset"),
		gtxt("Missing Value Treatment"),
		gtxt("Estimator Type"),
		gtxt("Results Dataset"),
		gtxt("Saved Workspace"),
		gtxt("Retain Workspace"),
		gtxt("Creation Date"),
		gtxt("Convergence"),
		gtxt("Cases"),
		gtxt("Valid Cases"),
		gtxt("Log-Likelihood"),
		gtxt("Degrees of Freedom"),
		gtxt("Pseudo R=Squared"),
		gtxt("AIC"),
		gtxt("BIC"),
		gtxt("Iterations"),
		gtxt("Fisher Scoring Iterations")
	)

	etype = switch(res$type, 
		"ML" = gtxt("Maximum Likelihood"),
		"BC" = gtxt("ML Bias Corrected"),
		"BR" = gtxt("ML Bias Reduced")
	)

	loglik = logLik(res2)
    vals = c(info["dep"],
		info["link"],
		info["plink"],
		ifelse(is.null(info[["offsetv"]]), gtxt("--NA--"), info["offsetv"]),
		info[["missingv"]],
		etype,
		ifelse(is.null(info[["dataset"]]), gtxt("--NA--"), info["dataset"]),
		ifelse(is.null(info[["rfile"]]), gtxt("--NA--"), info[["rfile"]]),
		ifelse(info[["retain"]], gtxt("Yes"), gtxt("No")),
		info[["modeldate"]],
		ifelse(res["converged"], gtxt("Yes"), gtxt("No")),
		info['ncases'],
		res["n"],
		round(loglik[[1]], 5),
		round(attr(loglik, "df"), 5),  #d.f.
		ifelse(is.na(res$pseudo.r.squared), gtxt("--NA--"), round(res$pseudo.r.squared,5)),
		round(AIC(res2), 5),
		round(BIC(res2),5),
		res["iterations"][[1]][[1]],
		res["iterations"][[1]][[2]]
    )
	
    # settings and result summary
    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="PROSUMMARY", outline=gtxt("Summary"),
		caption = gtxt("Computations done by R package betareg")
	)
     
	#coefficients table
	ctable = data.frame(coef(res)$mean)
	cnames = c(gtxt("Coefficient"), gtxt("Std. Error"), gtxt("Z Value"), gtxt("Sig."))
	names(ctable) = cnames
	spsspivottable.Display(ctable, title=gtxt("Mean Coefficients"),
		templateName="PROMEANCOEF", outline=gtxt("Mean Coefficients")
	)
	ctable = data.frame(coef(res)$precision)
	names(ctable) = cnames
	spsspivottable.Display(ctable, title=gtxt("Precision Coefficients"),
		templateName="PROPRECOEF", outline=gtxt("Precision Coefficients")
	)
	# plots

	plotreq = c(info['residso'],info['cooks'],info['leverage'],
		info['residslin'],info['halfnormal'],info['predobs'])
	plotseq = (1:6)[unlist(plotreq)]
	if (any(plotseq)) {
		plot(res2, which=plotseq, type=info[['presidtype']], sub.caption="", lwd=2)
	}

    warns$display(inproc=TRUE)
}

createdataset = function(res, res2, dta, details) {
    # Create residuals and/or predicted values dataset
    # Dataset name is known to be okay, and procedure state is ended
	# res is the summary results
	# res2 is the betareg object (non-summary)
	# dta is the data
	# iddata is the id variable (optional)
	# details is all the input arguments

    if (!is.null(details[["dataset"]])) {
		if (is.null(details[["fittedv"]]) && is.null(details[["resids"]])) {
			stop(gtxt("An output dataset was specified but no variables were requested"), call.=FALSE)
		}
	# construct a data frame with the requested variables
	# If the dataset has missing values, the predict and residuals output will have
	# No values for those, but they will carry the ID value, so we have to
	# pick up the ID values from the row names in the generated data frame and
	# copy them to the ID variable for sending back to Statistics
	# n.b. The doc for predict and residuals claims that these values are
	# present but NA, but the facts are otherwise.
	
		if (!is.null(details[["fittedv"]])) {
			ptype = list(responsemean="response",link="link",responsevar="variance")[[details[["fittedv"]]]]
			pred = predict(res2, type=ptype)
			theframe = data.frame(pred)
		}
		if (!is.null(details[["resids"]])) {
			resids = residuals(res2, details[["resids"]])
			if (is.null(details[["fittedv"]])) {
				theframe = data.frame(resids)
			} else {
				theframe[2] = resids
			}
		}
		theframe = data.frame(id=row.names(theframe), theframe)
        if (!is.null(details[["id"]])) {  # was an id variable provided
            vardict = spssdictionary.GetDictionaryFromSPSS(details["id"])
            loc = match(details[["id"]], vardict["varName",])
            idtype = as.integer(vardict[["varType", loc]])
            idformat = vardict[["varFormat", loc]]
            idlabel = vardict[["varLabel", loc]]
			if (idlabel == "") {  # if no label, use the id variable name as the label
				idlabel = details[["id"]]
			}
        } else {
            idtype = 0
            idformat = "F10.0"
            idlabel = ""
        }
		dictspec = list(c("ID", idlabel, idtype, idformat, "nominal"))
		if (!is.null(details["fittedv"])) {
			dictspec[2] = list(c("FittedValues", gtxtf("Fitted Values, Type=%s", details["fittedv"]),
				0, "F8.4", "scale"))
		}
		if (!is.null(details["resids"])) {
			dictspec[length(dictspec)+1] = list(c("Residuals", gtxtf("Residuals, Type=%s", details["resids"]),
				0, "F8.4", "scale"))
		}

        dict = spssdictionary.CreateSPSSDictionary(dictspec)
        spssdictionary.SetDictionaryToSPSS(details[["dataset"]], dict)
        spssdata.SetDataToSPSS(details[["dataset"]], theframe)
        spssdictionary.EndDataStep()
    }
}

dscheck = function(alldsspecs, warns) {
	# check dataset validation conditions
	
    if (!is.null(alldsspecs)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name if creating new datasets"),
                dostop=TRUE)
        }
        if (length(intersect(alldsspecs, alldatasets) > 0)) {
            warns$warn(gtxt("One or more specified output dataset names are already in use"),
                dostop=TRUE)
        }
    }
}

buildformula = function(dep, indep, indeppre, offsetv) {
	# return formula object
	# dep is the dependent variable name
	# indep is the list of independent variable names
	# indeppre is the list of variables in the precision model, if any
	# offset is the means offset or NULL
	
	main = paste(dep, paste(indep, collapse="+"), sep="~")
	if (!is.null(offsetv)) {
		main = paste(main, sprintf("offset(%s)", offsetv), sep="+")
	}
	if (!is.null(indeppre)) {
		main = paste(main, paste(indeppre, collapse="+"), sep="|")
	}
	return(as.formula(main))
}
	
# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_PROPOR_REGR"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_PROPOR_REGR"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

    if (is.null(msg) || dostop) {
        lcl$display(inproc)  # display messages and end procedure state
        if (dostop) {
            stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
        }
    }
}

    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

    if (lcl$msgnum == 0) {   # nothing to display
        if (inproc) {
            spsspkg.EndProcedure()
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        }
        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
            rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

    for (i in 1:lcl$msgnum) {
        rowcategory = spss.CellText.String(as.character(i))
        BasePivotTable.SetCategories(table,rowdim,rowcategory)
        BasePivotTable.SetCellValue(table,rowcategory, 
            spss.CellText.String(lcl$msglist[[i]]))
    }
    spsspkg.EndProcedure()   # implies display
} else { # can't produce a table
    for (i in 1:lcl$msgnum) {
        print(lcl$msglist[[i]])
    }
}
}
}
return(lcl)
}

Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=FALSE),
		spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        spsspkg.Template("INDEPPRE", subc="",  ktype="existingvarlist", var="indeppre", islist=TRUE),
		spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id", islist=FALSE),
		spsspkg.Template("LINK", subc="", ktype="str", var="link",
			vallist = list("logit","probit","cloglog","cauchit","log","loglog")),
		spsspkg.Template("PLINK", subc="", ktype="str", var="plink",
			vallist = list("log","sqrt","identity")),
        spsspkg.Template("OFFSET", subc="", ktype="existingvarlist", var="offsetv"),
		
		spsspkg.Template("ESTIMATOR", subc="OPTIONS", ktype="str", var="etype",
			vallist = list("ml", "bc", "br")),
        spsspkg.Template("MISSING", subc="OPTIONS",  ktype="str", 
            var="missingv", vallist=list("omit", "stop")),
        spsspkg.Template("MAXITER", subc="OPTIONS",  ktype="int", 
            var="maxiter", vallist=list(1)),
        spsspkg.Template("FSMAXITER", subc="OPTIONS",  ktype="int", 
            var="fsmaxiter", vallist=list(0)),
        spsspkg.Template("FSTOL", subc="OPTIONS",  ktype="float", 
            var="fstol", vallist=list(0)),
		spsspkg.Template("RESIDS", subc="OPTIONS", ktype="bool", var="residso"),
		spsspkg.Template("COOKS", subc="OPTIONS", ktype="bool", var="cooks"),
		spsspkg.Template("LEVERAGE", subc="OPTIONS", ktype="bool", var="leverage"),
		spsspkg.Template("RESIDSLIN", subc="OPTIONS", ktype="bool", var="residslin"),
		spsspkg.Template("HALFNORMAL", subc="OPTIONS", ktype="bool", var="halfnormal"),
		spsspkg.Template("PREDOBS", subc="OPTIONS", ktype="bool", var="predobs"),
		spsspkg.Template("PRESIDTYPE", subc="OPTIONS", ktype="str", var="presidtype",
			vallist=list("pearson", "sweighted2", "deviance", "response", "weighted", "sweighted")),
		
		spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
		spsspkg.Template("RESIDUALS", subc="SAVE", ktype="str", var="resids",
			vallist=list("pearson","sweighted2","deviance","response","weighted","sweighted")),
		spsspkg.Template("FITTED", subc="SAVE", ktype="str", var="fittedv",
			vallist=list("responsemean","link","responsevar")),
		spsspkg.Template("RFILE", subc="SAVE", ktype="literal", var="rfile"),
		spsspkg.Template("RETAIN", subc="SAVE", ktype="bool", var="retain")
    ))        
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    } else {
        res <- spsspkg.processcmd(oobj,args,"proext")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
