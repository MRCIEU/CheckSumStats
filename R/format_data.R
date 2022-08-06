#' format data
#'
#' Get the trait summary data ready for the QC checks. 
#'
#' @param dat the dataset to be formatted
#' @param trait the name of the trait.  
#' @param population describe the population ancestry of the dataset
#' @param ncase number of cases or name of the column specifying the number of cases
#' @param ncontrol number of controls or name of the column specifying the number of controls. If your summary data was generated in a linear model of a continuous trait, use ncontrol to indicate the total sample size. 
#' @param rsid name of the column containing the rs number or identifiers for the genetic variants
#' @param effect_allele name of the effect allele column 
#' @param other_allele name of the non-effect allele column 
#' @param beta name of the column containing the SNP effect sizes. Use this argument if your summary data was generated in a linear model of a continuous trait. 
#' @param se standard error for the beta. Use this argument if your summary data was generated in a linear model of a continuous trait.
#' @param lnor name of the column containing the log odds ratio. If missing, tries to infer it from the odds ratio
#' @param lnor_se name of the column containing the standard error for the log odds ratio. If missing, tries to infer it from 95% confidence intervals or pvalues
#' @param or name of column containing the odds ratio 
#' @param or_lci name of column containing the lower 95% confidence interval for the odds ratio 
#' @param or_uci name of column containing the upper 95% confidence interval for the odds ratio 
#' @param eaf name of the effect allele frequency column
#' @param p name of the pvalue columne
#' @param chr name of the column containing the chromosome number for each genetic variant
#' @param pos genomic position for the genetic variant in base pairs 
#' @param z_score effect size estimate divided by its standard error
#' @param drop_duplicate_rsids drop duplicate rsids? logical. default TRUE. duplicate rsids may for example correspond to triallelic SNPs. 
#'
#' @return data frame
#' @export

format_data<-function(dat=NULL,trait=NA,population=NA,ncase=NA,ncontrol=NA,rsid=NA,effect_allele=NA,other_allele=NA,beta=NA,se=NA,lnor=NA,lnor_se=NA,eaf=NA,p=NA,or=NA,or_lci=NA,or_uci=NA,chr=NA,pos=NA,z_score=NA,drop_duplicate_rsids=TRUE){
	# summary_set="FAsnps"

	# if(any(is.na(dat[,rsid]) |  dat[,rsid] ==".")) stop("rsid missing")
	if(rsid=="ID") 
	{ 
		names(dat)[names(dat) == rsid]<-"rsid"
		rsid<-"rsid"
	}

	# sometimes odds ratio and confidence intervals are reported
	if(!is.na(or) & is.na(lnor_se) & !is.na(or_uci))
	{
		dat$lnor<-log(as.numeric(dat[,or]))
		dat$lnor_se<-(log(as.numeric(dat[,or_uci]))-log(as.numeric(dat[,or_lci])))/(1.96*2)
		lnor<-"lnor"
	}

	# odds ratio and p value but no standard error or confidence intervals
	if(!is.na(or) & is.na(lnor_se) & is.na(or_uci))
	{
		dat$lnor<-log(as.numeric(dat[,or]))
		dat$z<-stats::qnorm(dat[,p]/2,lower.tail=F)
		dat$lnor_se<-abs(dat$lnor)/dat$z
	}

	if(!is.na(lnor) & is.na(lnor_se) & is.na(or_uci))
	{
		dat$z<-stats::qnorm(dat[,p]/2,lower.tail=F)
		dat$lnor_se<-abs(dat[,lnor])/dat$z
	}

	# sometimes Odds ratio and standard error of log odds ratio are reported
	if(!is.na(or) & !is.na(lnor_se))
	{
		dat$lnor<-log(as.numeric(dat[,or]))
		lnor<-"lnor"
	}

	if(is.na(p)  & is.na(z_score))
	{
		dat$z_score<-abs(as.numeric(dat[,lnor])/as.numeric(dat[,lnor_se]))
		dat$p<-stats::pnorm(dat$z_score ,lower.tail=F)*2
	}	

	if(is.na(p) & !is.na(z_score))
	{
		dat$p<-stats::pnorm(dat[,z_score] ,lower.tail=F)*2
    }


	Name_cols<-c("rsid","effect_allele","other_allele","beta","se","lnor","lnor_se","eaf","p","chr","pos","z_score")


	if(is.numeric(ncase) | is.numeric(ncontrol))
	{
		dat$ncase<-ncase
		dat$ncontrol<-ncontrol
	}else{
		Name_cols<-c(Name_cols,"ncase","ncontrol")
	}
	dat$trait<-trait
	dat$population<-population
	
	for(i in 1:length(Name_cols))
	{
		# print(i)
		if(!is.na(eval(parse(text=Name_cols[i]))))
		{
			names(dat)[names(dat) == eval(parse(text=Name_cols[i]))]<-Name_cols[i]
		}
	}


	dat$p<-as.numeric(dat$p)
	
	if("lnor" %in% names(dat))
	{
		dat$lnor<-as.numeric(dat$lnor)
		dat$lnor_se<-as.numeric(dat$lnor_se)
		study_id_temp<-paste(dat$rsid,dat$effect_allele,dat$other_allele,dat$lnor,dat$lnor_se)

		dat<-dat[!is.na(dat$lnor) & !is.na(dat$lnor_se),]
		if(any(is.na(dat$lnor_se)) & is.na(p))
		{ #is this redundant?
			dat<-dat[!is.na(dat$lnor_se), ]	
		}

		if(any(is.na(dat$lnor_se)) & !is.na(p))
		{		
			Dat1<-dat[is.na(dat$lnor_se),]
			if(sum(Dat1$p)!=0) stop("infer missing SE from p value")
		}
		dat<-dat[!is.na(dat$lnor_se),] #is this redundant?
		dat<-dat[which(dat$lnor_se != "Inf"),]
		dat<-dat[which(dat$lnor_se != 0), ]
	}
	
	if("beta" %in% names(dat))
	{
		dat$beta<-as.numeric(dat$beta)
	}

	if("se" %in% names(dat))
	{
		dat$se<-as.numeric(dat$se)
	}

	if(!is.na(eaf))
	{
		dat$eaf<-as.numeric(dat$eaf)
	}

	dat$effect_allele<-toupper(dat$effect_allele)
	if(!is.na(other_allele))
	{		
		dat$other_allele<-toupper(dat$other_allele)
	}
	
	# this script may not work when a study has a beta column that is set to NA and also a separate column called lnor that takes on numeric values. 
	if("beta" %in% names(dat))
	{
		study_id_temp<-paste(dat$rsid,dat$effect_allele,dat$other_allele,dat$beta,dat$se)
	}

	# remove duplicate rows. 
	# There are two types of duplicates. Those where the rsid and the results for the rsid are duplicated and those where only the rsid is duplicated (i.e. results vary across duplicate rsids). We first deal with the duplicates where results are also duplicated, retaining one of the duplicate rsids. Then we deal with the duplicates where only the rsid is duplicated. for the latter we drop the rsid entirely, i.e. we don't retain one of the duplicates. 

	# sometimes only signed Z scores are provided. Only drop the duplicates if the dataset contains lnor_se or se. This script may not work if a dataset has these columns present and they are set to NA. The objective is to drop duplicate rows (duplicated on rsid and the result/SNP-trait association). 
	if(any(c("lnor_se","se") %in% names(dat))) 
	{	
		dat<-dat[!duplicated(study_id_temp),]
	}
	
	# Drop duplicates rsids. These seem to usually correspond to trialleic SNPs
	if(drop_duplicate_rsids){
		Dups<-unique(dat$rsid[duplicated(dat$rsid)])	
		dat<-dat[!dat$rsid %in% Dups,]
	}
	

	# dat$summary_set<-summary_set
	# dat$all_summary_stats<-FALSE
	# if(all_summary_stats){
	# 	dat$all_summary_stats<-TRUE
	# 	dat$summary_set <- "full summary stats"
	# }

	return(dat)
}	

#' get efo
#'
#' Retrieve the experimental factor ontology (EFO) for some trait of interest. EFOs are retrieved from ZOOMA https://www.ebi.ac.uk/spot/zooma/
#'
#' @param trait the trait of interest
#'
#' @return list 
#' @export

get_efo<-function(trait=NULL){
	trait<-gsub(" ","+",trait)
	# system(paste0("wget \"https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=",trait,"\"  -O \"zooma2.txt\""))
	
	# curl.cmd<-paste0("https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=",trait) #this no longer works
	
	curl.cmd<-paste0("https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=",trait,"&filter=required:[gwas]")
	dat<-readLines(curl::curl(curl.cmd),warn=FALSE)
	dat<-gsub("\"","",dat)
	dat<-strsplit(dat,split=",")
	dat<-unlist(dat)
	EFO<-dat[grep("semanticTag",dat)]
	
	dat[grep("propertyValue",dat)]

	EFO<-EFO[grep("EFO",EFO)]
	Start<-unlist(gregexpr("EFO",EFO))
 	End<-unlist(lapply(gregexpr("([0-9])",EFO),FUN=function(x) 
		unlist(x)[length(unlist(x))]))
	EFO<-unique(substring(EFO,Start,End))
	confidence<-unique(dat[grep("confidence",dat)])
	if(length(EFO)==0) return("no EFO found")
	return(list("efo_id"=EFO,"confidence"=confidence))
}
