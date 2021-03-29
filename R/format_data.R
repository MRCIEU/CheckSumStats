#' format data
#'
#' Get the outcome summary data ready for the QC checks. The function assumes that the summary data was generated in a logistic regression model of a binary phenotype, e.g. cases and controls
#'
#' @param dat the the dataset to be formatted
#' @param outcome the name of the outcome trait
#' @param population describe the ancestry of the dataset
#' @param pmid pubmed identifier for a publication corresponding to the dataset
#' @param study name of the study
#' @param ncase number of cases or name of the column specifying the number of cases
#' @param ncontrol number of controls or name of the column specifying the number of controls
#' @param UKbiobank was the dataset generated in UK Biobank or in a dataset that overlaps with UK Biobank? default set to FALSE
#' @param rsid name of the column containing the rs number or identifiers for the genetic variants
#' @param effect_allele name of the effect allele column 
#' @param other_allele name of the non-effect allele column 
#' @param lnor name of the column containing the log odds ratio. If missing, tries to infer it from the odds ratio
#' @param se name of the column containing the standard error for the log odds ratio. If missing, tries to infer it from 95% confidence intervals or pvalues
#' @param or name of column containing the odds ratio 
#' @param lci name of column containing the lower 95% confidence interval for the odds ratio 
#' @param uci name of column containing the upper 95% confidence interval for the odds ratio 
#' @param eaf name of the effect allele frequency column
#' @param p name of the pvalue columne
#' @param info1 name of the column containing metrics of imputation quality, such as info or r2 scores 
#' @param info2 name of second column containing metrics of imputation quality, such as info or r2 scores 
#' @param info3 name of third column containing metrics of imputation quality, such as info or r2 scores name of the third info score column
#' @param info4 name of fourth column containing metrics of imputation quality, such as info or r2 scores 
#' @param HWEp name of the column containing pvalues for tests of Hardy-Weinberg equilibrium

#' @param phet name of the column containing pvalues for tests of between study heterogeneity

#' @param I2 name of the column containing I2 metric of between study heterogeneity
#' @param Q name of the column containing Cochran Q test for between study heterogeneity
#' @param Direction name of the column specifying the direction of effect in each individual study within a meta-analysis of multiple studies
#' @param effect_allele_confirmed have you confirmed the identity of the effect allele column? default set to FALSE
#' @param chr name of the column containing the chromosome number for each genetic variant
#' @param pos genomic position for the genetic variant in base pairs 
#' @param ID identifier for the dataset 
#' @param test_statistic log odds ratio divided by standard error of log odds ratio
#' @param all_summary_stats do you have the full set of summary data or only a subset of the summary data? Default set to FALSE 
#' @param open_gwas were the summary data downloaded from Open GWAS https://gwas.mrcieu.ac.uk/ ? Default set to FALSE
#' @param efo outcome trait of interest in the experimental factor ontology.  
#' @param efo_id ID for outcome trait of interest in the experimental factor ontology.  


#'
#' @return data frame
#' @export

format_data<-function(dat=NULL,outcome=NA,population=NA,pmid=NA,study=NA,ncase=NA,ncontrol=NA,UKbiobank=NA,rsid=NA,effect_allele=NA,other_allele=NA,lnor=NA,se=NA,eaf=NA,p=NA,info1=NA,info2=NA,info3=NA,info4=NA,HWEp=NA,phet=NA,I2=NA,Q=NA,Direction=NA,effect_allele_confirmed=FALSE,or=NA,lci=NA,uci=NA,chr=NA,pos=NA,ID=NULL,test_statistic=NA,all_summary_stats=FALSE,open_gwas=FALSE,efo=NA,efo_id=NA){
	# summary_set="FAsnps"

	# if(any(is.na(dat[,rsid]) |  dat[,rsid] ==".")) stop("rsid missing")
	if(rsid=="ID") { 
		names(dat)[names(dat) == rsid]<-"rsid"
		rsid<-"rsid"
	}
	# if(!is.null(ref)){
	# 	ref_dat<-read.table(ref,sep=" ",head=F,stringsAsFactors=F)		
	# 	Chr<-unlist(strsplit(ref_dat$V1,split="chr"))
	# 	Chr<-Chr[Chr!=""]
	# 	ref_dat$chr<-Chr
	# 	if(is.na(chr)){
	# 		Chr<-unlist(strsplit(dat[,rsid],split=":"))
	# 		dat$Chr<-as.numeric(Chr[seq(1,by=2,length(Chr))])
	# 		dat$Pos<-as.numeric(Chr[seq(2,by=2,length(Chr))])
	# 		chr<-"Chr"
	# 		pos<-"Pos"
	# 	}

	# 	dat<-merge(dat,ref_dat,by.x=c(chr,pos),by.y=c("chr","V2"),all.x=T)

	# 	if(!is.na(rsid)){
	# 		dat$V4[is.na(dat$V4)]<-dat[is.na(dat$V4),rsid] 
	# 	}

	# 	dat<-dat[!is.na(dat$V4),]
	# 	rsid<-"V4"		
	# 	# dat<-dat[grep("rs",dat$V4),]				
	# }

	if(!is.null(ID)){
		dat$ID <- ID
	}

	# sometimes odds ratio and confidence intervals are reported
	if(!is.na(or) & is.na(se) & !is.na(uci)){
		dat$lnor<-log(dat[,or])
		dat$se<-(log(dat[,uci])-log(dat[,lci]))/(1.96*2)
		lnor<-"lnor"
	}

	# odds ratio and p value but no standard error or confidence intervals
	if(!is.na(or) & is.na(se) & is.na(uci)){
		dat$lnor<-log(dat[,or])
		dat$z<-stats::qnorm(dat[,p]/2,lower.tail=F)
		dat$se<-abs(dat$lnor)/dat$z
	}

	if(!is.na(lnor) & is.na(se) & is.na(uci)){
		dat$z<-stats::qnorm(dat[,p]/2,lower.tail=F)
		dat$se<-abs(dat[,lnor])/dat$z
	}

	# sometimes Odds ratio and standard error of log odds ratio are reported
	if(!is.na(or) & !is.na(se)){
		dat$lnor<-log(as.numeric(dat[,or]))
		lnor<-"lnor"
	}

	if(is.na(p)  & is.na(test_statistic)){
		dat$test_statistic<-abs(as.numeric(dat[,lnor])/as.numeric(dat[,se]))
		dat$p<-stats::pnorm(dat$test_statistic ,lower.tail=F)*2
	}	

	if(is.na(p) & !is.na(test_statistic)){
		dat$p<-stats::pnorm(dat[,test_statistic] ,lower.tail=F)*2
    }



	Name_cols<-c("rsid","effect_allele","other_allele","lnor","se","eaf","p","info1","info2","info3","info4","HWEp","phet","I2","Q","Direction","chr","pos","test_statistic")


	if(is.numeric(ncase)){
		dat$ncase<-ncase
		dat$ncontrol<-ncontrol
	}else{
		Name_cols<-c(Name_cols,"ncase","ncontrol")
	}
	dat$pmid<-pmid
	dat$outcome<-outcome
	dat$population<-population
	dat$study<-study
	dat$UKbiobank<-UKbiobank
	dat$effect_allele_confirmed<-effect_allele_confirmed

	for(i in 1:length(Name_cols)){
		# print(i)
		if(!is.na(eval(parse(text=Name_cols[i])))){
			names(dat)[names(dat) == eval(parse(text=Name_cols[i]))]<-Name_cols[i]
		}
	}
	

	dat$p<-as.numeric(dat$p)
	dat$se<-as.numeric(dat$se)
	dat$lnor<-as.numeric(dat$lnor)
	
	if(!is.na(eaf)){
		dat$eaf<-as.numeric(dat$eaf)
	}

	dat$effect_allele<-toupper(dat$effect_allele)
	if(!is.na(other_allele)){		
		dat$other_allele<-toupper(dat$other_allele)
	}
	# dat<-dat[!duplicated(dat$rsid),]
	# drop duplicated SNPs

	# remove duplicate rows. 
	# There are two types of duplicates. Those where the rsid and the results for the rsid are duplicated and those where only the rsid is duplicated (i.e. results vary across duplicate rsids). We first deal with the duplicates where results are also duplicated, retaining one of the duplicate rsids. Then we deal with the duplicates where only the rsid is duplicated. for the latter we drop the rsid entirely, i.e. we don't retain one of the duplicates. 
	study_id_temp<-paste(dat$rsid,dat$effect_allele,dat$other_allele,dat$lnor,dat$se)
	dat<-dat[!duplicated(study_id_temp),]

	# Dups<-unique(study_id_temp[duplicated(study_id_temp)])
	# Pos.dups<-which(study_id_temp %in% Dups)
	# Temp<-dat[Pos.dups,]
	# head(Temp[order(Temp$rsid),])
	
	# Drop duplicates rsids. These seem to usually correspond to trialleic SNPS
	Dups<-unique(dat$rsid[duplicated(dat$rsid)])	
	dat<-dat[!dat$rsid %in% Dups,]
	
	dat<-dat[!is.na(dat$lnor) & !is.na(dat$se),]
	if(any(is.na(dat$se)) & is.na(p)){
		dat<-dat[!is.na(dat$se), ]	
	}

	if(any(is.na(dat$se)) & !is.na(p)){		
		Dat1<-dat[is.na(dat$se),]
		if(sum(Dat1$p)!=0) stop("infer missing SE from p value")
	}
	dat<-dat[!is.na(dat$se),]
	dat<-dat[which(dat$se != "Inf"),]
	dat<-dat[which(dat$se != 0), ]

	# dat$summary_set<-summary_set
	# dat$all_summary_stats<-FALSE
	# if(all_summary_stats){
	# 	dat$all_summary_stats<-TRUE
	# 	dat$summary_set <- "full summary stats"
	# }

	dat$open_gwas<-FALSE
	if(open_gwas){
		dat$open_gwas<-TRUE
	}
	dat$efo<-paste(efo,collapse="; ")
	dat$efo_id<-paste(efo_id,collapse="; ")
	return(dat)
}	


