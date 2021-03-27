
#' make a SNP list
#'
#' Create a list of rsids corresponding to "top hits" in the GWAS catalog, the 1000 genomes super popualtions and SNPs of specific interest to the user (e.g. genetic instruments/proxies for the exposure of interest).
#'
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param ref1000G_superops include reference SNPs from 1000 genomes super populations. Default=TRUE
#' @param snplist_user character vector of user specified rsids. 
#'
#' @return character vector
#' @export
#' @examples
#' snplist<-make_snplist(efo_id="EFO_0006859") 

make_snplist<-function(trait=NULL,efo_id=NULL,efo=NULL,ref1000G_superops=TRUE,snplist_user=NULL){
	requireNamespace("gwasrapidd", quietly=TRUE)
	
	if(!is.null(efo_id)){
		top_hits<-gwas_catalog_hits(efo_id=efo_id)	
	}

	if(!is.null(trait)){
		top_hits<-gwas_catalog_hits(trait=trait)	
	}
	if(!is.null(efo)){
		top_hits<-gwas_catalog_hits(efo=efo)
	}
	
	top_hits_rsids<-top_hits$rsid	
	
	# # snplist<-c(snplist,top_hits_rsids)
	utils::data("refdat_1000G_superpops",envir =environment())
	snplist<-c(top_hits_rsids,unique(refdat_1000G_superpops$SNP))	
	if(!is.null(snplist_user)){
		snplist<-c(snplist,snplist_user)
	}
	# snplist<-unique(snplist)		
	
	return(unique(snplist))
}

#' GWAS top hits 
#'
#' Extract results for top hits for the trait of interest from the NHGRI-EBI GWAS catalog
#'
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of intersest in the experimental factor ontology
#'
#' @return data frame
#' @export

gwas_catalog_hits<-function(trait=NULL,efo=NULL,efo_id=NULL){
	requireNamespace("gwasrapidd", quietly=TRUE)

	if(!is.null(efo)){
		efo<-trimws(unlist(strsplit(efo,split=";")))	
		gwas_studies<-gwasrapidd::get_studies(efo_trait = efo)		
		# unique(gwas_studies@studies$reported_trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for efo -",efo,"- returned 0 studies from the GWAS catalog"))
		}
	}

	if(!is.null(efo_id)){
		efo_id<-trimws(unlist(strsplit(efo_id,split=";")))	
		gwas_studies<-gwasrapidd::get_studies(efo_id = efo_id)		
		# unique(gwas_studies@studies$reported_trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for efo -",efo_id,"- returned 0 studies from the GWAS catalog"))
		}
	}
	
	if(!is.null(trait)){
		gwas_studies<-gwasrapidd::get_studies(reported_trait = trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for trait -",trait,"- returned 0 studies from the GWAS catalog"))
		}
	}
	
	if(nrow(gwas_studies@studies)!=0){
		ancestry_tab<-make_ancestry_table(gwas_studies=gwas_studies)		
		study_ids<-gwas_studies@studies$study_id	
		Dat<-NULL	
		for(i in 1:length(study_ids)){		
		# print(i)			
			gwas_associations<-gwasrapidd::get_associations(study_id = study_ids[i])
			if(nrow(gwas_associations@associations)!=0){
				associations<-data.frame(gwas_associations@associations,stringsAsFactors=F)
				risk_alleles<-data.frame(gwas_associations@risk_alleles)
				gwas_results<-merge(associations,risk_alleles,by="association_id")
				# p_values<-gwas_associations@associations$pvalue			
				# odds_ratios<-gwas_associations@associations$or_per_copy_number
				# risk_alleles<-gwas_associations@risk_alleles$risk_allele
				# rsids<-gwas_associations@risk_alleles$variant_id
				# eaf<-gwas_associations@risk_alleles$risk_frequency
			
				gwas_results$z_scores<-stats::qnorm(gwas_results$pvalue/2,lower.tail=F)
				gwas_results$log_odds_ratios<-log(gwas_results$or_per_copy_number)
				all(is.na(gwas_results$log_odds_ratios))
				gwas_results$log_odds_ratios<-gwas_results$beta_number
				gwas_results$standard_errors<-gwas_results$log_odds_ratios/gwas_results$z_scores
				gwas_results$study_id<-study_ids[i]								
				Dat[[i]]<-gwas_results	
			}			
		}
		if(!is.null(trait)) trait_efo<-trait
		if(!is.null(efo)) trait_efo<-efo

		if(is.null(Dat)){
			warning(paste0("no results found in GWAS catalog for ",trait_efo))
		}
		if(!is.null(Dat)){
			Dat<-do.call(rbind,Dat)
			Dat<-Dat[,c("variant_id","risk_allele","log_odds_ratios","standard_errors","risk_frequency","pvalue","z_scores","study_id")]
			names(Dat)<-c("rsid","Effect.Allele","lnor","se","eaf","p","test_statistic","study_id")		
			Dat<-merge(Dat,ancestry_tab,by="study_id")				
			return(Dat)
		}
	}			
}


make_ancestry_table<-function(gwas_studies=NULL){
	ancestry_tab<-data.frame(gwas_studies@ancestral_groups,stringsAsFactors=F)		
	ancestry_tab<-unique(ancestry_tab[,c("study_id","ancestral_group")])
	# ancestry_tab<-ancestry_tab[!is.na(ancestry_tab$ancestral_group),]
	study_ids<-unique(ancestry_tab$study_id)
	anc2<-NULL
	for(i in 1:length(study_ids)){
		anc1<-ancestry_tab[ancestry_tab$study_id == study_ids[i],]
		anc1$ancestral_group<-paste(anc1$ancestral_group,collapse="; ")		
		anc2[[i]]<-anc1[1,]
	}
	anc2<-do.call(rbind,anc2)
	return(anc2)
}

