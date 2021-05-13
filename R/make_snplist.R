
#' make a SNP list
#'
#' Create a list of rsids corresponding to "top hits" in the GWAS catalog, the 1000 genomes super popualtions and SNPs of specific interest to the user (e.g. genetic instruments/proxies for the exposure of interest).
#'
#' @param trait the name of the trait in the NHGRI-EBI GWAS catalog
#' @param efo_id experimental factor ontology ID for trait of interest 
#' @param efo experimental factor ontology for the trait of interest
#' @param ref1000G_superpops include reference SNPs from 1000 genomes super populations. Default=TRUE
#' @param snplist_user character vector of user specified rsids. 
#'
#' @return character vector
#' @export
#' @examples
#' snplist<-make_snplist(efo_id="EFO_0006859") 

make_snplist<-function(trait=NULL,efo_id=NULL,efo=NULL,ref1000G_superpops=TRUE,snplist_user=NULL){
	requireNamespace("gwasrapidd", quietly=TRUE)
	
	snplist<-""
	if(!is.null(efo_id)){
		top_hits<-gwas_catalog_hits(efo_id=efo_id)	
		snplist<-top_hits$rsid	
	}

	if(!is.null(trait)){
		top_hits<-gwas_catalog_hits(trait=trait)	
		if(!is.null(efo_id)){
			snplist1<-top_hits$rsid	
			snplist<-c(snplist,snplist1)
		}else{
			snplist<-top_hits$rsid	
		}
	}

	if(!is.null(efo)){
		top_hits<-gwas_catalog_hits(efo=efo)
		if(!is.null(efo_id) | !is.null(trait)){
			snplist1<-top_hits$rsid	
			snplist<-c(snplist,snplist1)
		}else{		
			snplist<-top_hits$rsid	
		}
		
	}
		
	# # snplist<-c(snplist,top_hits_rsids)
	if(ref1000G_superpops){
		utils::data("refdat_1000G_superpops",envir =environment())
		snplist<-unique(c(snplist,unique(refdat_1000G_superpops$SNP)))
	}
	
	if(!is.null(snplist_user)){
		snplist<-c(snplist,snplist_user)
	}
	# snplist<-unique(snplist)		
	snplist<-snplist[snplist!=""]	
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
		# head(gwas_studies@studies)
		# which(gwas_studies@studies$study_id=="GCST001633")
		 # gwasrapidd::get_studies(study_id = "GCST001633")
		# unique(gwas_studies@studies$reported_trait)
		if(class(unlist(gwas_studies)) == "character"){
			# if(nrow(gwas_studies@studies)==0){
			if(nrow(gwas_studies)==0){
				warning(paste("search for efo -",efo,"- returned 0 studies from the GWAS catalog"))
			}
		}
	}

	if(!is.null(efo_id)){
		efo_id<-trimws(unlist(strsplit(efo_id,split=";")))	
		gwas_studies<-gwasrapidd::get_studies(efo_id = efo_id)		
		# unique(gwas_studies@studies$reported_trait)
		if(class(unlist(gwas_studies)) == "character"){
			# if(nrow(gwas_studies@studies)==0){
			if(nrow(gwas_studies)==0){
				warning(paste("search for efo -",efo_id,"- returned 0 studies from the GWAS catalog"))
			}
		}
	}
	
	if(!is.null(trait)){
		gwas_studies<-gwasrapidd::get_studies(reported_trait = trait)
		if(class(unlist(gwas_studies)) == "character"){
			# if(nrow(gwas_studies@studies)==0){
			if(nrow(gwas_studies)==0){
				warning(paste("search for trait -",trait,"- returned 0 studies from the GWAS catalog"))
			}
		}
	}
	
	if(class(unlist(gwas_studies)) != "character"){
		if(nrow(gwas_studies@studies)!=0){
			ancestry_tab<-make_ancestry_table(gwas_studies=gwas_studies)	
			study_ids<-gwas_studies@studies$study_id	
			Dat<-NULL	
			# i<-4
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
					gwas_results$beta_gc<-log(gwas_results$or_per_copy_number)

					Test<-FALSE 
					if(all(is.na(gwas_results$beta_gc))){
						if(!all(is.na(gwas_results$beta_number))){
							gwas_results$beta_gc<-gwas_results$beta_number		
							Test<-any(is.na(gwas_results$beta_gc))
							if(Test){
								gwas_results2<-gwas_results[is.na(gwas_results$beta_gc),]
								gwas_results2$se_gc<-NA
							}
							gwas_results<-gwas_results[!is.na(gwas_results$beta_gc),]
							Pos1<-which(gwas_results$beta_direction == "increase")
							Pos2<-which(gwas_results$beta_direction == "decrease")
							if(!all(gwas_results$beta_gc>0)) stop("direction of beta_gc not always positive")
							gwas_results$beta_gc[Pos2]<-gwas_results$beta_gc[Pos2]*-1
							if(!all(unique(gwas_results$beta_direction) %in% c("increase","decrease"))) stop("beta_direction in gwas catalog not always increase or decrease")
						}
					}
					Pos<-which(!is.na(gwas_results$standard_error))			
					gwas_results$se_gc<-NA
					gwas_results$se_gc[Pos]<-gwas_results$standard_error[Pos]
					Pos<-which(is.na(gwas_results$se_gc))
					gwas_results$se_gc[Pos]<-abs(gwas_results$beta_gc[Pos]/gwas_results$z_scores[Pos])
					if(Test){
						gwas_results<-rbind(gwas_results,gwas_results2) 
					}

					# beta_gc still missing. This happens when there is at least one odds ratio amongst the rows. Sometimes the odds ratio colmn is missing beta_number corresponds to signed Z scores. 
					gwas_results$z.x<-gwas_results$beta_gc / gwas_results$se_gc 
					gwas_results1<-gwas_results[!is.na(gwas_results$beta_gc),]
					gwas_results2<-gwas_results[is.na(gwas_results$beta_gc),]
					# update to make allowance for presence of any z scores amongst the beta_units? 
					if(all(!is.na(gwas_results2$beta_unit) & gwas_results2$beta_unit == "z score")){
						gwas_results2$z.x<-gwas_results2$beta_number
						Pos<-which(gwas_results2$beta_direction=="decrease")
						gwas_results2$z.x[Pos]<-gwas_results2$z.x[Pos]*-1
						gwas_results<-rbind(gwas_results1,gwas_results2)
					}
					
					gwas_results$study_id<-study_ids[i]						
					Dat[[i]]<-gwas_results	
				}			
			}

			if(!is.null(efo_id)) trait_efo<-efo_id
			if(!is.null(efo)) trait_efo<-efo
			if(!is.null(trait)) trait_efo<-trait
			
			if(is.null(Dat)){
				warning(paste0("no results found in GWAS catalog for ",trait_efo))
			}
			if(!is.null(Dat)){
				Dat<-do.call(rbind,Dat)
				Dat<-Dat[,c("variant_id","risk_allele","beta_gc","se_gc","risk_frequency","pvalue","z_scores","study_id","z.x")]
				names(Dat)<-c("rsid","effect_allele","beta_gc","se_gc","eaf","p","test_statistic","study_id","z.x")		
				Dat<-merge(Dat,ancestry_tab,by="study_id")				
				return(Dat)
			}
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

