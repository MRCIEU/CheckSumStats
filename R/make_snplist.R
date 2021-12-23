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
#' snplist<-make_snplist(efo_id="EFO_0006859",ref1000G_superpops=FALSE) 
# CheckSumStats::make_snplist(efo_id="EFO_0006859") 
make_snplist<-function(trait=NULL,efo_id=NULL,efo=NULL,ref1000G_superpops=TRUE,snplist_user=NULL){
	requireNamespace("gwasrapidd", quietly=TRUE)
	
	snplist<-""
	 
	top_hits<-gwas_catalog_hits2(efo_id=efo_id,trait=trait,efo=efo)	
	
	if(class(top_hits)!="data.frame")
	{
		# if(nrow(top_hits))
		return(top_hits)
	}
			
	snplist<-top_hits$rsid

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
#' @param trait the trait of interest as reported in the GWAS catalog
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of intersest in the experimental factor ontology
#' @param map_association_to_study map associations to study in GWAS catalog. This supports matching of results on PMID and study ancestry, which increases accuracy of comparisons, but is slow when there are large numbers of associations. Default = TRUE
#'
#' @return data frame
#' @importFrom magrittr %>%
#' @export

gwas_catalog_hits2<-function(trait=NULL,efo=NULL,efo_id=NULL,map_association_to_study=TRUE)
{

	gwas_associations<-get_gwas_associations(reported_trait=trait,efo_trait=efo,efo_id=efo_id)
			

	if(class(gwas_associations) =="associations") 
	{
		
		if(nrow(gwas_associations@associations)!=0)	
		{
			
			associations<-data.frame(gwas_associations@associations,stringsAsFactors=F)
			risk_alleles<-data.frame(gwas_associations@risk_alleles)
			gwas_results<-merge(associations,risk_alleles,by="association_id")
			col.keep<-c("variant_id","association_id","risk_allele","beta_gc","se_gc","risk_frequency","pvalue","z_scores","z.x")
			names_col.keep<-c("rsid","association_id","effect_allele","beta_gc","se_gc","eaf","p","test_statistic","z.x")
		
						
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
			if(nrow(gwas_results2)>0){
				if(all(!is.na(gwas_results2$beta_unit) & gwas_results2$beta_unit == "z score")){
					gwas_results2$z.x<-gwas_results2$beta_number
					Pos<-which(gwas_results2$beta_direction=="decrease")
					gwas_results2$z.x[Pos]<-gwas_results2$z.x[Pos]*-1
					gwas_results<-rbind(gwas_results1,gwas_results2)
				}
			}
			
			if(!is.null(efo_id)) trait_efo<-efo_id
			if(!is.null(efo)) trait_efo<-efo
			if(!is.null(trait)) trait_efo<-trait
			
			if(is.null(gwas_results)){
				warning(paste0("no results found in GWAS catalog for ",trait_efo))
				return("no results found")
			}
			
			if(!is.null(gwas_results))			
			{
				if(map_association_to_study)
				{
					assoc2study<-map_association_to_study_id(associations=gwas_associations)		
					gwas_results<-merge(gwas_results,assoc2study,by="association_id")					
					ancestry_tab<-make_ancestry_table(association_id=gwas_results$association_id)					
					gwas_results<-merge(gwas_results,ancestry_tab,by="study_id")
					col.keep<-c(col.keep,"study_id")
					names_col.keep<-c(names_col.keep,"study_id")
				}					
				gwas_results<-gwas_results[,col.keep]
				names(gwas_results)<-c(names_col.keep)
				return(gwas_results)
			}
			
		}		
	}
	return(gwas_associations)	
}


map_association_to_study_id<-function(associations=NULL){
	association_ids <- associations@associations$association_id
	    names(association_ids) <- association_ids	
	
	studies <-purrr::map(association_ids, ~ gwasrapidd::get_studies(association_id = .x))
	
	
	association2study <-
	      purrr::imap_dfr(
	        studies,
	        ~ tibble::tibble(
	          association_id = .y,
	          study_id = .x@studies$study_id
	        )
	      )

	return(association2study)
}

make_ancestry_table<-function(association_id=NULL){
				
	studies<-gwasrapidd::get_studies(association_id = association_id)
	
	ancestries <-
      dplyr::left_join(studies@ancestries,
                       studies@ancestral_groups,
                       by = c('study_id', 'ancestry_id')) %>%
      dplyr::left_join(studies@countries_of_origin, by = c('study_id', 'ancestry_id')) %>%
      dplyr::rename(
        "co_country_name"="country_name",
        "co_major_area"="major_area",
        "co_region"="region"
      ) %>%
      dplyr::left_join(studies@countries_of_recruitment,
                       by = c('study_id', 'ancestry_id')) %>%      
      dplyr::rename(
        "cr_country_name"="country_name",
        "cr_major_area" ="major_area",
        "cr_region"="region"
      )
	ancestry_tab<-data.frame(ancestries,stringsAsFactors=F)
	ancestry_tab$ancestral_group[ancestry_tab$ancestral_group=="NA"]<-NA
	# if ancestral group is missing replace with major geographic area of recruitment
	ancestry_tab$ancestral_group[is.na(ancestry_tab$ancestral_group)]<-	ancestry_tab$cr_major_area[is.na(ancestry_tab$ancestral_group)]
	# ancestry_tab<-unique(ancestry_tab[,c("study_id","ancestral_group")])
	# ancestry_tab<-ancestry_tab[!is.na(ancestry_tab$ancestral_group),]
	study_ids<-unique(ancestry_tab$study_id)
	
	# for when ancestry varies within study, 
	anc2<-NULL
	for(i in 1:length(study_ids)){
		anc1<-ancestry_tab[ancestry_tab$study_id == study_ids[i],]
		Names<-names(anc1)
		for(j in 1:length(Names)){
			anc1[,Names[j]]<-paste(unique(anc1[,Names[j]]),collapse="; ")	
		}		
		if(nrow(unique(anc1))!=1) stop("unexpected number of rows in ancestry table")
		anc2[[i]]<-unique(anc1)
	}
	anc2<-do.call(rbind,anc2)
	return(anc2)
}


get_gwas_associations<-function(reported_trait=NULL,efo_trait=NULL,efo_id=NULL,verbose = FALSE,warnings = TRUE) 
{
  
	if(!is.null(reported_trait)) 
	{

		
		gwas_studies <- gwasrapidd::get_studies(reported_trait = reported_trait)
		reported_trait<-NULL 
		if(class(unlist(gwas_studies)) != "character")		
		{
			if(nrow(gwas_studies@studies)!=0)
			{	
				gwas_associations_by_reported_trait <- gwasrapidd::get_associations(study_id = gwas_studies@studies$study_id,verbose = verbose,warnings = warnings)
				reported_trait<-"results found"
				# sometimes a study is in the GWAS catalog but has no associations. 
				if(nrow(gwas_associations_by_reported_trait@associations)==0)
				{
					reported_trait<-NULL
				}
			}
		}
	}

	if(!is.null(efo_trait)) 
	{
		gwas_associations_by_efo_trait <- gwasrapidd::get_associations(efo_trait = efo_trait,verbose = verbose,warnings = warnings)
		# association objects are retrieved even if there are no associations in the GWAS catalog 
		if(nrow(gwas_associations_by_efo_trait@associations)==0)
		{
			efo_trait<-NULL
		}
	}

	if(!is.null(efo_id)) 
	{
		gwas_associations_by_efo_id <- gwasrapidd::get_associations(efo_id = efo_id,verbose = verbose,warnings = warnings)
		# sometimes association objects are retrieved even if there are no associations in the GWAS catalog 
		if(nrow(gwas_associations_by_efo_id@associations)==0)
		{
			efo_id<-NULL
		}
	}

	# it is redundant to specify both efo_trait and efo_id. 

	# if(!is.null(reported_trait) && !is.null(efo_trait))
	if(!is.null(reported_trait) && !is.null(efo_trait))
		return(gwasrapidd::union(gwas_associations_by_reported_trait, gwas_associations_by_efo_trait))

	# if(!is.null(reported_trait) && !is.null(efo_id) && is.null(efo_trait))
	if(!is.null(reported_trait) && !is.null(efo_id) && is.null(efo_trait))
		return(gwasrapidd::union(gwas_associations_by_reported_trait, gwas_associations_by_efo_id))

	# if(!is.null(reported_trait) && is.null(efo_trait) && is.null(efo_id))
	if(!is.null(reported_trait) && is.null(efo_trait) && is.null(efo_id))
		return(gwas_associations_by_reported_trait)

	# if(is.null(reported_trait) && !is.null(efo_trait))
	if(is.null(reported_trait) && !is.null(efo_trait))
		return(gwas_associations_by_efo_trait)

	# if(is.null(reported_trait) && is.null(efo_trait) && !is.null(efo_id))
	if(is.null(reported_trait) && is.null(efo_trait) && !is.null(efo_id))
		return(gwas_associations_by_efo_id)

	warning('no results found for the trait/efo in the GWAS catalog. Or perhaps you failed to specify `trait`, `efo` or `efo_id`? At least one of the three must be specified')
	return("no results found")
}
