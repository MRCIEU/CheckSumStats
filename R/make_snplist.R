make_snplist<-function(trait=NULL,efo_id=NULL,efo=NULL,population=NULL,Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/"){

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
	
	# snplist<-c(snplist,top_hits_rsids)
	load(paste0(Dir,"refdat_1000G_superpops.Rdata"))
	snplist<-c(top_hits_rsids,unique(refdat_1000G_superpops$SNP))	
	snplist<-unique(snplist)		
	
	return(snplist)

}


