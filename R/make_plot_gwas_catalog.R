
#' Compare the target study to the GWAS catalog
#'
#' Make a plot comparing signed Z scores, or effect allele frequency, between the target dataset and the GWAS catalog 
#'
#' @param dat the outcome dataset of interest
#' @param beta name of the column containing the SNP effect size
#' @param se name of the column containing the standard error for the SNP effect size. 
#' @param plot_type compare Z scores or effect allele frequency? For comparison of Z scores set plot_type to "plot_zscores". For comparison of effect allele frequency set to "plot_eaf". Default is set to "plot_zscores"
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param gwas_catalog_ancestral_group restrict the comparison to these ancestral groups in the GWAS catalog. Default is set to (c("European","East Asian") 
#' @param force_all_trait_study_hits force the plot to include GWAS hits from the outcome study if they are not in the GWAS catalog? This should be set to TRUE only if dat is restricted to GWAS hits for the trait of interest. This is useful for visualising whether the outcome/trait study has an unusually larger number of GWAS hits, which could, in turn, indicate that the summary statistics have not been adequately cleaned.
#' @param exclude_palindromic_snps should the function exclude palindromic SNPs? default set to TRUE. If set to FALSE, then conflicts with the GWAS catalog could reflect comparison of different reference strands. 
#' @param legend include legend in plot. Default TRUE
#' @param Title plot title
#' @param Title_size_subplot size of title 
#' @param Ylab label for Y axis 
#' @param Xlab label for X axis
#' @param Title_xaxis_size size of x axis title
#'
#' @return plot 
#' @export

make_plot_gwas_catalog<-function(dat=NULL,plot_type="plot_zscores",efo_id=NULL,efo=NULL,trait=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),legend=TRUE,Title="Comparison of Z scores between outcome study and GWAS catalog",Title_size_subplot=12,Ylab="Z score in outcome study",Xlab="Z score in GWAS catalog",Title_xaxis_size=12,force_all_trait_study_hits=FALSE,exclude_palindromic_snps=TRUE,beta="lnor",se="lnor_se"){

	# exclude the MAF 1k ref set. Causes problems if you force inclusion of SNPs missing from the GWAS catalog 
	utils::data("refdat_1000G_superpops",envir =environment())
	snps_exclude<-unique(refdat_1000G_superpops$SNP)
	dat<-dat[!dat$rsid %in% snps_exclude,]
	
	if(beta=="lnor"){
		if(!"lnor" %in% names(dat)) stop("name of beta column set to lnor but there is no column with that name")
	}

	if(!beta %in% names(dat)) stop(paste0("beta column set to '",beta,"' but there is no column with that name"))

	if(!se %in% names(dat)) stop(paste0("se column set to '",se,"' but there is no column with that name"))

	# if(!is.null(trait)){
		# gwas_catalog<-gwas_catalog_hits(trait=trait)
	# }
# gwas_catalog_ancestral_group="East Asian"
	if(!is.null(efo)){
		gwas_catalog<-gwas_catalog_hits(efo=efo)
		# gwas_catalog$ancestral_group
	}
	
	if(!is.null(efo_id)){
		gwas_catalog<-gwas_catalog_hits(efo_id=efo_id)
	}

	if(!is.null(efo) | !is.null(efo_id)){
		if(!is.null(trait)){
			gwas_catalog1<-gwas_catalog_hits(trait=trait)
			# gwas_catalog1$efo<-"trait"
			# gwas_catalog$efo<-"efo"
			gwas_catalog<-rbind(gwas_catalog,gwas_catalog1)
			gwas_catalog<-gwas_catalog[which(!duplicated(gwas_catalog[,c("study_id","rsid","test_statistic")])),]
		}
	}
	
	if(is.null(efo) & is.null(efo_id)){
		if(!is.null(trait)){
			gwas_catalog<-gwas_catalog_hits(trait=trait)
		}else{
			stop("efo_id, efo and trait are all null")
		}
	}

	# gwas_catalog[gwas_catalog$study_id == "GCST001633",]
	# gwas_catalog[which(gwas_catalog$rsid == "rs2736100"),]

	Dat.m<-merge(gwas_catalog,dat,by="rsid")	
	Dat.m<-Dat.m[!is.na(Dat.m$effect_allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$effect_allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$other_allele)==1,]
	Alleles<-paste0(Dat.m$effect_allele.y,Dat.m$other_allele)
	if(exclude_palindromic_snps){
		Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	}
	if(!is.null(gwas_catalog_ancestral_group)){
		# c("European","East Asian")
		Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]	
	}
	# Dat.m1<-Dat.m
	# Dat.m<-Dat.m1
	Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y	
	if(any(Pos)) {
		Dat.m1<-Dat.m[Pos,]
		Dat.m2<-Dat.m[!Pos,]
		Dat.m1<-flip_strand(dat=Dat.m1,allele1_col="effect_allele.x")
		# Dat.m1$effect_allele.x
		# Dat.m1$effect_allele.y
		# Dat.m1[,c("effect_allele.x","effect_allele.y","other_allele","rsid")]
		Dat.m<-rbind(Dat.m1,Dat.m2)
	}
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y
	if(any(Pos)){
		Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	}	
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y

	if(any(Pos)) {
		stop("effect alleles not fully harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m[,beta]/Dat.m[,se] 
	# Dat.m$z.x<-Dat.m$beta_gc/Dat.m$se_gc

	# Dat.m$z.y<-Dat.m$lnor.y/Dat.m$se.y
	# Dat.m$z.x<-Dat.m$lnor.x/Dat.m$se.x

	# head(Dat.m[,c("p.x","z.x","p.y","z.y")])
	# max(Dat.m$p.x)
	# dim(Dat.m)
	# Ylab<-""
	# Xlab<-""

	if("pmid" %in% names(dat)){
		gwas_studies<-gwasrapidd::get_studies(study_id=unique(Dat.m$study_id ))
		Publications<-gwas_studies@publications
		Publications<-Publications[!duplicated(Publications$study_id),]
		Dat.m<-merge(Dat.m,Publications,by="study_id")
	}

	if(plot_type=="plot_eaf"){
		Dat.m<-Dat.m[!is.na(Dat.m$eaf.x),]
		# ancestry2<-Dat.m$ancestral_group
		Pos1<-which(Dat.m$eaf.x<0.5 & Dat.m$eaf.y>0.5 | Dat.m$eaf.x>0.5 & Dat.m$eaf.y<0.5)
		EAF<-rep("black",nrow(Dat.m))
		EAF[Pos1]<-"blue"	 
		Pos2<-which(Dat.m$eaf.x<0.40 & Dat.m$eaf.y>0.60 | Dat.m$eaf.x>0.60 & Dat.m$eaf.y<0.40)
		EAF[Pos2]<-"red"
		Pos3<-which(Dat.m$pmid==Dat.m$pubmed_id)
		Pos4<-Pos1[Pos1 %in% Pos3] 
		EAF[Pos4]<-"red" #if there is a moderate eaf conflict (eaf close to 0.5) but both datasets are from the same study, then the conflict is upgraded to high

		labels_colour<-unique(EAF)
		values_colour<-unique(EAF)
		Dat.m$plot_x<-Dat.m$eaf.x
		Dat.m$plot_y<-Dat.m$eaf.y
		Dat.m$colour<-EAF 
		Name<-"EAF conflict"
		Ylab="EAF in outcome study"
		Xlab="EAF in GWAS catalog"
		Title="Comparison of EAF between outcome study and GWAS catalog"
	}

	if(plot_type=="plot_zscores"){
		if(force_all_trait_study_hits){
			if(any(!dat$rsid %in% gwas_catalog$rsid)){
				dat$rsid[!dat$rsid %in% gwas_catalog$rsid]
				dat2<-dat[!dat$rsid %in% gwas_catalog$rsid,]
				Dat.m2<-merge(gwas_catalog,dat2,by="rsid",all.y=TRUE)
				Dat.m2$z.y<-Dat.m2[,beta]/Dat.m2[,se]
				Dat.m2$z.x<-0
				# Dat.m$plot_x
				Dat.m2$ancestral_group<-unique(dat$population)
				Names<-names(Dat.m)[!names(Dat.m) %in% names(Dat.m2)]
				for(i in 1:length(Names)){
					Dat.m2[,Names[i]]<-NA
				}
				# Dat.m3<-Dat.m
				Dat.m<-rbind(Dat.m,Dat.m2)
			}
		}

		Z_scores<-rep("black",nrow(Dat.m))
		Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"blue"
		Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >= 3.890592 & abs(Dat.m$z.x) >= 3.890592 )]<-"red" # Z score of 3.890592 = 2 sided p value of 0.0001	
		Z_scores[which(Dat.m$pmid==Dat.m$pubmed_id & sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"red" #if the signs are different but Z.x and Z.y come from the same study, then there is a clear incompatability
		if(force_all_trait_study_hits){
			Z_scores[Dat.m$z.x==0]<-"red" #these SNPs are not in the GWAS catalog
		}
		# Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >=  4.891638  & abs(Dat.m$z.x) >=  4.891638 )]<-"red"
		labels_colour<-unique(Z_scores)
		values_colour<-unique(Z_scores)
		Dat.m$plot_x<-Dat.m$z.x
		Dat.m$plot_y<-Dat.m$z.y
		Dat.m$colour<-Z_scores
		Name<-"Effect size conflict"
	}

	labels_colour[labels_colour == "red"]<-"high"
	if(force_all_trait_study_hits & any(!dat$rsid %in% gwas_catalog$rsid)) {
		labels_colour[labels_colour == "high"]<-"high or not\npresent in GWAS catalog"
	}
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	# dim(unique(Dat.m[,c("rsid","pmid")]))
	# Dat.m2<Dat.m
	# length(unique(Dat.m$rsid))
	# unique(Dat.m$pmid)
	# Dat.m<-Dat.m[grep(";",Dat.m$ancestral_group,invert=T),]
	# Dat.m<-Dat.m[Dat.n$ancestral_group == "European"]
	ancestry1<-Dat.m$ancestral_group
	
	labels_shape<-unique(ancestry1)[order(unique(ancestry1))]
	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)
	# values_shape<-c(16,15,17,18)
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}

	
	# unique(ancestry1)[order(unique(ancestry1))]
	# 1:length(unique(ancestry1))
	# ancestry2<-c("European")
	
	# Shape<-ancestry1
	# Shape[Shape=="European"]<-15
	# Shape[Shape=="East Asian"]<-16
	# Shape<-as.numeric(Shape)
	# # Shape<-unique(Shape)[order(unique(Shape))]

	
	# Shape2<-Shape
	# Shape2<-unique(Shape2)[order(unique(Shape2))]
	

	Subtitle<-paste0(Dat.m$outcome," | ",Dat.m$population)

	# ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1))
	if(legend){
			Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=colour,shape=ancestry1)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"),
				)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
			 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		                     labels = labels_shape,
		                     # labels = unique(ancestry1)[order(unique(ancestry1))],
		                     # labels = c("European","East Asian"),
		                     values = values_shape) + 
		                     # values = 1:length(Shape2)) + 
		 	ggplot2::scale_colour_manual(name=Name,
			              labels=labels_colour,
			              values=values_colour)+
		 	ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
			ggplot2::theme(legend.text=ggplot2::element_text(size=8))
		}

	if(!legend){
		Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=colour,shape=ancestry1)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
	                    labels = labels_shape,	                     
	                     values = values_shape) + 
	 	ggplot2::scale_colour_manual(name=Name,
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8),plot.subtitle = ggplot2::element_text(size = 8),
	 		legend.position = "none")
	}
	 # ggplot2::scale_colour_manual(name="Z score conflict",
  #                     labels=unique(Z_scores)[order(unique(Z_scores))] ,
  #                     values=unique(Z_scores)[order(unique(Z_scores))]) 	 

	  # ggplot2::scale_colour_manual(name="Z score conflict",
                      # labels=c("none", "moderate","high"),
                      # values=c("black","blue", "red")) 	 
  	 
  	
	return(Plot)
}

harmonise_effect_allele<-function(dat=NULL,beta=beta){
	Pos<-which(dat$effect_allele.x!=dat$effect_allele.y)
	beta.y<-dat[,beta][Pos]*-1
	dat[,beta][Pos]<-beta.y
	oa<-dat$effect_allele.y[Pos]
	ea<-dat$other_allele[Pos]
	dat$effect_allele.y[Pos]<-ea
	dat$other_allele[Pos]<-oa
	eaf<-1-dat$eaf.y[Pos]
	dat$eaf.y[Pos]<-eaf		
	return(dat)
}