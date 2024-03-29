#' Plot comparing the test study to the GWAS catalog
#'
#' Make a plot comparing signed Z scores, or effect allele frequency, between the test dataset and the GWAS catalog, in order to identify effect allele meta data errors 
#'
#' @param dat the test dataset of interest
#' @param beta name of the column containing the SNP effect size
#' @param se name of the column containing the standard error for the SNP effect size. 
#' @param plot_type compare Z scores or effect allele frequency? For comparison of Z scores set plot_type to "plot_zscores". For comparison of effect allele frequency set to "plot_eaf". Default is set to "plot_zscores"
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param gwas_catalog_ancestral_group restrict the comparison to these ancestral groups in the GWAS catalog. Default is set to (c("European","East Asian") 
#' @param force_all_trait_study_hits force the plot to include GWAS hits from the outcome study if they are not in the GWAS catalog? This should be set to TRUE only if dat is restricted to GWAS hits for the trait of interest. This is useful for visualising whether the outcome/trait study has an unusually larger number of GWAS hits, which could, in turn, indicate that the summary statistics have not been adequately cleaned.
#' @param exclude_palindromic_snps should the function exclude palindromic SNPs? default set to TRUE. If set to FALSE, then conflicts with the GWAS catalog could reflect comparison of different reference strands. 
#' @param distance_threshold distance threshold for deciding if the GWAS hit in the test dataset is present in the GWAS catalog. For example, a distance_threshold of 25000 means that the GWAS hit in the test dataset must be within 25000 base pairs of a GWAS catalog association, otherwise it is reported as missing from the GWAS catalog.  
#' @param map_association_to_study map associations to study in GWAS catalog. This supports matching of results on PMID and study ancestry, which increases accuracy of comparisons, but is slow when there are large numbers of associations. Default = FALSE
#' @param gwas_catalog user supplied data frame containing results from the GWAS catalog for the trait of interest. If set to NULL then the function will retrieve results from the GWAS catalog.  
#' @param gc_dat output of compare_effect_to_gwascatalog2. This will typically be ignored by most users. Default NULL
#' @param legend include legend in plot. Default TRUE
#' @param Title plot title
#' @param Ylab label for Y axis 
#' @param Xlab label for X axis
#' @param nocolour if TRUE, effect size conflicts are illustrated using shapes rather than colours. Default FALSE
#' @param publication_quality produce a high resolution image e.g. for publication purposes. Default FALSE
#' @param return_dat if TRUE, the dataset used to generate the plot is returned to the user and no plot is made. 
#'
#' @return plot 
#' @export

make_plot_gwas_catalog<-function(dat=NULL,plot_type="plot_zscores",efo_id=NULL,efo=NULL,trait=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),legend=TRUE,Title="Comparison of Z scores between test dataset & GWAS catalog",Ylab="Z score in test dataset",Xlab="Z score in GWAS catalog",force_all_trait_study_hits=FALSE,exclude_palindromic_snps=TRUE,beta="beta",se="se",distance_threshold=25000,return_dat=FALSE,map_association_to_study=FALSE,gwas_catalog=NULL,nocolour=FALSE,publication_quality=FALSE,gc_dat=NULL){
	
	if(is.null(gc_dat)){
	Dat.m<-compare_effect_to_gwascatalog2(dat=dat,beta=beta,se=se,efo_id=efo_id,efo=efo,trait=trait,force_all_trait_study_hits=force_all_trait_study_hits,exclude_palindromic_snps=exclude_palindromic_snps,distance_threshold=distance_threshold,gwas_catalog=gwas_catalog,map_association_to_study=map_association_to_study )
	}else{
		Dat.m<-gc_dat
	}
	
	Names<-grep("beta",names(Dat.m))
	Names2<-grep("effect",names(Dat.m))
	Names3<-grep("eaf",names(Dat.m))
	Names4<-grep("ances",names(Dat.m))
	
	Dat.m$Z_scores[Dat.m$Z_scores=="high conflict"]<-"red"
	Dat.m$Z_scores[Dat.m$Z_scores=="not present in the GWAS catalog"]<-"red"
	Dat.m$Z_scores[Dat.m$Z_scores=="moderate conflict"]<-"blue"
	Dat.m$Z_scores[Dat.m$Z_scores=="no conflict"]<-"black"
	labels_colour<-unique(Dat.m$Z_scores)
	values_colour<-unique(Dat.m$Z_scores)
	Dat.m$plot_x<-Dat.m$z.x
	Dat.m$plot_y<-Dat.m$z.y
	Dat.m$colour<-Dat.m$Z_scores
	Name<-"Effect size conflict"

	if(plot_type=="plot_eaf")
	{
		Dat.m<-Dat.m[!is.na(Dat.m$eaf.x),]
		Dat.m$EAF[Dat.m$EAF=="high conflict"]<-"red"
		Dat.m$EAF[Dat.m$EAF=="moderate conflict"]<-"blue"
		Dat.m$EAF[Dat.m$EAF=="no conflict"]<-"black"
		labels_colour<-unique(Dat.m$EAF)
		values_colour<-unique(Dat.m$EAF)
		Dat.m$plot_x<-Dat.m$eaf.x
		Dat.m$plot_y<-Dat.m$eaf.y
		Dat.m$colour<-Dat.m$EAF
		Name<-"EAF conflict"
		Ylab="EAF in outcome study"
		Xlab="EAF in GWAS catalog"
		Title="Comparison of EAF between test dataset and GWAS catalog"
	}

	if(return_dat) return(Dat.m)
	
	labels_colour[labels_colour == "red"]<-"high"
	if(force_all_trait_study_hits & any(Dat.m$z.x ==0)) 
	{
		labels_colour[labels_colour == "high"]<-"high or not\npresent in GWAS catalog"
	}
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	ancestry1<-Dat.m$ancestral_group
	if(is.null(ancestry1))
	{
		ancestry1<-"Unknown"	
	}
	
	labels_shape<-unique(ancestry1)[order(unique(ancestry1))]

	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)
	if(any(is.na(values_shape))) {
		values_shape[is.na(values_shape)]<-17
	}
	
	# values_shape<-c(16,15,17,18)
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}

	################### 
	#nocolour argument: colour and shape values and labels for when nocolour argument is TRUE
	plot_shape_values<-Dat.m$Z_scores
	plot_shape_values[plot_shape_values=="red"]<-3
	plot_shape_values[plot_shape_values=="blue"]<-2
	plot_shape_values[plot_shape_values=="black"]<-1
	shape_map<-data.frame(cbind(c(3,2,1),c("high","moderate","none")))
	
	names(shape_map)<-c("values","labels")
	if(force_all_trait_study_hits & any(Dat.m$z.x ==0)){
		shape_map$labels[shape_map$labels == "high"]<-"high or not\npresent in GWAS catalog"
	}
	shape_map<-shape_map[order(shape_map$values),]
	shape_values<-as.numeric(shape_map$values)
	shape_labels<-shape_map$labels				

	colour_map<-data.frame(cbind(c("gray","black"),c("East Asian","European")))
	names(colour_map)<-c("values","labels")
	colour_map<-colour_map[order(colour_map$labels),]
	colour_values<-colour_map$values
	colour_labels<-colour_map$labels	
	
	Subtitle<-unique(paste0(Dat.m$trait," | reported ancestry in test dataset: ",Dat.m$population))

	Title_size1<-20
	Subtitle_size1<-10
	Legend_title_size1<-20
	Legend_text_size1<-10
	Axis.text_size1<-10
	Axis_title_size_x1<-10
	Axis_title_size_y1<-10		
	geom_point_size1<-2
	shape_width<-1	

	if(publication_quality){
		Title_size1<-50
		Subtitle_size1<-40
		Legend_title_size1<-32
		Legend_text_size1<-32
		Axis.text_size1<-32
		Axis_title_size_x1<-50
		Axis_title_size_y1<-50			
		geom_point_size1<-20
		shape_width<-3
	} 

	my_theme<-ggplot2::theme(
		plot.title = ggplot2::element_text(size = Title_size1,hjust = 0),
		plot.subtitle = ggplot2::element_text(size =Subtitle_size1),
		axis.title.x=ggplot2::element_text(size=Axis_title_size_x1),
		axis.title.y=ggplot2::element_text(size=Axis_title_size_y1),
		axis.text=ggplot2::element_text(size=Axis.text_size1),
		legend.title=ggplot2::element_text(size=Legend_title_size1),
		legend.text=ggplot2::element_text(size=Legend_text_size1))



	if(legend){
		Plot<-ggplot2::ggplot(Dat.m) + 
			ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=colour,shape=ancestry1),size=geom_point_size1,stroke=shape_width) + 
			ggplot2::ggtitle(Title) +
			ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + 						
			ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		                     labels = labels_shape,
		                     # labels = unique(ancestry1)[order(unique(ancestry1))],
		                     # labels = c("European","East Asian"),
		                     values = values_shape) + 
		                     # values = 1:length(Shape2)) + 
		 	ggplot2::scale_colour_manual(name=Name,
			              labels=labels_colour,
			              values=values_colour)+
		 	my_theme
		 	# ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"),
				# )+
		 	# ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		 	# ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
			# ggplot2::theme(legend.text=ggplot2::element_text(size=8))
			
			if(nocolour){					
				Plot<-ggplot2::ggplot(Dat.m) + 
					ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=ancestry1,shape=plot_shape_values),size=geom_point_size1,stroke=shape_width) +
					ggplot2::ggtitle(Title) +
					ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + 				
				ggplot2::scale_shape_manual(name = "Effect size conflict",
		                     labels = shape_labels,
		                     values = shape_values) + 		                  
		 		ggplot2::scale_colour_manual(name="GWAS catalog ancestry",
			              labels=colour_labels,
			              values=colour_values)+
		 		my_theme
		 		# ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		 		# ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		 		# ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
				# ggplot2::theme(legend.text=ggplot2::element_text(size=8))

			}
		}

	if(!legend){
		Plot<-ggplot2::ggplot(Dat.m) +
		ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=colour,shape=ancestry1),size=geom_point_size1,stroke=shape_width) +
		ggplot2::ggtitle(Title) +
		ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + 
		 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
	                    labels = labels_shape,	                     
	                     values = values_shape) + 
	 	ggplot2::scale_colour_manual(name=Name,
		              labels=labels_colour,
		              values=values_colour)+
	 	my_theme+
	 	ggplot2::theme(legend.position = "none")
	 	
	 	# ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		# ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
	 	# ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 	# 	legend.text=ggplot2::element_text(size=8),plot.subtitle = ggplot2::element_text(size = 8),
	 	# 	legend.position = "none")
	 	
	 	if(nocolour){					
				Plot<-ggplot2::ggplot(Dat.m) + 
					ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=ancestry1,shape=plot_shape_values),size=geom_point_size1,stroke=shape_width) +
					ggplot2::ggtitle(Title) +
					ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + 					
				ggplot2::scale_shape_manual(name = "Effect size conflict",
		                     labels = shape_labels,
		                     values = shape_values) + 		                  
		 		ggplot2::scale_colour_manual(name="GWAS catalog ancestry",
			              labels=colour_labels,
			              values=colour_values)+
				my_theme+
 				ggplot2::theme(legend.position = "none")
	 
		
		 		# ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
				# ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		}
	} 
  	
	return(Plot)
}


#' Compare the genetic effect sizes in the test dataset to the GWAS catalog
#'
#' Compare the direction of effects and effect allele frequency between the test dataset and the GWAS catalog, in order to identify effect allele meta data errors
#'
#' @param dat the test dataset of interest
#' @param beta name of the column containing the SNP effect size
#' @param se name of the column containing the standard error for the SNP effect size. 
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param gwas_catalog_ancestral_group restrict the comparison to these ancestral groups in the GWAS catalog. Default is set to (c("European","East Asian") 
#' @param force_all_trait_study_hits force the comparison to include GWAS hits from the test dataset if they are not in the GWAS catalog? This should be set to TRUE only if dat is restricted to GWAS hits for the trait of interest. This is useful for visualising whether the test trait study has an unusually larger number of GWAS hits, which could, in turn, indicate analytical issues with the summary statistics 
#' @param exclude_palindromic_snps should the function exclude palindromic SNPs? default set to TRUE. If set to FALSE, then conflicts with the GWAS catalog could reflect comparison of different reference strands. 
#' @param distance_threshold distance threshold for deciding if the GWAS hit in the test dataset is present in the GWAS catalog. For example, a distance_threshold of 25000 means that the GWAS hit in the test dataset must be within 25000 base pairs of a GWAS catalog association, otherwise it is reported as missing from the GWAS catalog.   
#'
#' @return dataframe
#' @export


compare_effect_to_gwascatalog<-function(dat=NULL,efo=NULL,efo_id=NULL,trait=NULL,beta=NULL,se=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),exclude_palindromic_snps=TRUE,force_all_trait_study_hits=FALSE,distance_threshold=distance_threshold)
{
	# exclude the MAF 1k ref set. Causes problems if you force inclusion of SNPs missing from the GWAS catalog 
	utils::data("refdat_1000G_superpops",envir =environment())
	snps_exclude<-unique(refdat_1000G_superpops$SNP)
	dat<-dat[!dat$rsid %in% snps_exclude,]

	if(beta=="lnor")
	{
		if(!"lnor" %in% names(dat)) stop("name of beta column set to lnor but there is no column with that name")
	}

	if(!beta %in% names(dat)) stop(paste0("beta column not found. Check you correctly specified the name of the beta column"))

	if(!se %in% names(dat)) stop(paste0("se column not found. Check you correctly specified the name of the se column"))

	if(is.null(efo) & is.null(efo_id) & is.null(trait)) stop("you must specify either efo, efo_id or trait")
	
	gwas_catalog<-gwas_catalog_hits(efo=efo,efo_id=efo_id,trait=trait)

	message_trait<-paste(c(efo,efo_id,trait),collapse="/")
	Dat.m<-merge(gwas_catalog,dat,by="rsid")	

	if(all(is.na(Dat.m$effect_allele.x))) stop(paste0("associations for ",message_trait," were found but all effect alleles are missing in the GWAS catalog. Therefore no comparison of effect size direction can be made"))
	Dat.m<-Dat.m[!is.na(Dat.m$effect_allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$effect_allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$other_allele)==1,]
	Alleles<-paste0(Dat.m$effect_allele.y,Dat.m$other_allele)
	if(exclude_palindromic_snps)
	{
		Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	}
	if(!is.null(gwas_catalog_ancestral_group))
	{
		# c("European","East Asian")
		Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]	
	}
	# Dat.m1<-Dat.m
	# Dat.m<-Dat.m1
	Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y	
	if(any(Pos)) 
	{
		Dat.m1<-Dat.m[Pos,]
		Dat.m2<-Dat.m[!Pos,]
		Dat.m1<-flip_strand(dat=Dat.m1,allele1_col="effect_allele.x")
		# Dat.m1$effect_allele.x
		# Dat.m1$effect_allele.y
		# Dat.m1[,c("effect_allele.x","effect_allele.y","other_allele","rsid")]
		Dat.m<-rbind(Dat.m1,Dat.m2)
	}
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y
	if(any(Pos))
	{
		Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	}	
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y

	if(any(Pos)) 
	{
		stop("effect alleles not fully harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m[,beta]/Dat.m[,se] 
	
	
	

	if("pmid" %in% names(dat))
	{
		gwas_studies<-gwasrapidd::get_studies(study_id=unique(Dat.m$study_id ))
		Publications<-gwas_studies@publications
		Publications<-Publications[!duplicated(Publications$study_id),]
		Dat.m<-merge(Dat.m,Publications,by="study_id")
	}


	#identifty eaf conflicts
	# ancestry2<-Dat.m$ancestral_group	
	Dat.m$EAF<-"no conflict"
	Dat.m$EAF[is.na(Dat.m$eaf.x)]<-NA
	# EAF<-rep("black",nrow(Dat.m))
	Pos1<-which(Dat.m$eaf.x<0.5 & Dat.m$eaf.y>0.5 | Dat.m$eaf.x>0.5 & Dat.m$eaf.y<0.5)
	Dat.m$EAF[Pos1]<-"moderate conflict"	 
	Pos2<-which(Dat.m$eaf.x<0.40 & Dat.m$eaf.y>0.60 | Dat.m$eaf.x>0.60 & Dat.m$eaf.y<0.40)
	Dat.m$EAF[Pos2]<-"high conflict"
	Pos3<-which(Dat.m$pmid==Dat.m$pubmed_id)
	Pos4<-Pos1[Pos1 %in% Pos3] 
	Dat.m$EAF[Pos4]<-"high conflict" #if there is a moderate eaf conflict (eaf close to 0.5) but both datasets are from the same study, then the conflict is upgraded to high

	# if(plot_type=="plot_zscores"){
	if(force_all_trait_study_hits)
	{
		gc_list<-find_hits_in_gwas_catalog(gwas_hits=dat$rsid,trait=trait,efo=efo,efo_id=efo_id,distance_threshold=distance_threshold)
			
		if(length(gc_list$not_in_gc)>0)
		{
		# if(any(!dat$rsid %in% gwas_catalog$rsid)){
			# dat$rsid[!dat$rsid %in% gwas_catalog$rsid]
			dat2<-dat[dat$rsid %in% gc_list$not_in_gc,] #the snps not in the GWAS catalog. Genomic coordinates for SNPs associated with trait/efo in the GWAS catalog did not overlap with these SNPs (including +/- 250 kb) 
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

	Dat.m$Z_scores<-"no conflict"
	# Z_scores<-rep("black",nrow(Dat.m))
	Dat.m$Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"moderate conflict"
	Dat.m$Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >= 3.890592 & abs(Dat.m$z.x) >= 3.890592 )]<-"high conflict" # Z score of 3.890592 = 2 sided p value of 0.0001	
	Dat.m$Z_scores[which(Dat.m$pmid==Dat.m$pubmed_id & sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"high conflict" #if the signs are different but Z.x and Z.y come from the same study, then there is a clear incompatability
	if(force_all_trait_study_hits){
		Dat.m$Z_scores[Dat.m$z.x==0]<-"high conflict" #these SNPs are not in the GWAS catalog
	}
		# Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >=  4.891638  & abs(Dat.m$z.x) >=  4.891638 )]<-"red"
	return(Dat.m)
}

#' Compare the genetic effect sizes in the test dataset to the GWAS catalog
#'
#' Compare the direction of effects and effect allele frequency between the test dataset and the GWAS catalog, in order to identify effect allele meta data errors
#'
#' @param dat the test dataset of interest
#' @param beta name of the column containing the SNP effect size
#' @param se name of the column containing the standard error for the SNP effect size. 
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param gwas_catalog_ancestral_group restrict the comparison to these ancestral groups in the GWAS catalog. Default is set to (c("European","East Asian") 
#' @param exclude_palindromic_snps should the function exclude palindromic SNPs? default set to TRUE. If set to FALSE, then conflicts with the GWAS catalog could reflect comparison of different reference strands. 
#' @param map_association_to_study map associations to study in GWAS catalog. This supports matching of results on PMID and study ancestry, which increases accuracy of comparisons, but is slow when there are large numbers of associations. Default = FALSE.
#' @param gwas_catalog user supplied data frame containing results from the GWAS catalog for the trait of interest. If set to NULL then the function will retrieve results from the GWAS catalog. 
#' @param force_all_trait_study_hits force the comparison to include GWAS hits from the test dataset if they are not in the GWAS catalog? This should be set to TRUE only if dat is restricted to GWAS hits for the trait of interest. This is useful for visualising whether the test trait study has an unusually larger number of GWAS hits, which could, in turn, indicate analytical issues with the summary statistics 
#' @param distance_threshold distance threshold for deciding if the GWAS hit in the test dataset is present in the GWAS catalog. For example, a distance_threshold of 25000 means that the GWAS hit in the test dataset must be within 25000 base pairs of a GWAS catalog association, otherwise it is reported as missing from the GWAS catalog.  
#'
#' @return dataframe
#' @export


compare_effect_to_gwascatalog2<-function(dat=NULL,efo=NULL,efo_id=NULL,trait=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),exclude_palindromic_snps=TRUE,map_association_to_study=FALSE,beta="beta",se="se",gwas_catalog=NULL,force_all_trait_study_hits=FALSE,distance_threshold=distance_threshold)
{
	
	if(!beta %in% names(dat)) stop("could not find the effect size column; please specify the name of your effect size column using the beta argument")

	if(is.null(gwas_catalog))
	{
		gwas_catalog<-gwas_catalog_hits(efo=efo,efo_id=efo_id,trait=trait,map_association_to_study=map_association_to_study)
	}
	
	if(length(gwas_catalog) ==1)
	{
		if(gwas_catalog=="no results found") return(gwas_catalog)
	}
	
	message_trait<-paste(c(efo,efo_id,trait),collapse="/")
	Dat.m<-merge(gwas_catalog,dat,by="rsid")

	if(exclude_palindromic_snps)
	{
		Alleles<-paste0(Dat.m$effect_allele.y,Dat.m$other_allele)
		Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	}	


	if(all(is.na(Dat.m$effect_allele.x)) | nrow(Dat.m)==0) return(paste0("associations for ",message_trait," were found in the GWAS catalog but all effect alleles were missing or all SNPs were palindromic. Therefore no comparison of effect size direction could be made"))
	
	Dat.m<-Dat.m[!is.na(Dat.m$effect_allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$effect_allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$other_allele)==1,]
		
	if(!is.null(gwas_catalog_ancestral_group) & map_association_to_study)
	{
		if(!is.null(Dat.m$ancestral_group))
		{
			Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]	
		}
	}
	# Dat.m1<-Dat.m
	# Dat.m<-Dat.m1
	Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y	
	if(any(Pos)) 
	{
		Dat.m1<-Dat.m[Pos,]
		Dat.m2<-Dat.m[!Pos,]
		Dat.m1<-flip_strand(dat=Dat.m1,allele1_col="effect_allele.x")
		# Dat.m1$effect_allele.x
		# Dat.m1$effect_allele.y
		# Dat.m1[,c("effect_allele.x","effect_allele.y","other_allele","rsid")]
		Dat.m<-rbind(Dat.m1,Dat.m2)
		if(nrow(Dat.m)==0) return("no SNP associations found in GWAS catalog")
	}
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y
	if(any(Pos))
	{
		Dat.m<-harmonise_effect_allele(dat=Dat.m,beta=beta)
	}	
	Pos<-Dat.m$effect_allele.x!=Dat.m$effect_allele.y

	if(any(Pos)) 
	{
		stop("effect alleles not fully harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m[,beta]/Dat.m[,se] 
	
	if("pmid" %in% names(dat) & map_association_to_study)
	{
		gwas_studies<-gwasrapidd::get_studies(study_id=unique(Dat.m$study_id ))
		Publications<-gwas_studies@publications
		Publications<-Publications[!duplicated(Publications$study_id),]
		Dat.m<-merge(Dat.m,Publications,by="study_id")
	}
	#identifty eaf conflicts
	# ancestry2<-Dat.m$ancestral_group	
	Dat.m$EAF<-"no conflict"
	Dat.m$EAF[is.na(Dat.m$eaf.x) | is.na(Dat.m$eaf.y)]<-NA
	# EAF<-rep("black",nrow(Dat.m))
	Pos1<-which(Dat.m$eaf.x<0.5 & Dat.m$eaf.y>0.5 | Dat.m$eaf.x>0.5 & Dat.m$eaf.y<0.5)
	Dat.m$EAF[Pos1]<-"moderate conflict"	 
	Pos2<-which(Dat.m$eaf.x<0.40 & Dat.m$eaf.y>0.60 | Dat.m$eaf.x>0.60 & Dat.m$eaf.y<0.40)
	Dat.m$EAF[Pos2]<-"high conflict"
	Pos3<-which(Dat.m$pmid==Dat.m$pubmed_id)
	Pos4<-Pos1[Pos1 %in% Pos3] 
	Dat.m$EAF[Pos4]<-"high conflict" #if there is a moderate eaf conflict (eaf close to 0.5) but both datasets are from the same study, then the conflict is upgraded to high

	if(force_all_trait_study_hits)
	{#this step assumes that dat contains GWAS hits defined in the test dataset 
		gc_list<-find_hits_in_gwas_catalog(gwas_hits=dat$rsid,trait=trait,efo=efo,efo_id=efo_id,distance_threshold=distance_threshold)
			
		if(length(gc_list$not_in_gc)>0)
		{
		# if(any(!dat$rsid %in% gwas_catalog$rsid)){
			# dat$rsid[!dat$rsid %in% gwas_catalog$rsid]
			dat2<-dat[dat$rsid %in% gc_list$not_in_gc,] #the snps not in the GWAS catalog. Genomic coordinates for SNPs associated with trait/efo in the GWAS catalog did not overlap with these SNPs (including +/- 250 kb) 
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

	Dat.m$Z_scores<-"no conflict"
	Dat.m$Z_scores[is.na(Dat.m$z.x) | is.na(Dat.m$z.y)]<-NA

	# Z_scores<-rep("black",nrow(Dat.m))
	Dat.m$Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"moderate conflict"
	Dat.m$Z_scores[which(as.numeric(Dat.m$z.x)==0) ]<-"not present in the GWAS catalog"
	Dat.m$Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >= 3.890592 & abs(Dat.m$z.x) >= 3.890592 )]<-"high conflict" # Z score of 3.890592 = 2 sided p value of 0.0001	
	Dat.m$Z_scores[which(Dat.m$pmid==Dat.m$pubmed_id & sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"high conflict" #if the signs are different but Z.x and Z.y come from the same study, then there is a clear incompatability	
	return(Dat.m)
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


#' Are hits in the GWAS catalog? 
#'
#' Identify GWAS hits in the test dataset and see if they overlap with GWAS hits in the GWAS catalog. 
#'
#' @param gwas_hits the "GWAS hits" in the test dataset (e.g. SNP-trait associations with P<5e-8)
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param distance_threshold distance threshold for deciding if the GWAS hit in the test dataset is present in the GWAS catalog. For example, a distance_threshold of 25000 means that the GWAS hit in the test dataset must be within 25000 base pairs of a GWAS catalog association, otherwise it is reported as missing from the GWAS catalog.   
#'
#' @return list 
#' @export

find_hits_in_gwas_catalog<-function(gwas_hits=NULL,trait=NULL,efo=NULL,efo_id=NULL,distance_threshold=25000){

	utils::data("refdat_1000G_superpops",envir =environment())
	snps_exclude<-unique(refdat_1000G_superpops$SNP)
	gwas_hits<-gwas_hits[!gwas_hits %in% snps_exclude]

	ensembl<-get_positions_biomart(gwas_hits=gwas_hits)
	if(!is.null(efo)) efo<-trimws(unlist(strsplit(efo,split=";")))	
	if(!is.null(efo_id)) efo_id<-trimws(unlist(strsplit(efo_id,split=";")))	
	if(!is.null(trait)) trait<-trimws(unlist(strsplit(trait,split=";")))	
 
	gwas_variants<-get_gwas_associations(reported_trait=trait,efo_trait=efo,efo_id=efo_id)	

	if(class(gwas_variants) != "associations")
	{
		if(gwas_variants == "no results found") return(gwas_variants)
	}

	
	# gwas_variants<-gwasrapidd::get_variants(efo_trait = efo,efo_id=efo_id,reported_trait=trait)		
	if(class(unlist(gwas_variants)) == "character")
	{
		if(nrow(gwas_variants)==0)
		{
			warning(paste("search returned 0 variants from the GWAS catalog"))
		}
	}

	
	if(is.null(trait) & is.null(efo) & is.null(efo_id))
	{
		genomic_range<-list(chromosome=as.character(ensembl$chr_name),start=ensembl$chrom_start - distance_threshold,end=ensembl$chrom_start + distance_threshold)
		gwas_variants<-gwasrapidd::get_variants(genomic_range=genomic_range)
		gwas_variants<-data.frame(gwas_variants@variants)
		
		ens.m<-merge(ensembl,gwas_variants,by.x="chr_name",by.y="chromosome_name",all.x=TRUE)
		
		Pos<-abs(ens.m$chrom_start.x-ens.m$chrom_start.y)<distance_threshold
		# Pos<-which(ens.m$chromosome_position>ens.m$bp_minus & ens.m$chromosome_position<ens.m$bp_plus) 
		gwashit_in_gc<-unique(ens.m$refsnp_id[Pos])
		gwashit_notin_gc<-unique(ens.m$refsnp_id[!ens.m$refsnp_id %in% gwashit_in_gc])
		return(list("not_in_gc"=gwashit_notin_gc,"in_gc"=gwashit_in_gc))
	}
	
	if(!(is.null(trait) & is.null(efo) & is.null(efo_id))){
	
		# for now use ensembl/biomart to determine positions for GWAS catalog and test variants. Both are in GRCh38 so could also use GWAS catalog positions for GWAS catalog variats (maybe this would be faster too) but there is the risk that the reference build could diverge over time between biomart/ensembl and GWAS catalog. might update this so that chromosome positions could be based on GWAS catalog instead	
		# if(positions_biomart)
		# {
		# gwas_variants<-data.frame(gwas_variants@variants)
		
		# gwas_hits %in% gwas_variants@variants$variant_id 
		# ensembl2<-get_positions_biomart(gwas_hits=unique(gwas_variants$variant_id))

		ensembl2<-get_positions_biomart(gwas_hits=gwas_variants@risk_alleles$variant_id)
			# }
		gwashit_in_gc<-NA
		if(any(ensembl$chr_name %in% ensembl2$chr_name))
		{
			gwashit_notin_gc<-ensembl$refsnp_id[!ensembl$chr_name %in% ensembl2$chr_name]
			ens.m<-merge(ensembl,ensembl2,by="chr_name")	
			# ens.m[which(ens.m$refsnp_id.x  =="rs12239737"),c("chrom_start.x","chrom_start.y")]

			Test<-any(abs(ens.m$chrom_start.x-ens.m$chrom_start.y)<distance_threshold)
			if(Test)
			{
				Pos<-abs(ens.m$chrom_start.x-ens.m$chrom_start.y)<distance_threshold
				# Pos<-ens.m$chrom_start.x>ens.m$bp_minus.y & ens.m$chrom_start.x<ens.m$bp_plus.y		
				gwashit_in_gc<-unique(ens.m$refsnp_id.x[Pos])
				ens.m<-ens.m[!ens.m$refsnp_id.x %in% gwashit_in_gc,]
				Pos<-abs(ens.m$chrom_start.x-ens.m$chrom_start.y)<distance_threshold			
				# Pos<-ens.m$chrom_start.x>ens.m$bp_minus.y & ens.m$chrom_start.x<ens.m$bp_plus.y		
				gwashit_notin_gc<-c(gwashit_notin_gc,unique(ens.m$refsnp_id.x[!Pos]))
			}
			if(!Test)
			{
				gwashit_notin_gc<-c(gwashit_notin_gc,unique(ens.m$refsnp_id.x))
			}
		}else{
			gwashit_notin_gc<-unique(ensembl$refsnp_id)
			gwashit_in_gc<-NA
		}
		return(list("not_in_gc"=gwashit_notin_gc,"in_gc"=gwashit_in_gc))
	}
}

get_positions_biomart<-function(gwas_hits=NULL){
	# library(biomaRt)

	# Get chromosomal positions and genes names from ENSEMBL. Should be build 38. Version object contains version ID for genome build used
	Mart <- biomaRt::useMart(host="https://www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Version<-biomaRt::listDatasets(Mart)[ biomaRt::listDatasets(Mart)$dataset=="hsapiens_snp","version"]
	message(paste0("Using ",Version," of human genome from ensembl for genomic coordinates"))
	Attr<-biomaRt::listAttributes(Mart)

	ensembl<-biomaRt::getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters="snp_filter",values=gwas_hits,mart=Mart)
	ensembl<-ensembl[order(ensembl$refsnp_id),]
	ensembl<-ensembl[nchar(ensembl$chr_name)<3,]
	ensembl$chr_name<-as.numeric(ensembl$chr_name)
	# ensembl$bp_minus<-ensembl$chrom_start - bp_down
	# ensembl$bp_plus<-ensembl$chrom_start + bp_up
	return(ensembl)
}

#' Flag conflicts with the GWAS catalog
#'
#' Flag conflicts with the GWAS catalog through comparison of reported effect alleles and reported effect allele frequency.  
#'
#' @param dat the test dataset of interest
#' @param beta name of the column containing the SNP effect size
#' @param se name of the column containing the standard error for the SNP effect size. 
#' @param trait the trait of interest
#' @param efo_id ID for trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param gwas_catalog_ancestral_group restrict the comparison to these ancestral groups in the GWAS catalog. Default is set to (c("European","East Asian") 
#' @param exclude_palindromic_snps should the function exclude palindromic SNPs? default set to TRUE. If set to FALSE, then conflicts with the GWAS catalog could reflect comparison of different reference strands. 
#'
#' @return list
#' @export

flag_gc_conflicts<-function(dat=NULL,beta="lnor",se="lnor_se",efo=NULL,trait=NULL,efo_id=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),exclude_palindromic_snps=TRUE){
	
	gc_dat<-compare_effect_to_gwascatalog(dat=dat,efo=efo,trait=trait,efo_id=efo_id,beta=beta,se=se,gwas_catalog_ancestral_group=gwas_catalog_ancestral_group,exclude_palindromic_snps=exclude_palindromic_snps)

	effect_size_conflict<-gc_dat$Z_scores
	gc_conflicts<-c("high conflict","moderate conflict","no conflict")
	es_conflicts_list<-lapply(1:length(gc_conflicts),FUN=function(x)
		length(effect_size_conflict[which(effect_size_conflict==gc_conflicts[x])]))
	total<-length(which(!is.na(gc_dat$Z_scores)))
	es_conflicts_list<-c(es_conflicts_list,total)
	names(es_conflicts_list)<-c(gc_conflicts,"n_snps")

	eaf_conflicts<-gc_dat$EAF
	eaf_conflicts_list<-lapply(1:length(gc_conflicts),FUN=function(x)
		length(eaf_conflicts[which(eaf_conflicts==gc_conflicts[x])]))
	total<-length(which(!is.na(gc_dat$EAF)))
	eaf_conflicts_list<-c(eaf_conflicts_list,total)
	names(eaf_conflicts_list)<-c(gc_conflicts,"n_snps")	

	# gc_ancestries<-paste(unique(gc_dat$ancestral_group),collapse="; ")	

	all_conflicts_list<-list("effect_size_conflicts"=es_conflicts_list,"eaf_conflicts"=eaf_conflicts_list)
	return(all_conflicts_list)
}


#' Flag conflicts with the GWAS catalog
#'
#' Flag conflicts with the GWAS catalog through comparison of reported effect alleles and reported effect allele frequency.  
#'
#' @param gc_dat dataset generated by compare_effect_to_gwascatalog2()
#'
#' @return list
#' @export

flag_gc_conflicts2<-function(gc_dat=NULL){
	effect_size_conflict<-gc_dat$Z_scores
	gc_conflicts<-c("high conflict","moderate conflict","no conflict")
	es_conflicts_list<-lapply(1:length(gc_conflicts),FUN=function(x)
		length(effect_size_conflict[which(effect_size_conflict==gc_conflicts[x])]))
	total<-length(which(!is.na(gc_dat$Z_scores)))
	es_conflicts_list<-c(es_conflicts_list,total)
	names(es_conflicts_list)<-c(gc_conflicts,"n_snps")

	eaf_conflicts<-gc_dat$EAF
	eaf_conflicts_list<-lapply(1:length(gc_conflicts),FUN=function(x)
		length(eaf_conflicts[which(eaf_conflicts==gc_conflicts[x])]))
	total<-length(which(!is.na(gc_dat$EAF)))
	eaf_conflicts_list<-c(eaf_conflicts_list,total)
	names(eaf_conflicts_list)<-c(gc_conflicts,"n_snps")	

	# gc_ancestries<-paste(unique(gc_dat$ancestral_group),collapse="; ")	

	all_conflicts_list<-list("effect_size_conflicts"=es_conflicts_list,"eaf_conflicts"=eaf_conflicts_list)
	return(all_conflicts_list)
}

