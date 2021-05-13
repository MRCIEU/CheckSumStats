
#' MAF plot 
#'
#' Make a plot comparing minor allele frequency between outcome and reference studies.
#'
#' @param ref_dat user supplied reference dataset. data frame. optional 
#' @param ref_1000G if ref_dat is NULL, the user should indicate the 1000 genomes reference study of interest. options are: AFR, AMR, EAS, EUR, SAS or ALL. Default is to make plots for all super populations
#' @param snp_reference rsid column in ref_dat
#' @param target_dat  the outcome dataset of interest. Data frame. 
#' @param snp_target rsid column in target_dat
#' @param eaf name of the effect allele frequency column in target_dat 
#' @param ref_dat_maf name of the minor allele frequency column in the reference dataset.  Only necessary if ref_dat is specified
#' @param ref_dat_minor_allele name of the minor allele column in the reference dataset. Only necessary if ref_dat is specified
#' @param ref_dat_major_allele name of the major allele column in the reference dataset. Only necessary if ref_dat is specified 
#' @param target_dat_effect_allele name of the effect allele column in target_dat
#' @param target_dat_other_allele name of the non-effect allele column in target_dat 
#' @param target_dat_population population ancestry of target_dat 
#' @param ref_dat_population name of column describing population ancestry of reference dataset. Only necessary if ref_dat is specified
#' @param outcome name of the outcome trait corresponding to target_dat 
#' @param ID ID for the outcome trait of interest 
#' @param target_study column in target_dat indicating name of target study 
#' @param ref_study column in reference study indicating name of reference study #'. Only necessary if ref_dat is specified
#' @param Title_xaxis_size size of title on x axis #' 
#' @param Title_size size of plot title #' 
#' @param Title plot title 
#' @param Ylab Y label 
#' @param Xlab X label 
#' @param cowplot_title title of overall plot 
#'
#' @return plot 
#' @export

make_plot_maf<-function(ref_dat=NULL,ref_1000G=c("AFR","AMR", "EAS", "EUR", "SAS","ALL"),target_dat=NULL,eaf="eaf",snp_target="rsid",snp_reference="SNP",ref_dat_maf="MAF",target_dat_effect_allele="effect_allele",target_dat_other_allele="other_allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",ID=NULL,target_dat_population="population",ref_dat_population="population",target_study="study",ref_study="study",Title_xaxis_size=12,Title_size=12,Title="Comparison of MAF between outcome study and reference study",Ylab="Allele frequency in outcome study",Xlab="MAF in reference study",cowplot_title="Comparison of MAF between outcome study and 1000 genomes project"){	

	# should exclude palindromic SNPs which can cause apparent conflicts when target and reference datasets on different strands for some SNPs. drop palindromic SNPs or show in different shape?
	# ref_dat<-load_refdata(refstudy=refstudy,Dir=Dir)
	utils::data("refdat_1000G_superpops",envir =environment())
	if(is.null(ref_dat)){
		ref_dat<-refdat_1000G_superpops[refdat_1000G_superpops$population %in% ref_1000G,]
		ref_dat$study <- paste0("1000G ",ref_1000G)
		strand1<-c("A","G","T","C")
		strand2<-c("T","C","A","G")
	}

	if(!any(names(ref_dat) == "minor_allele")) stop("minor_allele column missing")
	if(!any(names(ref_dat) == "major_allele")) stop("major_allele column missing")


	ref_dat2<-flip_strand(dat=ref_dat,allele1_col="minor_allele",allele2_col="major_allele")
	ref_dat$minor_allele2<-ref_dat2$minor_allele
	ref_dat$major_allele2<-ref_dat2$major_allele

	# c("AFR","AMR","EAS","EUR","SAS", "ALL")

	names(target_dat)[names(target_dat) == target_dat_population]<-"target_dat_population"
	names(ref_dat)[names(ref_dat) == ref_dat_population]<-"ref_dat_population"
	names(target_dat)[names(target_dat) == target_study]<-"target_study"
	names(ref_dat)[names(ref_dat) == ref_study]<-"ref_study"
	names(ref_dat)[names(ref_dat) == ref_dat_maf]<-"maf_ref"

	if(any(names(ref_dat) %in% c(target_dat_effect_allele,target_dat_other_allele,target_dat_effect_allele))) warning("effect allele, other allele or eaf present in refererence dataset with same name as in target dataset")	

	dat.m<-merge(ref_dat,target_dat,by.x=snp_reference,by.y=snp_target)
	# dat.m[dat.m$SNP == "rs1298999",c("SNP",eaf,"maf_ref",target_dat_effect_allele,target_dat_other_allele,"minor_allele","major_allele","minor_allele2","major_allele2")]
	Pos<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	
	# harmonise study to ref dataset minor allele. need to clean this up as a new function with flip_strand and harmonise_allele functions like in the GWAS catalog functions
	dat.m[,eaf]<-as.numeric(dat.m[,eaf])
	dat.m[,"maf_ref"]<-as.numeric(dat.m[,"maf_ref"])
	dat.m[,eaf][Pos]<-1-dat.m[,eaf][Pos]
	EA<-dat.m[,target_dat_effect_allele][Pos]
	OA<-dat.m[,target_dat_other_allele][Pos]
	dat.m[,target_dat_effect_allele][Pos]<-OA
	dat.m[,target_dat_other_allele][Pos]<-EA
	
	# harmonise SNPs on different strands
	Pos1<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele]) #make the assumption that if the alleles still do not match they are on different strands
	Pos2<-which(dat.m[,target_dat_effect_allele] == dat.m[,ref_dat_minor_allele])
	dat.m1<-dat.m[Pos1,]
	dat.m2<-dat.m[Pos2,]
	
	# dat.m1[,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	Pos<-which(dat.m1[,target_dat_effect_allele] != dat.m1[,"minor_allele2"])
	# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	dat.m1[,eaf][Pos]<-1-dat.m1[,eaf][Pos]
	EA<-dat.m1[,target_dat_effect_allele][Pos]
	OA<-dat.m1[,target_dat_other_allele][Pos]
	dat.m1[,target_dat_effect_allele][Pos]<-OA
	dat.m1[,target_dat_other_allele][Pos]<-EA
	Pos1<-which(dat.m1[,target_dat_effect_allele] != dat.m1[,"minor_allele2"])
	if(sum(Pos1) > 0 ) stop("there are still mismatched alleles after flipping the strands")

	dat.m<-rbind(dat.m1,dat.m2)
	
	# dat.m[,c("SNP",eaf,"maf_ref",target_dat_effect_allele,"minor_allele","minor_allele2")]

	# dat.m[,c("Effect.Allele","Other.Allele","minor_allele","major_allele","eaf","maf")]

# Pos<-which(dat.m1$Effect.Allele  != dat.m1$minor_allele2 &  dat.m1$Effect.Allele  != dat.m1$minor_allele)
# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele","minor_allele2")]
	dat.m$alleles<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	dat.m$alleles2<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	
	dat.m.test<-dat.m
	outcome_plot<-outcome
	if(is.null(outcome)){
		outcome_plot<-""
	}

	if(!is.null(outcome) & !is.null(target_study) & !is.null(ID)){
		outfile_name<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID=",dat.m.test[,ID]))
		outfile_name<-gsub("\\|","",outfile_name)
		outfile_name<-gsub(" ","_",outfile_name)
		outfile_name<-gsub("__","_",outfile_name)
		outfile_name<-gsub("=","",outfile_name)
		outfile_name<-gsub("/","_",outfile_name)
		outcome_plot<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		target_study<-unique(paste0(dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		# outcome_plot2<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test[,study]," | ID: ",dat.m.test[,ID]))
	}

	Plot_list<-NULL
	dat.m.test<-dat.m.test[order(dat.m.test$ref_dat_population),]
	Pops<-unique(dat.m.test$ref_dat_population)
	# if(length(Pops)>1){
		Ylab<-""
	# }
	for(i in 1:length(Pops)){
		# print(pop)
		# dat1<-dat[dat$ref_dat_population==pop,]			
		dat1<-dat.m.test[dat.m.test$ref_dat_population==Pops[i], ]
		pop2<-c("European","East Asian","African","American","South Asian","Global")
		Pops2<-c("EUR","EAS","AFR","AMR","SAS","ALL")
		j<-which(Pops2 %in% Pops[i])
		# if(Xlab==""){
		# Xlab<-paste0(pop2[j]," MAF")
		Xlab<-pop2[j]

		# }
		# Title<-pop2[i]
		# Title<-gsub("ALL","Global pop",Title)
		# if(Title == ""){
		# 	Title<-target_study
		# }

		# if(subtitle_off) Title<-""

		Colour<-rep("black",nrow(dat1))
		Colour[which(dat1[,eaf]>0.5)]<-"blue"
		# Colour[dat1[,target_dat_effect_allele]!=dat1[,ref_dat_minor_allele]]<-"red"	
		
		# fix harmonisation functions above so that efffect allele strand flipped to minor_allele (not harmonised with minor_allele2)
		Colour[dat1$Effect.Allele!=dat1$minor_allele & dat1$Effect.Allele!=dat1$minor_allele2]<-"red"	

		Diff<-abs(dat1[,eaf]-dat1[,"maf_ref"])
		Colour[which(Diff>0.10)]<-"red"
		Shape<-rep(19,nrow(dat1))
		Shape[which(dat1$alleles %in% c("AT","TA","GC","CG"))]<-1
		dat1$eaf<-dat1[,eaf]
		dat1$maf<-dat1[,"maf_ref"]

		Title_size1<-Title_size
		Subtitle_size1<-8
		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			Title_size1<-0
			Subtitle_size1<-0
		}

		
		# Temp<-dat1[dat1$eaf > 0.5,c("Effect.Allele","minor_allele","Other.Allele","major_allele","eaf","maf","alleles")]
		# dim(Temp)
		# length(which(Temp$alleles %in% c("CG","GC","TA","AT")))

		# snps<-dat1$rsid[dat1$eaf>0.55]
		# dat1[dat1$eaf>0.55,c("minor_allele","major_allele","Effect.Allele","Other.Allele","maf","eaf")]
		# dat2<-dat[dat$rsid %in% snps, ]
		# dat2[,c("rsid","Effect.Allele","Other.Allele","eaf")]
		# head(ref_dat)
		# dat2.m<-merge(dat2,ref_dat,by="rsid")
		# head(dat2.m[,c("rsid","Effect.Allele","Other.Allele","eaf", "maf", "minor_allele2","major_allele2" )])
		
		Subtitle<-unique(paste0("Reported population: ",dat1$target_dat_population))
		
		Plot<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(colour=Colour) +ggplot2::ggtitle(Title) +ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size1,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = Subtitle_size1)) 

			if(length(Pops)> 1){
				Plot_list[[i]]<-Plot
			}
		}	

		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			if(is.null(cowplot_title)){
				cowplot_title<-target_study
			}	
			Plot<-make_cow_plot(Plot_list=Plot_list,Title=cowplot_title,Xlab="",Ylab="",return_plot=TRUE,Title_size=Title_size,Subtitle=Subtitle,Subtitle_size=8)
		}
		# Title_axis_size
	return(Plot)
}


flip_strand<-function(dat=NULL,allele1_col=NULL,allele2_col=NULL){
	# Pos<-dat[,allele1]!=dat[,allele2]	
	strand1<-c("A","T","G","C")
	strand2<-c("T","A","C","G")
	# lnor.y<-dat$lnor.y[Pos]*-1
	# dat$lnor.y[Pos]<-lnor.y
	allele1<-dat[,allele1_col]
	if(!is.null(allele2_col)){
		allele2<-dat[,allele2_col]
		dat[,allele2_col]<-strand2[match(allele2,strand1)]				
	}
	dat[,allele1_col]<-strand2[match(allele1,strand1)]	
	return(dat)
}

make_cow_plot<-function(Plot_list=NULL,Title="",Xlab="",Ylab="",out_file=NULL,return_plot=TRUE,width=1000,height=1000,Title_size=0,Title_axis_size=10,Subtitle="",Subtitle_size=0){
	
	Plot<-cowplot::plot_grid(plotlist=Plot_list)

	if(Title!="") { 
		title <- cowplot::ggdraw() + 
				cowplot::draw_label(
					Title,
					# fontface = 'bold',
					fontface = 'plain',
					x = 0,
					hjust = 0,
					size=Title_size)  +
				ggplot2::theme(
				# add margin on the left of the drawing canvas,
				# so title is aligned with left edge of first plot
					plot.margin = ggplot2::margin(0, 0, 0, 7)
					)

		subtitle <- cowplot::ggdraw() + 
		cowplot::draw_label(
			Subtitle,
			# fontface = 'bold',
			fontface = 'plain',
			x = 0,
			hjust = 0,
			# element = "plot.subtitle",
			size=Subtitle_size)  +
		ggplot2::theme(
		# add margin on the left of the drawing canvas,
		# so title is aligned with left edge of first plot
			plot.margin = ggplot2::margin(0, 0, 0, 7)
			)
		# subtitle <- ggdraw() +
  # 						draw_label_theme("By census tract, 2016",
  #                  theme = theme_georgia(), 
  #                  element = "plot.subtitle",

  #                  x = 0.05, hjust = 0, vjust = 1)

		Plot<-cowplot::plot_grid(title,subtitle, Plot,ncol = 1,rel_heights = c(0.05,0.05, 1))
	}
	y.grob <- grid::textGrob(Ylab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size), rot=90)

	x.grob <- grid::textGrob(Xlab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size))

	if(!return_plot){
		grDevices::png(out_file, width = width, height = height)
			gridExtra::grid.arrange(gridExtra::arrangeGrob(Plot, left = y.grob, bottom = x.grob))
		grDevices::dev.off()	
	}
	if(return_plot){
		return(Plot)
	}

}
