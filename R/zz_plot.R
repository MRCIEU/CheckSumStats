#' ZZ plot
#'
#' Calculate Z scores from the reported P values (Zp) and the reported log odds ratios (Zlnor). Construct a scatter plot of Zp and Zlnor
#'
#' @param dat the target dataset of interest
#' @param beta the name of the column containing the SNP effect size
#' @param se the name of the column containing the standard error for the SNP effect size
#' @param Title plot title
#' @param Ylab label for Y axis 
#' @param Xlab label for X axis
#' @param exclude_1000G_MAF_refdat exclude rsids from the 1000 genome MAF reference dataset. 
#'
#' @return plot 
#' @export

zz_plot<-function(dat=NULL,Title="ZZ plot",Ylab="Z score inferred from p value",Xlab="Z score inferred from effect size and standard error",beta="lnor",se="lnor_se",exclude_1000G_MAF_refdat=TRUE){

	if(exclude_1000G_MAF_refdat){
		utils::data("refdat_1000G_superpops",envir =environment())
		snps_exclude<-unique(refdat_1000G_superpops$SNP)
		dat<-dat[!dat$rsid %in% snps_exclude,]
	}

	dat<-dat[abs(dat$p)<=1,]
	dat$z.p<-stats::qnorm(dat$p/2,lower.tail=F)
	dat<-dat[dat$z.p!="Inf" ,]
	# dat[dat$z.p == "NaN", ]
	# dat$z<-dat$test_statistic
	dat$z.beta<-abs(dat[,beta]/dat[,se])
	Colour<-"black"
	Cor<-stats::cor(dat$z.p,dat$z.beta)
	if(Cor<0.99){
		Colour<-"red"
	}
	subtitle<-paste0("Pearson correlation coefficient=",round(Cor,3))
	# +ggplot2::labs(y= Ylab, x =Xlab,)
	# plot(dat$z.p,dat$z.beta)
	dat$plot_x <- dat$z.p
	dat$plot_y <- dat$z.beta
	my_theme<-ggplot2::theme(
		plot.title = ggplot2::element_text(size = 50,hjust = 0),
		plot.subtitle = ggplot2::element_text(size =40),
		axis.title.x=ggplot2::element_text(size=50),
		axis.title.y=ggplot2::element_text(size=50),
		axis.text=ggplot2::element_text(size=32),
		legend.title=ggplot2::element_text(size=32),
		legend.text=ggplot2::element_text(size=32))

	geom_point_size1<-20

	plot<-ggplot2::ggplot(dat, ggplot2::aes(x=plot_x, y=plot_y)) + 
		ggplot2::geom_point(colour=Colour,size=geom_point_size1) + 
		ggplot2::ggtitle(Title) + 		
		ggplot2::labs(y=Ylab, x =Xlab,subtitle=subtitle) + 
		my_theme
		# ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		# ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size),plot.subtitle = ggplot2::element_text(size = 8))
	
	return(plot)
}
