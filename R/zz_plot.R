#' ZZ plot
#'
#' Calculate Z scores from the reported P values (Zp) and the reported log odds ratios (Zlnor). Construct a scatter plot of Zp and Zlnor
#'
#' @param dat the target dataset of interest
#' @param Title_size size of title
#' @param Title plot title
#' @param Ylab label for Y axis 
#' @param Xlab label for X axis
#' @param Title_xaxis_size size of x axis title
#'
#' @return plot 
#' @export

zz_plot<-function(dat=NULL,Title_size=12,Title="ZZ plot",Ylab="Z score inferred from p value",Xlab="Z score inferred from log odds ratio and standard error",Title_xaxis_size=12){
	dat<-dat[abs(dat$p)<=1,]
	dat$z.p<-stats::qnorm(dat$p/2,lower.tail=F)
	# dat[dat$z.p == "NaN", ]
	# dat$z<-dat$test_statistic
	dat$z.lnor<-abs(dat$lnor/dat$se)
	Colour<-"black"
	Cor<-stats::cor(dat$z.p,dat$z.lnor)
	if(Cor<0.99){
		Colour<-"red"
	}
	subtitle<-paste0("Pearson correlation coefficient=",round(Cor,2))
	# +ggplot2::labs(y= Ylab, x =Xlab,)
	# plot(dat$z.p,dat$z.lnor)
	dat$plot_x <- dat$z.p
	dat$plot_y <- dat$z.lnor
	plot<-ggplot2::ggplot(dat, ggplot2::aes(x=plot_x, y=plot_y)) + ggplot2::geom_point(colour=Colour) + ggplot2::ggtitle(Title) + ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+ggplot2::labs(y=Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size),plot.subtitle = ggplot2::element_text(size = 8))
	
	return(plot)
}
