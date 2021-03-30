#' Make cow plot
#'
#' Combine all plots into a single plot using the cowplot package
#'
#' @param Plot_list list of plots to combine
#' @param out_file filepath to save the plot
#' @param return_plot logical argument. If TRUE, plot is returned and is not save to out_file 
#' @param width width of plot
#' @param height height of plot 
#' @param bycols logical argument. If true, plots by column
#' @param Title plot title
#' @param Xlab label for X axis
#' @param Ylab label for Y axis
#' @param Title_size size of title
#' @param Title_axis_size size of x axis title
#' @param ncol number of columns
#' @param nrow number of rows
#'
#' @return plot 
#' @export

make_cow_plot2<-function(Plot_list=NULL,out_file=NULL,return_plot=FALSE,width=1000,height=1000,Title="",Xlab="",Ylab="",Title_size=0,Title_axis_size=10,bycols=TRUE,ncol=2,nrow=4){

	# Plot<-cowplot::plot_grid(plotlist=Plot_list[[1]])
	if(bycols){
		Plot<-cowplot::plot_grid(plotlist=Plot_list,nrow=nrow,ncol=ncol)
	}
	if(!bycols){
		Plot<-cowplot::plot_grid(plotlist=Plot_list)
	}
	if(Title!="") { 
		title <- cowplot::ggdraw() + 
				cowplot::draw_label(
					Title,
					fontface = 'bold',
					# fontface = 'plain',
					x = 0,
					hjust = 0,
					size=Title_size)  +
				ggplot2::theme(
				# add margin on the left of the drawing canvas,
				# so title is aligned with left edge of first plot
					plot.margin = ggplot2::margin(0, 0, 0, 7)
					)

		Plot<-cowplot::plot_grid(title, Plot,ncol = 1,rel_heights = c(0.05, 1))
	}

	y.grob <- grid::textGrob(Ylab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size), rot=90)

	x.grob <- grid::textGrob(Xlab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size))

	if(!return_plot){
		grDevices::png(out_file, width = width, height = height,)
			gridExtra::grid.arrange(gridExtra::arrangeGrob(Plot, left = y.grob, bottom = x.grob))
		grDevices::dev.off()	
	}
	if(return_plot){
		return(Plot)
	}

}
