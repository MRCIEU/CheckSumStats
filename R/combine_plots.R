#' Make cow plot
#'
#' Combine all plots into a single plot using the cowplot package
#'
#' @param Plot_list plots to combine. Can either be vector of character strings giving the names of plot objects or a list of plot objects. 
#' @param out_file filepath to save the plot
#' @param return_plot logical argument. If TRUE, plot is returned and is not save to out_file 
#' @param width width of plot
#' @param height height of plot 
#' @param by2cols logical argument. If true, forces plot to have 2 columns
#' @param Title plot title
#' @param Xlab label for X axis
#' @param Ylab label for Y axis
#' @param Title_size size of title
#' @param Title_axis_size size of x axis title
#' @param Ncol number of columns
#'
#' @return plot 
#' @export

combine_plots<-function(Plot_list=NULL,out_file=NULL,return_plot=FALSE,width=800,height=1000,Title="",Xlab="",Ylab="",Title_size=0,Title_axis_size=0,by2cols=TRUE,Ncol=2){

	if(is.null(Plot_list)){
		Plot_list2<-ls()[grep("Plot[0-9]",ls())] 
		Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
	}
		
	if(is.character(Plot_list)){
		Plot_list<-lapply(1:length(Plot_list),FUN=function(x) eval(parse(text=Plot_list2[x])))
	}	

	# Plot<-cowplot::plot_grid(plotlist=Plot_list[[1]])
	if(by2cols){
		Nrow<-round(length(Plot_list)/2)
		Plot<-cowplot::plot_grid(plotlist=Plot_list,nrow=Nrow,ncol=Ncol)
	}
	if(!by2cols){
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

make_plot_list<-function(){
	Plot_list2<-ls()[grep("Plot[0-9]",ls())] 
	Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))	
	return(Plot_list)
}
