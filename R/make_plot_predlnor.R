
#' Predicted versus reported effect sizes
#'
#' Make a plot comparing the predicted effect sizes to the reported effect sizes. 
#'
#' @param dat the target dataset of interest
#' @param pred_beta name of column containing the predicted effect size
#' @param pred_beta_se name of column containing the standard error for the predicted effect size
#' @param beta name of column containing the reported effect size
#' @param se name of column containing the standard error for the reported effect size
#' @param sd_est the standard deviation of the phenotypic mean. Can either be a numeric vector of length 1 or name of the column in dat containing the standard deviation value (in which case should be constant across SNPs). Only applicable for continuous traits. If not supplied by the user, the standard deviation is approximated using sd_est, estimated by the predict_beta_sd() function. 
#' @param bias logical argument. If TRUE, plots the % deviation of the predicted from the reported effect size on the Y axis against the reported effect size on the X axis.  
#' @param subtitle subtitle
#' @param maf_filter minor allele frequency threshold. If not NULL, genetic variants with a minor allele frequency below this threshold are excluded
#' @param legend logical argument. If true, includes figure legend in plot
#' @param Title_size size of title
#' @param Title plot title
#' @param Ylab label for Y axis 
#' @param Xlab label for X axis
#' @param Title_xaxis_size size of x axis title
#' @param standard_errors logical argument. If TRUE, plots the predicted versus the reported standard errors for the effect sizes
#'
#' @return plot 
#' @export

make_plot_pred_effect<-function(dat=NULL,Xlab="Reported effect size",Ylab="Predicted effect size",subtitle="",maf_filter=FALSE,bias=FALSE,Title_size=12,Title="Predicted versus reported effect size",Title_xaxis_size=12,legend=TRUE,standard_errors=FALSE,pred_beta="lnor_pred",pred_beta_se="lnor_se_pred",beta="lnor",se="lnor_se",sd_est="sd_est"){

	if("ncase" %in% names(dat)){
		dat<-format_data_predlnor_sh(dat=dat)
	}
	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))

	utils::data("refdat_1000G_superpops",envir =environment())
	snps_exclude<-unique(refdat_1000G_superpops$SNP)
	dat<-dat[!dat$rsid %in% snps_exclude,]

	# summary(dat$bias)

	# for(i in 1:ncol(dat)){
	# 	dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
	# }
	
	# dat<-dat[complete.cases(dat),]
	
	if(is.null(Title)){
	# if(!is.null(outcome_name)){
		Title<-outcome_name
	}

	if(maf_filter){
		maf<-dat$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		dat<-dat[maf>=maf_filter,]
	}
	
	MAF<-rep("black",nrow(dat))
	MAF<-dat$eaf
	MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	MAF[MAF<=0.10]<-"0.01-0.10"
	MAF[MAF>0.10 & MAF<=0.20]<-"0.11-0.20"
	MAF[MAF>0.20 & MAF<=0.30]<-"0.21-0.30"
	MAF[MAF>0.30 & MAF<=0.40]<-"0.31-0.40"
	MAF[MAF>0.40 & MAF<=0.50]<-"0.41-0.50"
	Shape<-rep(19,nrow(dat))

	dat$plot_y<-dat[,pred_beta]
	dat$plot_x<-dat[,beta]
	
	if(sd_est %in% names(dat)){
		dat$plot_x<-dat[,beta]/dat[,sd_est]
	}

	if(pred_beta!="lnor_pred"){
		if(!sd_est %in% names(dat) ){
			if(is.na(as.numeric(sd_est))){
				warning("no SD value supplied")
			}
			if(!is.na(as.numeric(sd_est))){
				if(length(unique(sd_est))!=1) warning("more than one SD value has been supplied when only one is expected")
				dat$plot_x<-dat[,beta]/sd_est
			}
		}
	}
	
	dat$bias<-(dat$plot_y-dat$plot_x )/dat$plot_x*100	

	# Xlab<-"Reported log odds ratio"
	# Ylab<-"Predicted log odds ratio"
	if(standard_errors){
		dat$plot_y<-dat[,pred_beta_se]
		dat$plot_x<-dat[,se]
		Xlab<-"Reported standard error"
		Ylab<-"Predicted standard error"
	}

	if(bias){
		dat$plot_y<-dat$bias
		Med<-round(summary(dat$bias)[3],1)
		p25<-round(summary(dat$bias)[2],1)
		p75<-round(summary(dat$bias)[5],1)
		Min<-round(summary(dat$bias)[1],1)
		Max<-round(summary(dat$bias)[6],1)
		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
		Ylab<-"% deviation of predicted from reported effect size"
	}

	# dat$X<-dat[,1]
	# dat$Y<-dat[,2]
	# if(bias){	
	# 	linear_regression<-FALSE
	# }

	Values<-c("red","orange","purple","blue","black")
	Labels<-c("0.01-0.10","0.11-0.20","0.21-0.30","0.31-0.40","0.41-0.50")

	if(!bias){	
		Model<-summary(stats::lm(plot_y~plot_x,dat))
		int<-Model$coefficients[1,1]
		slope<-Model$coefficients[2,1]
		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
	}		


	if(!legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=maf))+ggplot2::theme(legend.position = "none")  + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))
		# +
		# ggplot2::scale_colour_manual(name = "maf",
	 #        labels = Labels,
	 #        values = Values)
	}

	if(legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=plot_x, y=plot_y,colour=maf))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))
				# +
				# ggplot2::scale_colour_manual(name = "maf",
			        # labels = Labels,
			        # values = Values)
	}

	return(Plot)
}

# make_plot_predlnor2<-function(dat=NULL,Xlab="",Ylab="",linear_regression=TRUE,subtitle="",maf_filter=NULL,bias=FALSE,Title_size=0,Title=NULL,Title_xaxis_size=10){

# 	# outfile_name<-unique(paste(dat$outcome,dat$study,dat$ID,sep="_"))
# 	# outfile_name<-gsub(" ","_",outfile_name)
# 	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))

# 	dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
# 	dat$bias<-(dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs*100	
# 	# Bias1<-((dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs)*100	
# 	# Bias2<-(dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs*100	

# 	for(i in 1:ncol(dat)){
# 		dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
# 	}

# 	dat<-dat[complete.cases(dat),]
	
# 	if(is.null(Title)){
# 	# if(!is.null(outcome_name)){
# 		Title<-outcome_name
# 	}

# 	if(!is.null(maf_filter)){
# 		maf<-dat$eaf
# 		maf[maf>0.5]<-1-maf[maf>0.5]
# 		dat<-dat[maf>maf_filter,]
# 	}
	
# 	Colour<-rep("black",nrow(dat))
# 	Colour<-dat$eaf
# 	Colour[Colour>0.5]<-1-Colour[Colour>0.5]
# 	Colour[Colour<=0.10]<-"red"
# 	Colour[Colour>0.10 & Colour<=0.20]<-"orange"
# 	Colour[Colour>0.20 & Colour<=0.30]<-"green"
# 	Colour[Colour>0.30 & Colour<=0.40]<-"blue"
# 	Colour[Colour>0.40 & Colour<=0.50]<-"black"
# 	Shape<-rep(19,nrow(dat))

# 	dat$Y<-dat$lnor_pred2
# 	dat$X<-dat$lnor_obs

# 	if(bias){
# 		dat$X<-dat$bias
# 		Med<-round(summary(dat$bias)[3],1)
# 		p25<-round(summary(dat$bias)[2],1)
# 		p75<-round(summary(dat$bias)[5],1)
# 		Min<-round(summary(dat$bias)[1],1)
# 		Max<-round(summary(dat$bias)[6],1)
# 		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
# 	}

# 	# dat$X<-dat[,1]
# 	# dat$Y<-dat[,2]
# 	if(bias){	
# 		linear_regression<-FALSE
# 	}

# 	if(linear_regression){	
# 		Model<-summary(lm(Y~X,dat))
# 		int<-Model$coefficients[1,1]
# 		slope<-Model$coefficients[2,1]
# 		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
# 	}
# 		# Title_xaxis_size
# 	Plot<-ggplot2::ggplot(dat, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_point(colour=Colour,shape=Shape) + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
# 	ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
# 	ggplot2::scale_colour_manual(name = "MAF",
# 		# labels = c("0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
#         # values = c("orange","yellow","blue","black"))
#         labels = c("0-0.10","0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
#         values = c("red","orange","yellow","blue","black"))
# 	return(Plot)
# }


# predict_lnor<-function(lnor=NULL,se=NULL,n=NULL,p=NULL,cases=NULL,controls=NULL){
# 	t.stat<-lnor/se
# 	lnor_pred1<-sqrt(t.stat^2 / (t.stat^2 + n) / 2 / p*(1-p) * pi^2/3)
# 		  # sqrt(t^2 / (t^2 + n) / 2 / p*(1-p) * pi^2/3)
# 	lnor_pred1<-lnor_pred1*sign(t.stat)	

# # 	exposed cases = (1-(1-maf)^2)*ncases
# # unexposed cases = (1-maf)^2*ncases
# # exposed controls = (1-(1-maf)^2)*ncontrols
# # unexposed controls = (1-maf)^2*ncontrols

# 	# n1<-(p^2+2*p*(1-p))*cases
# 	n1<- (1-(1-p)^2)*cases #exposed cases
# 	n2<- (1-p)^2*cases #unexposed cases
# 	n3<- (1-(1-p)^2)*controls #exposed controls
# 	n4<- (1-p)^2*controls #unexposed controls
# 	lnse_pred2<-sqrt(1/n1+1/n2+1/n3+1/n4)	
# 	lnor_pred2<-t.stat*lnse_pred2
# 	lnor_pred2<-lnor_pred2/1.4 #method 2 overestimates log odds ratio by 40% on average
# 	Dat<-data.frame(do.call(cbind,list(lnor_obs=lnor,lnor_pred1=lnor_pred1,lnor_pred2=lnor_pred2,eaf=p)))
# 	return(Dat)
# }


format_data_predlnor_sh<-function(dat=NULL){
	dat$MAC_case<-dat$maf*dat$ncase*2
	dat$MAC_control<-dat$maf*dat$ncontrol*2
	dat<-dat[dat$MAC_case>=50 & dat$MAC_control>=50,] 
	dat<-dat[dat$lnor_pred <= 1.999 & dat$lnor_pred>= -1.999,] #lnor_sh ==1.999 is an artifiact. 
	return(dat)
}


#' Predicted log odds ratio
#'
#' Predict the log odds ratio, using the Harrison approach. https://seanharrisonblog.com/2020/. The log odds ratio is inferred from the reported number of cases and controls, Z scores and minor allele frequency 
#'
#' @param dat the outcome dataset of interest
#'
#' @return data frame
#' @export
	
predict_lnor_sh<-function(dat=NULL){
	# load("~/harmonised_data.Rdata")
	# dat<-Dat
	# dat<-Dat[!is.na(Dat$eaf),]
	# Pos<-nrow(dat)/2
	# dat1<-dat[1:Pos,]
	# Pos<-Pos+1
	# dat<-dat[Pos:nrow(dat),]
	if(!any(names(dat) == "z")){
		dat$z<-dat$lnor/dat$lnor_se
	}
	
	log_or<-NULL
	log_or_se<-NULL
	# ID<-NULL
	# snp<-NULL
	snp_n<-nrow(dat)
	
	if(!"maf" %in% names(dat)){
		maf<-dat$eaf
		Pos<-which(maf>0.5)
		maf[Pos]<-1-maf[Pos]
		dat$maf<-maf
	}
	
	for(i in 1:snp_n)
	# for(i in 1:10)
	{
		# i<-1
		print(paste0("Analysing SNP " ,i, " of ",snp_n))
		# qui {

		#Odds of the outcome	
		n_ncase<-dat$ncase[i]
		n_total<-dat$ncase[i]+dat$ncontrol[i]
		odds <- n_ncase/(n_total-n_ncase)
		
		maf<-dat$maf[i]
		Z<-dat$z[i]

		# Ns given the MAF

		# p^2 + 2pq + q^2 = 1

		n0 <- n_total*maf^2
		n1 <- n_total*2*maf*(1-maf)
		n2 <- n_total*(1-maf)^2
		N <- n_total
		z <- Z

		#Simulate values of the log-OR, and estimate the Z score

		# preserve
		# clear
		n<-1:1000000
		if(z >= 0) 
		{
			# gen x = _n*0.000001
			x <- n*0.000001
		}else{
			x <- n*-0.000001
		}

		# gen n = _n 

		p0 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)))
		p1 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)-x))
		p2 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)-2*x))

		a <- n0*p0*(1-p0)+n1*p1*(1-p1)+n2*p2*(1-p2)
		b <- n1*p1*(1-p1)+4*n2*p2*(1-p2)
		c <- (n1*p1*(1-p1)+2*n2*p2*(1-p2))^2

		se <- sqrt(a/(a*b-c))
		y = abs(x/se-z)

		y<- abs(x/sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))-z)

		complete<- 0

		j <- 0

		while(complete == 0) 
		{

			# qui su _y
			# Mean<-as.numeric(summary(y)[4])
			# qui su n if _y == r(min)
			Min<-as.numeric(summary(y))[1] #finds the minimum value of _y
			n<-which(y == Min)
			
			# If the minimum isn’t the last observation
			# if r(mean) < 1000000 
			if(n < 1000000)
			{
				complete <- 1
			}else{# If the minimum is the last observation, it hasn’t reached the actual minimum yet, so increase/decrease the observations
				j <- j+1
				if(z >= 0)
				{
					x <- n*0.000001 + (1000000-100)*0.000001*j
				}else{
					x <- -(n*0.000001 + (1000000-100)*0.000001*j)
				}
			}
		 
			# Update the difference between the estimated and observed Z scores

			y <- abs(x/sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))-z)

			

		}

		#While loop complete, so minimum difference found

		# su _y
		# su _x if _y == r(min)
		# local x = r(mean)

		Min<-as.numeric(summary(y)[1])
		x_local<-summary(x[which(y==Min)])[4]
		
		y_se <- sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))

		# su _y
		# su _y_se if _y == r(min)
		Min<-as.numeric(summary(y)[1])
		se_local<-summary(y_se[which(y==Min)])[4]

		# restore
		# End of simulation, bring back the original data and update the results

		# replace _log_or = `x’ in `i’
		# replace _log_or_se = `se’ in `i’

		log_or[[i]] <- x_local
		log_or_se[[i]] <- se_local
		# ID[[i]]<-dat$ID[i]
		# snp[[i]]<-dat$rsid[i]		
	}
	dat$lnor_pred <- as.numeric(unlist(log_or))
	dat$lnor_se_pred <- as.numeric(unlist(log_or_se))
	return(dat)
}


#' Predicted standardised beta
#'
#' Predict the standardised beta using sample sise, Z score and minor allele frequency. Returns the predicted standardised beta, proportion of phenotypic variance explained by the SNP (r2) and F statistic for each SNP
#'
#' @param dat the outcome dataset of interest
#' @param beta the effect size column
#' @param se the standard error column
#' @param eaf the effect allele frequency column
#' @param sample_size the sample size column
#' @param pval name of the p value column
#' @param var the variance of the sample mean of the trait of interest. Default assumed to be 1. 
#'
#' @return data frame with predicted standardised beta, r2 and F stat statistics and estimated standard deviation 
#' @export


predict_beta_sd<-function(dat=NULL,beta="beta",se="se",eaf="eaf",sample_size="ncontrol",var=1,pval="p"){
	z <- dat[,beta]/dat[,se]
	Pos<-which(z==0)
	# if z is zero when calculate from beta and se, infer from the P value
	if("p" %in% names(dat) & sum(Pos) !=0) {		
		z2<-stats::qnorm(dat$p[Pos]/2,lower.tail=FALSE)	
		z[Pos]<-z2
	}	
	
	# Pos<-which(Z1!="Inf")
	# plot(abs(z[Pos]),Z1[Pos])
	if(!"maf" %in% names(dat)){
		maf<-dat[,eaf]
		Pos<-which(maf>0.5)
		maf[Pos]<-1-maf[Pos]
		dat$maf<-maf
	}	
	maf<-dat$maf
	n<-dat[,sample_size]
	dat$beta_sd = z/sqrt(2*maf*(1-maf)*(n+z^2))
 	dat$se_sd<-abs(dat$beta_sd/z)
 	Pos<-which(dat$se_sd=="NaN")
 	if(sum(Pos) != 0 ) warning("se_sd is NaN for some SNPs")
 	k<-1
	var<-1
	dat$r2<-2*dat$beta_sd^2*maf*(1-maf)/var
	dat$F_stat<-dat$r2*(n-1-k)/((1-dat$r2)*k )

	estimated_sd <- dat[,beta] / dat$beta_sd
  	estimated_sd <- estimated_sd[!is.na(estimated_sd)]
  # estimate variance for Y from summary data using method 2
  	dat$sd_est <- stats::median(estimated_sd)
	return(dat)
}

