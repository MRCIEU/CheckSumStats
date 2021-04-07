
#' Extract SNPs
#'
#' Exract the summary data for the rsids of interest from a target study. This only works on linux/
#'mac operating systems. Will not work on Windows. 
#'
#' @param snplist a list of rsids of interest, either a character vector or path_to_target_file with the list of 
#' rsids 
#' @param path_to_target_file path to the target file This contains the summary data for the trait of #' interest
#' @param path_to_target_file_sep column/field separator. Default assumes that data is tab separated 
#' @param exact_match search for exact matches. Default TRUE 
#' @param Test.gz is the target data a gz file? Default set to FALSE
#' @param Head Does the file have a header ? Default set to TRUE
#' @param fill argument from read.table. logical. If ‘TRUE’ then in case the rows have unequal length, blank fields are implicitly added.  Default is FALSE
#' @param Comment comment to pass to comment.char in read.table. default = "#"
#'
#' @return data frame
#' @export

extract_snps<-function(snplist=NULL,path_to_target_file=NULL,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=FALSE,fill=FALSE,Comment = "#",Head=TRUE){

    if(is.null(snplist)) stop("you must specifiy a list of rsids, either a character vector or a path_to_target_file containing the list of rsids")

	if(length(snplist)==1){
		snplist<-readLines(snplist)
		if(any(duplicated(snplist))) message("duplicate SNPs present in snplist")
	}
	# if(length(snplist)>1){

    f <- file.path(tempdir(), "temp.txt")

# paste(out_dir,"temp.txt",sep="")
    path_to_snpslist<-file.path(tempdir(), "temp.txt")
	utils::write.table(unique(snplist),path_to_snpslist,col.names=F,row.names=F,quote=F)
	snplist<-path_to_snpslist
    # snplist<-paste(out_dir,"temp.txt",sep="")
	# }
	
    path_to_outfile<-file.path(tempdir(), "output.txt")
    path_to_filehead<-file.path(tempdir(), "filehead.txt")
    path_to_outfile_plus_filehead<-file.path(tempdir(), "output_head.txt")
    # out_dir,"output.txt"
	if(Test.gz){
		system(paste("zcat ", path_to_target_file," | fgrep -wf ",snplist," > ", path_to_outfile,sep="")) 
		# system(paste("zgrep -wf ",snplist," ", path_to_target_file," > ",out_dir,"output.txt",sep="")) 
        # out_dir,"path_to_target_filehead.txt"
        system(paste("zcat ",path_to_target_file," | head -1 > ",path_to_filehead,sep=""))
        # out_dir,"path_to_target_filehead.txt "
        # out_dir,"output.txt"
        # ",out_dir,"output_head.txt"
        system(paste("cat",path_to_filehead,path_to_outfile,">", path_to_outfile_plus_filehead,sep=" "))        
	}

	if(exact_match & !Test.gz){
        # system(paste("grep -w ",snplist," ", path_to_target_file," > ",out_dir,"output.txt",sep="")) 
        system(paste("fgrep -wf ",snplist," ", path_to_target_file," > ",path_to_outfile,sep="")) 
        system(paste("head -1 ",path_to_target_file," > ",path_to_filehead,sep=""))
        system(paste("cat",path_to_filehead,path_to_outfile,">",path_to_outfile_plus_filehead,sep=" ")) 
        # system("cat /projects/MRC-IEU/users/ph14916/fatty_acids_summary/output_head.txt")
    }

    if(!exact_match & !Test.gz){
        # system(paste("head",snplist))
        # system(paste("grep 11_61736411",path_to_target_file))
        # system("cat /projects/MRC-IEU/users/ph14916/fatty_acids_summary/temp.txt")
        system(paste("fgrep -f ",snplist," ", path_to_target_file," > ",path_to_outfile,sep="")) 

        system(paste("head -1 ",path_to_target_file," > ",path_to_filehead,sep=""))
        system(paste("cat",path_to_filehead,path_to_outfile,">",path_to_outfile_plus_filehead,sep=" "))
    }
    
  	if(fill){
	    Res<-utils::read.table(path_to_outfile_plus_filehead,sep=path_to_target_file_sep,head=Head,stringsAsFactors=F,fill=T,comment.char = Comment)
	}else{
		Res<-utils::read.table(path_to_outfile_plus_filehead,sep=path_to_target_file_sep,head=Head,stringsAsFactors=F, comment.char = Comment)
	}
    path_to_target_file_name<-unlist(strsplit(path_to_target_file,split="/"))
    if(nrow(Res)!=0) Res$path_to_target_file<-path_to_target_file_name[length(path_to_target_file_name)]

    if(ncol(Res)==2){
    	warning("Had trouble extracting SNPs, possibly because the column separator was correctly specified. Gonna try again using alternative method.")
    	# print(names(Res))
    	Res<-read_plink(path_to_target_file=path_to_outfile_plus_filehead)
    }
    return(Res)
}    

#' Transform betas
#'
#' Transform betas from a linear model to a log odds ratio scale. Assumes betas have been derived from a linear model of case-control status regressed on SNP genotype (additively coded). 
#'
#' @param dat the target dataset
#' rsids 
#' @param effect the column containing the beta. We wish to transform this to a log odds ratio scale
#' @param effect.se standard error for the beta 
#'
#' @return data frame
#' @export


transform_betas<-function(dat=NULL,effect="lnor",effect.se="se"){
    # formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE  
    beta<-dat[,effect]
    se<-dat[,effect.se]
    u<-dat$ncase/(dat$ncase+dat$ncontrol)
    dat[,effect] <- beta / (u * (1 - u))
    dat[,effect.se]<-se / (u * (1 - u))     
    return(dat)
}



read_plink<-function(path_to_target_file=NULL){
    ref<-readLines(path_to_target_file)
    Dat<-NULL
    for(i in 1:length(ref)){
        # print(i)
        A<-unlist(strsplit(ref[i],split=" "))
        A<-A[A!=""]
        Dat[[i]]<-A
    }
    Dat2<-data.frame(do.call(rbind,Dat),stringsAsFactors=F)
    names(Dat2)<-Dat2[1,]
    Dat2<-Dat2[2:nrow(Dat2),]
    return(Dat2)
}

