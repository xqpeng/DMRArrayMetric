#' Calculation of differential methylation values for DMRs
#'
#' This function calculates the differential methylation value of the DMR and the number of probes and CpGs included and the length of the DMR
#'
#' @param DMR the converted DMR that converted by function preinput
#' @param beta the beta matrix with case group and control group(The case group is in the front, the group group is in the back)
#' @param human_ref Human genome reference strand from different cancer tissue or different tissue
#' @param Case_name the colnames of case group
#' @param Control_name the colnames of control group
#' @export
#' @examples
#' beta_file = system.file("extdata","example_beta.RDS",package = 'DMRArrayMetric')
#' example_beta = readRDS(beta_file)
#' case = colnames(example_beta)[1:5]
#' control = colnames(example_beta)[6:10]
#' res_DMR = compute_diff(DMR = convert_DMR,beta = example_beta, human_ref = 'Cancer', Case_name=case, Control_name=control)
#'


compute_diff <-
  function(DMR, beta, human_ref=c('Cancer',"Tissue"), Case_name, Control_name){
    res <- as.data.frame(matrix(NA, 0, 7))
    human_ref <- match.arg(human_ref)
    cat(sprintf("[%s] # Calculate the number of cpg and the difference methylation #\n", Sys.time()))
    if(human_ref == 'Cancer'){
      ref = DMRArrayMetric::RefCancer
    }else{
      ref = DMRArrayMetric::RefTissue
    }
    cat("Calculate the difference methylation value for each probe...\n")
    prb_diff <- apply(beta[,Case_name], 1, mean) - apply(beta[,Control_name], 1, mean)
    cat("Calculate the difference methylation value for each DMR...\n")
    for(i in 1:nrow(DMR)){
      match <- intersect(strsplit(DMR[i,5], split=', ')[[1]], rownames(beta))
      DMR_beta <- beta[match, ]
      probe <- rownames(DMR_beta)
      manchr <- ref[probe,]
      numprbs <- length(probe)
      manchr$diff <- abs(prb_diff[probe])
      manchr <- manchr[order(manchr$CHR,manchr$MAPINFO),]

      chr <- DMR[i,1]
      start <- DMR[i,2]
      end <- DMR[i,3]
      DMRlength<- DMR[i,4]
      molecular <- sum(manchr[,3]*manchr[,5])
      Denominator <- sum(manchr[,4])
      diff <- molecular/Denominator
      res=rbind(res,c(chr,start,end,numprbs,Denominator,DMRlength,diff))
    }
    colnames(res) <- c("chr","Start", "End", "prb","cpg", "length", "Diff")
    res$Start <- as.integer(res$Start)
    res$End <- as.integer(res$End)
    res$prb <- as.integer(res$prb)
    res$cpg <- as.integer(res$cpg)
    res$length <- as.integer(res$length)
    res$Diff <- as.numeric(res$Diff)
    res = res[!is.na(res$Diff),]
    cat(sprintf("[%s] Done\n", Sys.time()))
    res
}



