#' Calculation metric DMRn for each DMRs
#'
#' This function calculates the metric DMRn
#'
#' @param DMR the DMR is obtained by compute_diff()
#' @export
#' @examples
#' DMRn = compute_metric(DMR = res_DMR)
#'


compute_metric <- function(DMR){

  cat(sprintf("[%s] # Calculate the metric DMRn #\n", Sys.time()))
  interval <- as.data.frame(matrix(NA,10,3))
  rownames(interval) <- c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                          "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
  colnames(interval) <- c("Tprb","Tcpg", "Tlength")
  interval$Tprb<-0
  interval$Tcpg<-0
  interval$Tlength<-0
  t = max(DMR[,7])
  count = 0
  while(t < 1){
    t = t*10
    count = count + 1
  }
  cat("The group differences of DMR range from 0 to", paste0(1 * 0.1^(count-1)),",so the threshold will be",as.character(seq(0,1*0.1^(count-1),0.1^count)[1:10]),"...\n")

  cat("Calculate the total number of cpg, probe and DMR length for each interval...\n")
  for(i in c(1:10)){
    t1 <- (i-1)*(0.1^count)
    t2 <- i*(0.1^count)
    tmp <- DMR[which(DMR[,7] > t1 & DMR[,7] <= t2),]
    interval[i,1] <- sum(tmp[,4])
    interval[i,2] <- sum(tmp[,5])
    interval[i,3] <- sum(tmp[,6])

  }
  cat("Calculate the DMRn at the threshold from", as.character(seq(0,1*0.1^(count-1),0.1^count)[1]),"to",as.character(seq(0,1*0.1^(count-1),0.1^count)[10]),"\n")
  w <- c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)
  res_Qn <- c()
  for(i in c(1:10)){
    interval$weight <- c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)
    tmp <- interval[i:dim(interval)[1], ]
    tmp <- tmp[tmp$Tprb!=0,]
    if(dim(tmp)[1]==0){
      res_Qn <- append(res_Qn,0)
      next;
    }
    Qn_tmp<-sum(tmp$weight*tmp$Tprb^3*tmp$Tcpg/tmp$Tlength)/sum(w)
    Qn_tmp <- log2(Qn_tmp+1)
    res_Qn <- append(res_Qn,Qn_tmp)
  }
  cat(sprintf("[%s] Done\n", Sys.time()))
  names(res_Qn) <- as.character(seq(0,1*0.1^(count-1),0.1^count)[1:10])
  res_Qn
}

