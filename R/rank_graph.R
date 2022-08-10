#' Draw a line chart for DMRn.
#'
#' This function will draw a line chart for DMRn.
#'
#' @param DMRfilePath The text file containing the metric DMRn file path and the name of DMR detection method.
#' @param threshold The methylation difference threshold t, default 0-0.9.
#' @param Title The title of rank graph, default NULL.
#' @export
#' @import ggplot2
#' @examples
#' library("ggplot2")
#' DMRn_file = read.table(system.file('extdata','rank_DMRn.txt',package='DMRArrayMetric'))
#' for(i in 1:7){
#'    DMRn_file[i,1] = eval(parse(text=DMRn_file[i,1]))
#'
#' }
#' rank_graph(DMRn_file,threshold=seq(0,0.9,0.1),Title='GroupA-GroupB')
rank_graph <- function(DMRfilePath, threshold = seq(0,0.9,0.1), Title=NULL){
  cat(sprintf("[%s] # Draw a line chart for DMRn. #\n", Sys.time()))
  if(ncol(DMRfilePath) > 2){
    message("The text file must only have two column, DMRn file path and the name of DMR detection method.\n")
  }
  group = table(DMRfilePath[,2])
  if(sum(group>1)){
    message("A method appears more than once in the file, please ensure that each method appears only once.\n")
  }
  cat("There are", length(group),'methods in the file...\n')
  DMR_data = as.data.frame(matrix(0,0,3))
  if(is.null(Title)){
    Title='Title'
    cat("The title is NULL...\n")
  }
  for(i in 1:nrow(DMRfilePath)){
    file_path = DMRfilePath[i,1]
    DMR_name = DMRfilePath[i,2]
    DMRn = read.table(file_path)
    DMR_data = rbind(DMR_data,cbind(DMRn=as.numeric(DMRn),Group=DMR_name,threshold=threshold))

  }
  DMR_data$DMRn = as.numeric(DMR_data$DMRn)
  p = ggplot(data = DMR_data, aes(x = threshold, y = DMRn)) + geom_line(aes(group =
                                                                          Group, color = Group)) + geom_point(aes(color = Group, shape = Group)) +
    scale_shape_manual(values = 0:nrow(DMRfilePath)) +
    theme(panel.background = element_blank(), axis.line = element_line(), axis.text.y=element_text()) +
    ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5)) +
    xlab("methylation difference threshold t") + ylab(expression(paste('log'[2],"(",italic('DMRn'),")")))
  cat("The plot is stored in the path:",getwd(),"as pdf file...\n")
  cat(sprintf("[%s] Done\n", Sys.time()))
  ggsave(p,file='rank_graph.pdf')
}

