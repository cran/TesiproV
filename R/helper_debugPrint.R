#'
#' @title internal  Helper function to debug more easy
#' @param infoLevel If 0, no Output (just Errors), if 1 little output, if 2 bigger output
#'
#' @param flag Parse additonal info
#' @param values If you check variables then post this into values
#' @param msg here add some extra msg
#' @param type Type can be "INFO" or "ERROR"
#'
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'



debug.print <- function(infoLevel,flag="",values,msg="",type="INFO"){
  t<-format(Sys.time(),"%H:%M:%OS3")

  if (infoLevel==1 && type=="INFO"){
    warning("#",type," - ",t,"\t",msg,values,"\n")
  }else if(infoLevel>1){
    warning("#",type," FROM ", flag, " - ",t,"\t",msg,values,"\n")
  }else if(infoLevel==0 && type=="ERROR"){
    warning("#!",type, flag," - ",t,"\t",values,"!\n\n")
  }

}

debug.print.getDebugLevel <- function(){
  return(0)
}


info.print <- function(tag,verbose,varnames, values){
  if(verbose==2){
    cat("\r#",tag,": ",paste(varnames, ":", values, "\t"))
  }
}
