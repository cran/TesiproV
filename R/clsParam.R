
#' @title Object for parametric Studies
#' @description Object to create probabilistic problems in parametric studies context. There are no changes how to use compared with SYS_PROB
#'
#' @field beta_params Outputfield: See the beta values of the studie
#' @field res_params Outputfield: See the the full result output of each run
#'
#'
#' @export SYS_PARAM
#' @exportClass SYS_PARAM
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau


SYS_PARAM <- setRefClass(
  Class="SYS_PARAM",
  contains="SYS_PROB",
    fields = list(
      beta_params="list",
      res_params="list",
      nParams="vector",
      nLsfs="numeric",
      nMachines="numeric"
  )
)


SYS_PARAM$methods(
  list(
      runMachines = function(){
        #Check integrity of basisvar´s first
        for(i in 1:length(sys_input)){
          sys_input[[i]]$check()
        }

        res_full_loc <- list()
        p<-0
        paramVar_i <- 0
        runNames <- vector()
        nLsfs <<- length(sys_input)
        nMachines <<- length(probMachines)
        for(m in 1:length(sys_input)){
          for(n in 1:length(sys_input[[m]]$vars)){
            var <- sys_input[[m]]$vars[[n]]
            if(class(var)[1]=="PARAM_BASEVAR"){
              paramVar_i <- paramVar_i + 1
              nParams[paramVar_i] <<- length(var$ParamValues)
              for(u in 1:length(var$ParamValues)){
                 p <- p+1

                runName <- sprintf("Paramrun No. %d. Paramname %s with value(Type): %f (%s)",p, var$Name, var$getCurrentParam(),var$ParamType)
                runNames <- c(runNames,runName)

                beta_loc<-matrix(nrow=length(probMachines),ncol=length(sys_input))
                colnames(beta_loc)<-vector("character",length(sys_input))
                rownames(beta_loc)<-vector("character",length(probMachines))
                #Loop through each Problem in the System
                for(i in 1:length(sys_input)){
                  if(!isempty(sys_input[[i]]$expr)){
                    lsf_expr <- sys_input[[i]]$expr
                  }

                  lsf<-sys_input[[i]]$getLSF()
                  #lsf<-sys_input[[i]]$func
                  if(isempty(sys_input[[i]]$name)){sys_input[[i]]$name<<-"Unknown Problem Name"}
                  colnames(beta_loc)[i]<-sys_input[[i]]$name

                  #Build basic distribution list for probMachines
                  distr<-list()
                  for(k in 1:length(sys_input[[i]]$vars)){
                    v<-sys_input[[i]]$vars[[k]]
                    distr[[k]]<-v$getlDistr()
                  }

                  #Loop through each probMachine
                  res_machine <- list()
                  for (j in 1:length(probMachines)){
                    if(probMachines[[j]]$fCall == "MVFOSM" ){
                      options <- probMachines[[j]]$options
                      if(options$isExpression == "TRUE"){
                        params<-list("lsf"=lsf_expr,"lDistr"=distr)
                      }else{
                        params<-list("lsf"=lsf,"lDistr"=distr)
                      }
                    }else{
                      params<-list("lsf"=lsf,"lDistr"=distr)
                    }

                    if(!isempty(probMachines[[j]]$options)){
                      params<-c(params,probMachines[[j]]$options)
                    }

                    res<-do.call(probMachines[[j]]$fCall,params)
                    rownames(beta_loc)[j]<-probMachines[[j]]$name
                    beta_loc[j,i]<-res$beta
                    res_machine[[j]]<-res
                  }
                  res_full_loc[[i]]<-res_machine
                }

                beta_params[[p]]<<-beta_loc
                res_params[[p]]<<-res_full_loc

                # set the parametric basisvariable one step further
                var$nextParam()

                # check if there are determenistic vars that also needs to be shifted
                for(z in 1:length(sys_input[[m]]$vars)){
                  if(class(sys_input[[m]]$vars[[z]])[1] =="PARAM_DETVAR"){
                    sys_input[[m]]$vars[[z]]$nextParam()
                  }
                }

              } #next parameter iteration

            }
          }
        }
        names(beta_params) <<- runNames
        names(res_params) <<- runNames
      },
      printResults = function(path=""){

        if(!path==""){
          sink(paste(path,".txt",sep = ""),type="output")
          cat("Ergebnisausdruck TesiproV Berechnung\n")
          cat(date())
          cat("\n")
        }

        n_sys <- length(sys_input)
        n_machines <- length(probMachines)
        cat("\n -----------  1. Berechnungsergebnisse (Zusammenfassung)    ---------\n\n")

        cat("\n1.1 Ergebnisse je Loesungsalgorithmus\n")
        for (i in 1:n_machines) {
          cat("\n____________________________________________\n")
          for(j in 1:nParams){
            loc_res <- res_params[[j]][[1]][[i]]
            cat(sprintf("Durchlauf Nr. %.0f \tBeta: %.4f\tPf: %f\tMethod: %s\tRuntime: %s\n",j,loc_res$beta,loc_res$pf,loc_res$method, loc_res$runtime))
          }
        }

        cat("\n\n-------------------------  2. SYSTEM BESCHREIBUNG ------------------\n")

        cat(sprintf("2.1 Das System umfasst eine Gleichung.\n\n"))
        for(i in 1:n_sys){
          cat("\n___________________________________\n")
          cat(sprintf("2.1.%d Gleichung Nr: %d \t-%s\n",i,i,sys_input[[i]]$name))
          if(!isempty(sys_input[[i]]$expr)){
            cat("Folgender symbolischer Ausdruck wurde definiert:\n")
            print(sys_input[[i]]$expr)
          }
          if(!isempty(sys_input[[i]]$func)){
            cat("Folgende Funktion wurde definiert:\n")
            print(sys_input[[i]]$func)
          }

          n_vars <- length(sys_input[[i]]$vars)
          cat(sprintf("\nBasisvariablen (%d gegeben):\n",n_vars))
          for(k in 1:n_vars){
            v<-sys_input[[i]]$vars[[k]]

            if(class(v)[1]=="PARAM_BASEVAR"){
              cat("\n____________________")
              cat("\nParametervariable:\n")
              cat(sprintf("%d.\tName (ID): %s (%d)\tPackage::Verteilungstyp: %s::%s\tMean: %.2f\tSd: %.2f\tCov: %.2f\tx0: %.2f\tVerteilungsparameter: %.5f\t%.5f\n",
                          k,v$Name,v$Id,v$Package,v$DistributionType,v$Mean,v$Sd,v$Cov,v$x0,v$DistributionParameters[1],v$DistributionParameters[2]))
              cat(sprintf("ParamType: %s\n",v$ParamType))
              cat("ParamValues:")
              cat(v$ParamValues)
              cat("\n____________________\n\n")
              }
            else if(class(v)[1]=="PARAM_DETVAR"){
              cat("\n____________________")
              cat("\nDeterministische, von Parametervariable anhaengige Variable \n")
              cat(sprintf("%d.\tName (ID): %s (%d)\n",k,v$Name,v$Id))
              cat("ParamValues:")
              cat(v$ParamValues)
              cat("\n____________________\n\n")
            }
            else
              {
                cat(sprintf("%d.\tName (ID): %s (%d)\tPackage::Verteilungstyp: %s::%s\tMean: %.3f\tSd: %.3f\tCov: %.3f\tx0: %.3f\tVerteilungsparameter: %.5f\t%.5f\n",
                        k,v$Name,v$Id,v$Package,v$DistributionType,v$Mean,v$Sd,v$Cov,v$x0,v$DistributionParameters[1],v$DistributionParameters[2]))
            }
          }
        }
        cat("\n\n ")

        cat("\n ------------------------  3. METHODEN BESCHREIBUNG ------------------\n")
        if(length(n_machines)>0){
          cat(sprintf("\n\n3.1 Die Berechnung wurde von %d Methoden analysiert.\nIm folgenden werden die Eingangsparameter der Methoden beschrieben.",n_sys))
          for (i in 1:n_machines) {
            cat("\n______________________________________\n")
            cat(sprintf("3.%d Machine Name (Type): %s (%s)\n" ,i,probMachines[[i]]$name, probMachines[[i]]$fCall))
            if(!is.null(probMachines[[i]]$options)){
              cat("Die folgenden Optionen waren von den Standartwerten abweichend:\n")
              for (j in 1:length(probMachines[[i]]$options)) {
                cat(sprintf("\t-%s\t=\t%s\n",names(probMachines[[i]]$options[j]),probMachines[[i]]$options[j]))
              }
            }else{
              cat("Ausschliesslich Standartoptionen verwendet!\n")
            }
          }

        }


        cat("--------------------- ENDE -----------------")


        if(!path==""){
          sink()
        }
      }
    )
  )




#' @title System Limit State Functions
#' @description Interface for LSF through PROB_LSF. No changes.
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @export PARAM_LSF
#' @exportClass PARAM_LSF
#'
PARAM_LSF<-setRefClass(
  Class="PARAM_LSF",
  contains="SYS_LSF"
)



#' @title Object for parametric variable
#' @description Object to create parametric basic variables
#'
#' @field ParamValues A vector of values of the parametric studie (e.g. c(1,3,5,7) or seq(1,10,2))
#' @field ParamType A field to determine what should be parametric. Possible is: "Mean", "Sd", "DistributionType"
#'
#'
#' @export PARAM_BASEVAR
#' @exportClass PARAM_BASEVAR
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'

PARAM_BASEVAR <- setRefClass(
  Class="PARAM_BASEVAR",
  contains = "PROB_BASEVAR",
  fields = list(
    ParamValues = "vector",
    ParamType = "character",
    pos="numeric"
  )
)



PARAM_BASEVAR$methods(
  list(
    initialize = function(...){
      initFields(...)
      pos <<- 1
      nextParam()
    },
    nextParam = function(){
      i <- pos
      if(i>=length(ParamValues)){
        i<-length(ParamValues)
        pos <<- 1
      }else{
          pos <<- pos+1
        }

      DistributionParameters <<- c(0,0)

      if(ParamType=="Mean"){
        Mean <<- ParamValues[i]
        if(Mean==0){
          Cov <<- 0
        }else{
          Cov <<- Sd/Mean
        }
      }else if (ParamType=="Sd"){
        Sd <<- ParamValues[i]
        Cov <<- Sd/Mean
      }else if (ParamType =="DistributionType"){
        DistributionType <<- ParamValues[i]
      }else{
        warning("Parametric type not implemented!")
      }

      prepare()
    },
    getCurrentParam = function(){
      if(pos == 1){
        return(ParamValues[length(ParamValues)])
      }else{
        return(ParamValues[pos-1])
      }
    }
  )
)



#' @title Object for parametric deterministic variable
#' @description Object to create parametric deterministic variables
#'
#' @field ParamValues A vector of values. The first element goes with the first run, second element with second run and so on.
#'
#' @export PARAM_DETVAR
#' @exportClass PARAM_DETVAR
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'

PARAM_DETVAR <- setRefClass(
  Class="PARAM_DETVAR",
  contains = "PROB_BASEVAR",
  fields = list(
    ParamValues = "vector",
    pos="numeric"
  )
)



PARAM_DETVAR$methods(
  list(
    initialize = function(...){
      initFields(...)
      DistributionType <<- "norm" #Eigentlich nicht notwendig, da PROB_BASE var bei DistributionTyp = "" auf "norm" setzt
      pos <<- 1
      nextParam()

    },
    nextParam = function(){
      i <- pos
      if(i>=length(ParamValues)){
        i<-length(ParamValues)
        pos <<- 1
      }else{
        pos <<- pos+1
      }

      DistributionParameters <<- c(0,0)
      Mean <<- ParamValues[i]
      Sd <<- Mean/1e7
      if(Mean==0){
        Cov <<- 0
      }else{
        Cov <<- Sd/Mean
      }
      prepare()
    },
    getCurrentParam = function(){
      if(pos == 1){
        return(ParamValues[length(ParamValues)])
      }else{
        return(ParamValues[pos-1])
      }
    }
  )
)




