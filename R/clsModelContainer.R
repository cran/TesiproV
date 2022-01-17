#' @title System Probabliation Solution Object
#' @description Object to create probabilistic problems. Including Equation, List of Basisvariable, and Solutionmachines
#'
#' @field sys_input List of SYS_LSFs
#' @field sys_type determining serial or parallel system, not implemented yet
#' @field probMachines list of PROB_MACHINES
#' @field res_single grab results after .runMachines()
#' @import gridExtra
#' @import ggplot2
#' @importFrom methods "new"
#'
#' @examples
#' ps <- SYS_PROB(
#' sys_input=list(SYS_LSF(),SYS_LSF()),
#' probMachines=list(PROB_MACHINE()),
#' sys_type="serial")
#' \dontrun{
#' ps$runMachines()
#' ps$beta_sys
#' ps$res_sys
#' ps$printResults("example_1")
#' ps$saveProject(4,"example_1")
#' }
#'
#' @export SYS_PROB
#' @exportClass SYS_PROB
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
SYS_PROB <- setRefClass(
  Class="SYS_PROB",
  fields = list(
    sys_input = "list", #List of SYS_LSFs
    sys_type="character", #determining serial or parallel system, not implemented yet
    probMachines="list", #list of PROB_MACHINES
    res_single="list",
    res_sys="list",
    beta_single="matrix",
    beta_sys="matrix",
    params_sys="list",
    debug.level="numeric"
  )
)


SYS_PROB$methods(
  list(
    runMachines = function(){
      "Starts solving all given problems (sys_input) with all given algorithms (probMachines). After that one can access via $res...1"
      for(i in 1:length(sys_input)){
        sys_input[[i]]$check()
      }

      if(isempty(debug.level)){
        debug.level <<- 0
      }

      beta_single<<-matrix(nrow=length(probMachines),ncol=length(sys_input))
      colnames(beta_single)<<-vector("character",length(sys_input))
      rownames(beta_single)<<-vector("character",length(probMachines))
      #Loop through each Problem in the System
      for(i in 1:length(sys_input)){
        if(!isempty(sys_input[[i]]$expr)){
          lsf_expr <- sys_input[[i]]$expr
        }

        lsf<-sys_input[[i]]$getLSF()
        #lsf<-sys_input[[i]]$func
        if(isempty(sys_input[[i]]$name)){sys_input[[i]]$name<<-"Unknown problem name"}
        colnames(beta_single)[i]<<-sys_input[[i]]$name

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
                params<-list("lsf"=lsf_expr,"lDistr"=distr,"debug.level"=debug.level)
              }else{
                params<-list("lsf"=lsf,"lDistr"=distr,"debug.level"=debug.level)
              }
            }else{
              params<-list("lsf"=lsf,"lDistr"=distr,"debug.level"=debug.level)
            }

            if(!isempty(probMachines[[j]]$options)){
              params<-c(params,probMachines[[j]]$options)
            }

            res<-do.call(probMachines[[j]]$fCall,params)
            rownames(beta_single)[j]<<-probMachines[[j]]$name
            beta_single[j,i]<<-res$beta
            res_machine[[j]]<-res
          }
          res_single[[i]]<<-res_machine
      }
    },
    calculateSystemProbability = function(calcType="simpleBounds",params=list()){
    "Calculates the system probablity if more than one lsf is given and a system_type (serial or parallel) is set.
      If calcType is empty (or simpleBounds), only simpleBounds are applied to further calculation of single soultions.
      If calcType is MCIS, than a Monte Carlo Importance Sampling Method is used (only for parallel systems available).
      If calcType is MCC, than a Crude Monte Carlo Simulation is used.
      If calcType is MCSUS, than the Subset Sampling Algorithm ll be used.
      You can pass arguments to methods via the params field, while the argument has to be a named list (for example check the vignette)."
      params_sys[[length(params_sys)+1]] <<- list("calcType"=calcType,params)

      if(isempty(debug.level)){
        debug.level <<- 0
      }
      # -------------------- Read and set parameters
      if(exists("cov_user", where=params)){
        cov_user <- params$cov_user
      }else{
        cov_user <- 0.05
      }
      if(exists("n_batch", where=params)){
        n_batch <- params$n_batch
      }else{
        n_batch <- 64
      }
      if(exists("n_max", where=params)){
        n_max <- params$n_max
      }else{
        n_max <- 5e6
      }
      if(exists("use_threads", where=params)){
        use_threads <- params$use_threads
      }else{
        use_threads <- 8
      }
      if(exists("beta_l", where=params)){
        beta_l <- params$beta_l
      }else{
        beta_l <- 8
      }

      if(exists("Nsubset", where=params)){
        Nsubset <- params$Nsubset
      }else{
        Nsubset <- 1e5
      }

      if(exists("variance", where=params)){
        variance <- params$variance
      }else{
        variance <- "uniform"
      }

      if(exists("p0", where=params)){
        p0 <- params$p0
      }else{
        p0 <- 0.1
      }

      if(exists("densityType", where=params)){
        densityType <- params$densityType
      }else{
        densityType <- "norm"
      }

      if(exists("dps", where=params)){
        dps <- params$dps
      }else{
        dps <- NULL
      }


      if(calcType=="simpleBounds"){

        min_beta <- vector("numeric",nrow(beta_single))
        max_beta <- vector("numeric",nrow(beta_single))

        if(sys_type=="parallel"){
          for(m in 1:nrow(beta_single)){ # m are the machines used
            min_beta[m] <- max(beta_single[m,])
            max_beta[m] <- -1* qnorm(pnorm(sum(-beta_single[m,])))
          }
        }else if (sys_type=="serial"){

          for(m in 1:nrow(beta_single)){ # m are the machines used
            min_beta[m] <- -1*qnorm(1-prod(pnorm(beta_single[m,])))
            max_beta[m] <- min(beta_single[m,])
          }
        }

        if(nrow(beta_sys)<1){
          beta_sys <<- matrix(nrow=2,ncol=1)
          beta_sys[1,1] <<- min_beta
          beta_sys[2,1] <<- max_beta
          rownames(beta_sys) <<- c("SB min","SB max")
        }else{
          beta_sys <<- rbind(beta_sys,"SB min"=min_beta)
          beta_sys <<- rbind(beta_sys,"SB max"=max_beta)
        }

        res <- list(
          "method"="Simple Bounds",
          "system_mode"=sys_type,
          "single_betas"=beta_single,
          "min_beta"=min_beta,
          "max_beta"=max_beta
        )

        res_sys[[length(res_sys)+1]] <<- res
        info.print("SimpleBounds",debug.level,c("system_mode","beta_min","beta_max"),c(sys_type,min_beta,max_beta))

      }else if(calcType=="MCIS"){
        lsfs <- list()
        distr <- list()
        for(i in 1:length(sys_input)){
          sys_input[[i]]$check()
          lsfs[[i]]<-sys_input[[i]]$getLSF()
          subdistr <- list()
          for(k in 1:length(sys_input[[i]]$vars)){
            v<-sys_input[[i]]$vars[[k]]
            subdistr[[k]]<-v$getlDistr()
          }
          distr[[i]] <- subdistr
        }


        # ----------------- Start Simulation
        res <- TesiproV::MC_IS(
          lsf=lsfs,
          lDistr = distr,
          sys_type=sys_type,
          cov_user = cov_user,
          n_batch = n_batch,
          n_max = n_max,
          use_threads = use_threads,
          beta_l=beta_l,
          densityType=densityType,
          dps = dps,
          debug.level = debug.level
        )
        if(nrow(beta_sys)<1){
          beta_sys <<- matrix(nrow=1,ncol=1)
          beta_sys[1,1] <<- res$beta
          rownames(beta_sys) <<- "MC IS"
        }else{
          beta_sys <<- rbind(beta_sys,"MC IS"=res$beta)
        }
        res_sys[[length(res_sys)+1]] <<- res

      }else if(calcType=="MCC" | calcType =="MCSUS"){

        lsfs <- list()
        distr <- list()
        n <- 1
        n_vars <- vector("numeric",length(sys_input))
        for(i in 1:length(sys_input)){
          sys_input[[i]]$check()
          lsfs[[i]]<-sys_input[[i]]$getLSF()
          n_vars[i] <-length(sys_input[[i]]$vars)
          for(k in 1:n_vars[i]){
            v<-sys_input[[i]]$vars[[k]]
            distr[[n]]<-v$getlDistr()
            n <- n+1
          }
        }

        if(sys_type=="parallel"){
          sys_lsf <- function(x){
            n_lsfs <- length(lsfs)
            res <- vector("numeric",n_lsfs)
            k <- 1
            for(i in 1:n_lsfs){
              res[i] <- lsfs[[i]](x[k:(k+n_vars[i]-1)])
              k <- k+n_vars[i]
            }
            return(max(res))
          }
        }else if(sys_type=="serial"){
          sys_lsf <- function(x){
            n_lsfs <- length(lsfs)
            res <- vector("numeric",n_lsfs)
            k <- 1
            for(i in 1:n_lsfs){
              res[i] <- lsfs[[i]](x[k:(k+n_vars[i]-1)])
              k <- k+n_vars[i]
            }
            return(min(res))
          }
        }

        if(calcType=="MCC"){
          res <- TesiproV::MC_CRUDE(
            lsf=sys_lsf,
            lDistr = distr,
            cov_user = cov_user,
            n_batch = n_batch,
            n_max = n_max,
            use_threads = use_threads,
            debug.level=debug.level
          )
          if(nrow(beta_sys)<1){
            beta_sys <<- matrix(nrow=1,ncol=1)
            beta_sys[1,1] <<- res$beta
            rownames(beta_sys) <<- "MC Crude"
          }else{
            beta_sys <<- rbind(beta_sys,"MC Crude"=res$beta)
          }
          res_sys[[length(res_sys)+1]] <<- res

        }else{
          res <- TesiproV::MC_SubSam(
            lsf=sys_lsf,
            lDistr = distr,
            Nsubset=Nsubset,
            variance=variance,
            p0=p0,
            debug.level = debug.level
          )
          if(nrow(beta_sys)<1){
            beta_sys <<- matrix(nrow=1,ncol=1)
            beta_sys[1,1] <<- res$beta
            rownames(beta_sys) <<- "MC SubSam"
          }else{
            beta_sys <<- rbind(beta_sys,"MC SubSam"=res$beta)
          }
          res_sys[[length(res_sys)+1]] <<- res
        }
      }
    },
    printResults = function(path=""){
    "TesiproV can create a report file with all the necessary data for you. If you provide a path (or filename, without ending) it will store
      the data there, otherwise it will report to the console. Set the path via setwd() or check it via getwd()."
      if(!path==""){
        sink(paste(path,".txt",sep = ""),type="output")
        cat("Ergebnisausdruck TesiproV Berechnung\n")
        cat(date())
        cat("\n")
      }

      n_sys <- length(sys_input)
      n_machines <- length(probMachines)
      cat("\n -----------  1. Berechnungsergebnisse (Zusammenfassung)    ---------\n\n")

      cat("1.1 Versagensindicies:\n")
      if (nrow(beta_single) > 0){
        cat("\nSicherheitsindices bei Einzelversagen:")
        print(beta_single)
      }

      if (nrow(beta_sys) > 0){
        cat("\nSicherheitsindices bei Systemversagen:")
        print(beta_sys)
      }

      cat("\n1.2.1 Ergebnisse je Problemstellung und Algorithmus\n")

      for (i in 1:n_sys) {
        cat("\n____________________________________________\n")
        cat(sprintf("1.2.1.%d Ergebnisse fuer Problem Nr: %d \t-%s\n",i,i,sys_input[[i]]$name))
        for(j in 1:n_machines){
          cat(sprintf("Beta: %.4f\tPf: %f\t",res_single[[i]][[j]]$beta,res_single[[i]][[j]]$pf))
          cat(sprintf("%s (%s)\t" ,probMachines[[j]]$name, probMachines[[j]]$fCall))
          print(res_single[[i]][[j]]$runtime)
        }
      }

      cat("\n1.2.2 Ergebnisse je Systemberechnung und Simulationsmethode\n")
      if (nrow(beta_sys) > 0){
        for(i in 1:length(res_sys)){
        cat("\n____________________________________________\n")
        cat(sprintf("1.2.2.%d Ergebnisse Systemsimulation %d \t %s\n",i,i,res_sys[[i]]$method))
          if(res_sys[[i]]$method =="Simple Bounds"){
            cat(sprintf("Min Beta: %.4f\t Max Beta: %f\t",res_sys[[i]]$min_beta,res_sys[[i]]$max_beta))
          }else{
            cat(sprintf("Beta: %.4f\tPf: %f\t",res_sys[[i]]$beta,res_sys[[i]]$pf))
          }
        }
      }

      cat("\n\n-------------------------  2. SYSTEM BESCHREIBUNG ------------------\n")

      cat(sprintf("2.1 Das System umfasste %d Gleichungen die in %s Beziehung stehen.\n\n",n_sys, sys_type))
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
          #cat(sprintf("Name (ID): %s (%d)\tBeschreibung",v$Name,v$Id))

          cat(sprintf("%d.\tName (ID): %s (%d)\tPackage::Verteilungstyp: %s::%s\tMean: %.3f\tSd: %.3f\tCov: %.3f\tx0: %.3f\tVerteilungsparameter: %.5f\t%.5f\n",
                      k,v$Name,v$Id,v$Package,v$DistributionType,v$Mean,v$Sd,v$Cov,v$x0,v$DistributionParameters[1],v$DistributionParameters[2]))
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

      if(length(params_sys)>0){
      cat(sprintf("\n\n3.2 Parameter der Systemsimulationen."))
        for (i in 1:length(params_sys)) {
          cat("\n______________________________________\n")
          cat(sprintf("3.%d Simulation Method: %s\n" ,i,params_sys[[i]]$calcType))
          # if(!isempty(params_sys[[i]][[2]])){
          #   cat("Die folgenden Optionen waren von den Standartwerten abweichend:\n")
          #   for (j in 1:length(params_sys[[i]][[2]])) {
          #     cat(sprintf("\t-%s\t=\t%s\n",names(params_sys[[i]][[j]]),params_sys[[i]][[j]]))
          #   }
          #   cat("\n")
          # }else{
          #   cat("Ausschliesslich Standartoptionen (Vgl. Helpdatei) verwendet.\n")
          # }
          print(params_sys)
        }
      }



      cat("\n -------------------------  AUSFUEHRLICHE ERGEBNISSE ------------------\n")
      cat("\n4.1 Einzelprobleme\n")
      if(length(n_sys)>0){
        for (i in 1:n_sys) {
          cat("\n_____________________________________\n")
          cat(sprintf("Ergebnisse fuer Problem Nr: %d \t-%s\n",i,sys_input[[i]]$name))
          for(j in 1:n_machines){
            cat(sprintf("Machine Name (Type): %s (%s)\n" ,probMachines[[j]]$name, probMachines[[j]]$fCall))

            n_erg <- length(res_single[[i]][[j]])
            for(k in 1:n_erg)
              if(names(res_single[[i]][[j]][k])=="data"){
              }else if(names(res_single[[i]][[j]][k])=="runtime"){
                cat("\t -runtime:")
                print(res_single[[i]][[j]][k])
              }else{
                cat(sprintf("\t-%s\t=\t%s\n",names(res_single[[i]][[j]][k]),res_single[[i]][[j]][k]))
              }
          }
        }
      }

      if(length(res_sys)>0){
      cat("\n\n4.2 Systemprobleme\n")
        for(i in 1:length(res_sys)){
          cat("\n_____________________________________\n")
          cat(sprintf("Simulationmethode: %s\n" ,res_sys[[i]]$method))

          n_erg <- length(res_sys[[i]])
          for(k in 1:n_erg){
            if(names(res_sys[[i]][k])=="data"){
            }else if(names(res_sys[[i]][k])=="runtime"){
              cat("\t -runtime:")
              print(res_sys[[i]][k])
            }else{
              cat(sprintf("\t-%s\t=\t%s\n",names(res_sys[[i]][k]),res_sys[[i]][k]))
            }
          }
        }
      }


      cat("--------------------- ENDE -----------------")


      if(!path==""){
        sink()
      }
    },
    saveProject = function(level,filename="tesiprov_project"){
      "You can save your calculation project with saveProject().
      There are four different levels of detail to save
      1st Level: Only the beta values
      2nd Level: The result Objects of single or systemcalculation
      3th Level: All The Probablity System Object, including limit state functions, machines and solutions
      4th Level: An image of your entire workspace"
      if(level==1){
        beta <- list("beta_single"=beta_single,
                    "beta_sys"=beta_sys)
        filename_loc <- paste(filename,"_beta.rds",sep = "")
        saveRDS(beta,filename_loc)
        cat(paste("Beta values successfully stored in\n", getwd(),"/",filename_loc,".\nOpen with readRDS().",sep=""))
      }else if(level==2){
        res <- list("res_single"=res_single,
                    "res_sys"=res_sys)
        filename_loc <- paste(filename,"_res.rds",sep = "")
        saveRDS(res,filename_loc)
        cat(paste("Results successfully stored in\n", getwd(),"/",filename_loc,".\nOpen with readRDS().",sep=""))
      }else if(level==3){
        ps <- list("sys_input"=sys_input,
                   "sys_type"=sys_type,
                   "probMachines"=probMachines,
                   "res_single"=res_single,
                   "res_sys"=res_sys,
                   "beta_single"=beta_single,
                   "beta_sys"=beta_sys)
        filename_loc <- paste(filename,"_ps.rds",sep = "")
        saveRDS(ps,filename_loc)
        cat(paste("Prob. Sys. object successfully saved in\n", getwd(),"/",filename_loc,".\nOpen with readRDS().",sep=""))
      }else if(level==4){
        filename_loc <- paste(filename,".RData",sep = "")
        save.image(filename_loc)
        cat(paste("Entire workspace successfully saved in\n", getwd(),"/",filename_loc,".\nOpen with load().",sep=""))
      }else{
        cat(paste("Level needs to be between 1 (less output) and 4 (entire workspace).",sep=""))
      }

    }
    ,
    plotGraph = function(plotType="sim.performance"){
    "not finally implemented. Do not use."
      for (i in 1:length(sys_input)) {
        sysName <- sys_input[[i]]$name
        for(j in 1:length(probMachines)){
          pmName <- probMachines[[j]]$name
          pmfCall <- probMachines[[j]]$fCall
          if(pmfCall=="MC_IS"){

            if(plotType=="sim.performance"){
              df <- res_single[[i]][[j]]$data

              g1 <- ggplot(df,  aes(x=n_sim))+
                geom_line(aes(y=cov))+
                labs(x="Simulationen (-)",
                     y="CoV (-)",
                     title=paste(pmName,"fuer",sysName))+
                theme_bw()

              g2 <- ggplot(df,  aes(x=n_sim))+
                geom_line(aes(y=time))+
                labs(x="Simulationen (-)",
                     y="Zeit (s)")+
                theme_bw()

              ga<- grid.arrange(g1,g2,ncol=1)
              return(ga)

            }else if(plotType=="sim.beta"){
              df <- res_single[[i]][[j]]$data

              g1 <- ggplot(df,  aes(x=n_sim))+
                geom_line(aes(y=-qnorm(pf)),size=1)+
                geom_hline(yintercept = res_single[[i]][[j]]$beta, linetype=2)+
                geom_text(x=0, y=res_single[[i]][[j]]$beta,vjust="bottom",hjust="left",label=paste("Beta =",round(res_single[[i]][[j]]$beta,4)))+
                labs(x="Simulationen (-)",
                     y="Beta (-)",
                     title=paste(pmName,"fuer",sysName))+
                theme_bw()
              return(g1)

            }
          }
        }
      }

      if(plotType=="sim.hv"){
        hv_func <- res_single[[1]]$hv_func
        plots <- list()
        for(i in 1:length(hv_func)){
          rval <- 0.0001
          samplerate <- 1e5
          x.min <- hv_func[[i]]$qfun(rval)
          x.max <- hv_func[[i]]$qfun(1-rval)
          x <- seq(x.min,x.max, length.out=samplerate)
          dhv <- hv_func[[i]]$dfun(x)
          rhv <- hv_func[[i]]$rfun(samplerate)
          df<- data.frame(x=x,d=dhv)
          plots[[i]] <- ggplot(data=df,aes(x=x,y=d))+
            geom_line()+
            labs(x="Groesse",
                 y="Wahrscheinlichkeit (-)",
                 title=paste("Verbunddichte hv_",i,sep=""))+
            theme_bw()
        }
        return(plots)
      }
    }
  )
)

#' @title System Limit State Functions
#' @description Object that represents a limit state function
#' @field expr prepared for expression like SYS_LSF$expr <- expression(f_ck - d_nom)...
#' @field func prepared for objective functions like SYS_LSF$func <- function(x){return(x[1] + x[2])}
#' @field vars needs list of PROB_BASEVAR-Object
#' @field name Can be added for better recognition. Otherwise the problem will be called "Unkown Problem"
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @examples
#' list_of_vars <- list(PROB_BASEVAR(),PROB_BASEVAR())
#' lsf1 <- SYS_LSF(name="my first lsf", vars=list_of_vars)
#' lsf1$func <- function(var1,var2){var1-var2}
#'
#' @export SYS_LSF
#' @exportClass SYS_LSF
#'
SYS_LSF<-setRefClass(
  Class="SYS_LSF",
  fields=list(
    name="character",
    expr="expression",
    func="function",
    vars="list"
  ),
  methods = list(
    ExpressionToFunction = function(){
      "Transforms a valid expression into a objective function. Need the set of Variables with correct spelled names and IDs"
      s<-deparse(expr)
      s<-substr(s,12,nchar(s)-1)
      for (i in 1:length(vars)){
        symb<-vars[[i]]$Name
        vec<-paste(" x[",vars[[i]]$Id,"] ",sep="")
        s<-gsub(symb,vec,s)
      }
      func<<-eval(parse(text = paste('f <- function(x) {' , s , '}', sep='')))
    },
    FunctionToExpression = function(){
      # dont know if this even is possible
    },check = function(){
      "Checks all variables. You dont need to execute this, since the system object will do anyway."
      for (i in 1:length(vars)) {
        vars[[i]]$prepare()
      }
    },getLSF = function(){
      # internal
      lsf.local <- function(x){
        vars.local <- as.list(x)
        vars.names <- vector()
        for(i in 1:length(vars)){
          vars.names <- c(vars.names,vars[[i]]$Name)
        }
        names(vars.local) <- vars.names
        return(do.call(func,vars.local))
      }
      return(lsf.local)
    }
  )
)

#' @title Object to store the distribution model for base vars
#' @description Object to store the distribution model for base vars...
#'
#' @field Id Place in vector of objective functional expression function(x){x[id]}
#' @field Name name like f_ck, used in the limit state function as input name
#' @field Description Used for better understanding of vars
#' @field DistributionType Distributiontypes like "norm", "lnorm", "weibull", "t", "gamma", etc...
#' @field Package The name of the package the Distribution should be taken from (e.g. "evd")
#' @field Mean The Mean Value of this Basisvariable
#' @field Sd The SD Value of this Basisvariable
#' @field Cov The Cov fitting to Mean and Sd.
#' @field x0 Shiftingparameter
#' @field DistributionParameters Inputparameters of the distribution, may be calculated internally
#'
#' @examples
#' var1 <- PROB_BASEVAR(Name="var1", Description="yield strength",
#' DistributionType="norm", Mean=500, Sd=60)
#' var1$prepare()
#'
#' var2 <- PROB_BASEVAR(Name="var2", Description="Load",
#' DistributionType="gumbel",Package="evd",Mean=40, Sd=3)
#' var2$prepare()
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @export PROB_BASEVAR
#' @exportClass PROB_BASEVAR
PROB_BASEVAR <- setRefClass(
  Class="PROB_BASEVAR",
  fields=list(
    Id="numeric",
    Name="character",
    Description="character",
    DistributionType="character",
    Mean="numeric",
    Sd="numeric",
    Cov="numeric",
    x0="numeric",
    Package="character",
    DistributionParameters="vector"
  ),
  methods=list(
    getlDistr = function(){
      prepare()
      if(isempty(Name)){Name<<-"Unknown Basevar Name"}

      if(isempty(Package)){
        Package<<-"stats"
      }

      d <- function(x){
        if(DistributionType == "lnorm"){
          d <- do.call(paste("d",DistributionType,sep=""),envir=loadNamespace(Package), c(list(x),as.list(DistributionParameters)))
        }else{
          d <- do.call(paste("d",DistributionType,sep=""),envir=loadNamespace(Package), c(list(x),as.list(DistributionParameters)))
        }
        return(d)
      }

      p <- function(x){
        p <- do.call(paste("p",DistributionType,sep=""),envir=loadNamespace(Package),  c(list(x),as.list(DistributionParameters)))
        return(p)
      }

      q <- function(x){
        q <- do.call(paste("q",DistributionType,sep=""),envir=loadNamespace(Package), c(list(x),as.list(DistributionParameters)))
        return(q)
      }

      r <- function(x){
        r <- do.call(paste("r",DistributionType,sep=""),envir=loadNamespace(Package), c(list(x),as.list(DistributionParameters)))
        return(r)
      }

      # d<-get(paste("d",DistributionType,sep=""),envir=loadNamespace(Package))
      # p<-get(paste("p",DistributionType,sep=""),envir=loadNamespace(Package))
      # q<-get(paste("q",DistributionType,sep=""),envir=loadNamespace(Package))
      # r<-get(paste("r",DistributionType,sep=""),envir=loadNamespace(Package))

      return(list(list("d"=d,"p"=p,"q"=q,"r"=r,"name"=Name,"mean"=Mean,"Sd"=Sd,"X0"=x0,
                       "Cov"=Cov,"DistributionType"=DistributionType,"DistributionPackage"=Package),
                        DistributionParameters))

    },prepare = function(){
    "Runs the transformations (from mean, sd -> parameters or the other way round) and checks COV, MEAN and SD fitting together.
    If distribution is not available an error ll be thrown."
      if(isempty(x0)){x0<<-0}

      if(is.null(DistributionParameters)){
        DistributionParameters <<- vector("numeric",2)
      }

      # Check if minimal inputs are given
      if(isempty(Mean)&(isempty(Sd)|isempty(Cov))&((DistributionParameters[1]==0)&(DistributionParameters[2]==0))){
        warning(sprintf("No mean defined for the PROB_BASEVAR Object: %s with ID: %s!\nPlease define at least a mean and CoV or SD!",Name, as.character(Id)))
      }
      # Check if distribution and package is given and known
      if(isempty(DistributionType)){DistributionType <<- "norm"}
      if(isempty(Package)){Package <<- "stats"}
      knownDistributions <- c("norm","lnorm","gumbel", "gamma", "exp", "t", "binom", "unif","weibull","lt")
      knownPackages <- c("stats", "evd")
      if(!(DistributionType %in% knownDistributions & Package %in% knownPackages)){
        warning(sprintf("\tID:%s, Name: %s.
                        The given DistributionType: \"%s\" from \"%s\" package is not known.
                        Thats possible for new or \"exoctic\" Distributions. Please double check.",DistributionType, Package,DistributionType,Package))
      }


      # Transforming happens here
      if(!(isempty(Cov)&isempty(Mean))){
        if(isempty(Sd)){
          #Wenn Cov und Mean gegeben, aber SD leer ist, errechne SD
          Sd <<- Mean*Cov
        }else if (isempty(Cov)){

          if(Mean==0){
            Cov <<- 0
          }else{
            Cov <<- Sd/Mean
          }

        }else{

          if((Mean*Cov - Sd)>0.01){
            warning(sprintf("ID:%s, Name: %s.
                             The given Mean (%f), Coefficent of Variance (%f) and Standard Deviation (%f) doesnt fit together.
                             Please check!",Id, Name,Mean,Cov,Sd))
          }
        }
      }

      if((DistributionParameters[1]==0)&(DistributionParameters[2]==0)){
        if(DistributionType == "norm"){

          DistributionParameters[1]<<-Mean
          DistributionParameters[2]<<-Sd

        }else if(DistributionType=="lnorm"){
          sn <- sqrt(log(1+(Sd/(Mean))^2))
          mn <- log((Mean)/sqrt(1+(Sd/(Mean))^2))

          DistributionParameters[1]<<-mn
          DistributionParameters[2]<<-sn

        }else if(DistributionType == "gumbel"){
          #https://en.wikipedia.org/wiki/Gumbel_distribution
          # digamma(1) bringt negativen Wert, daher wird -1*digamma() gerechnet
          scale <- (Sd*sqrt(6))/pi
          location <- Mean+(digamma(1)*1/scale)

          DistributionParameters[1]<<-location
          DistributionParameters[2]<<-1/scale

        }else if(DistributionType == "gamma"){
          # https://en.wikipedia.org/wiki/Gamma_distribution, See right sided table, column: Method of Moments
          shape <- Mean^2/Sd^2
          scale <- (Sd^2)/Mean

          DistributionParameters[1]<<-shape
          DistributionParameters[2]<<-1/scale

        }else if(DistributionType == "weibull"){
          #https://de.wikipedia.org/wiki/Weibull-Verteilung
          #k is shapeparam, lambda is so called scale param
          #https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
          k <- (Sd/Mean)^(-1.086)
          lambda <- Mean/(gamma(1+1/k))

          DistributionParameters[1]<<-k
          DistributionParameters[2]<<-lambda

        }else{
          warning(sprintf("\tID: %s, Name: %s
          Transformation from Mean and COV into DistributionParameters for %s is not implemented yet!
          Please paste the Parameters into DistributionParameters Field of PROB_BASEVAR Object!
          First Element of vector will be parsed as first argument for Distributionfunction",Id, Name, DistributionType))
        }
      }else if((isempty(Mean))&(isempty(Sd))&!((DistributionParameters[1]==0)&(DistributionParameters[2]==0))){
        # Transformation from Parameters to Mean and Sd
        if(DistributionType == "norm"){

          Mean <<- DistributionParameters[1]
          Sd <<- DistributionParameters[2]
          Cov <<- Sd/Mean
        }else if(DistributionType=="lnorm"){
          mn <- DistributionParameters[1]
          sn <- DistributionParameters[2]

          Mean <<- exp(mn + (sn^2)/2)
          Sd <<- exp(mn + (sn^2)/2)*sqrt(exp((sn^2))-1)
          Cov <<- Sd/Mean
        }else if(DistributionType == "lt"){

          m <- DistributionParameters[1]
          s <- DistributionParameters[2]
          n <- DistributionParameters[3]
          nue <- DistributionParameters[4]

          s_dach <- s/sqrt(n/(n+1))
          loc_sd <- sqrt((s_dach^2)*nue/(nue-2))

          Mean <<- exp(m + (loc_sd^2)/2)
          Sd <<- exp(m + (loc_sd^2)/2)*sqrt(exp((loc_sd^2))-1)
        }else if(DistributionType == "gumbel"){
          #https://en.wikipedia.org/wiki/Gumbel_distribution
          # digamma(1) bringt negativen Wert, daher wird -1*digamma() gerechnet
          location <- DistributionParameters[1]
          scale <- DistributionParameters[2]
          Mean <<- location-digamma(1)*scale
          Sd <<- pi/sqrt(6)*scale

        }else if(DistributionType == "gamma"){
          # https://en.wikipedia.org/wiki/Gamma_distribution, See right sided table, column: Method of Moments
          shape <- DistributionParameters[1]
          scale <- DistributionParameters[2]

          Mean <<- shape*scale
          Sd <<- sqrt(shape * scale^2)
          Cov <<- Sd/Mean
        }else if(DistributionType == "weibull"){
          #https://de.wikipedia.org/wiki/Weibull-Verteilung
          #k is shapeparam, lambda is so called scale param
          #https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
          warning(sprintf("Transformationf or %s from Parameters to Mean & SD for weibull not included",Name))

        }else{
          warning(sprintf("\tID: %s, Name: %s
          Transformation from Mean and COV into DistributionParameters for %s is not implemented yet!
          Please paste the Parameters into DistributionParameters Field of PROB_BASEVAR Object!
          First Element of vector will be parsed as first argument for Distributionfunction",Id, Name, DistributionType))
        }
      }
    }
  )
)


#' @title Object to store a deterministic model for base vars
#' @description Object to store a deterministic model for base vars
#'
#' @field Id Place in vector of objective functional expression function(x){x[id]}
#' @field Name readable name like f_ck, used for transform expression to objective function
#' @field Description - Used for better understanding of vars
#' @field Value - The deterministic value that sould be used (as mean for the normal distribution with infinite small sd)
#'
#' @examples
#' form_rf<-PROB_MACHINE(name="FORM RF",fCall="FORM",options=list("n_optim"=20,
#' "loctol"=0.001, "optim_type"="rackfies"))
#' sorm <- PROB_MACHINE(name="SORM",fCall="SORM")
#' mcis<-PROB_MACHINE(name="MC IS",fCall="MC_IS",options=list("cov_user" = 0.05, "n_max"=300000))
#' mcsus<-PROB_MACHINE(name="MC SuS",fCall="MC_SubSam")
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @export PROB_DETVAR
#' @exportClass PROB_DETVAR
PROB_DETVAR <- setRefClass(
  Class="PROB_DETVAR",
  contains = "PROB_BASEVAR",
  fields = c(
    Value="numeric"
    ),
  methods = list(
    initialize = function(...){
      initFields(...)
      Mean <<- Value
      Sd <<- Mean/1e7
      prepare()
      }
    )
  )



#' @title Object to store prob machines
#'
#' @field name individual name
#' @field fCall Function Call of the method. Possible is: "MVFOSM","FORM", "SORM", "MC_Crude", "MC_IS", "MC_SubSam"
#' @field options additional options for the method provided as a list. For form e.g. options=list("optim_type"="rackfies").
#' To get insight of all available settings of each method open the help with ?FORM, ?SORM, ?MC_IS etc.
#'
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @export PROB_MACHINE
#' @exportClass PROB_MACHINE
PROB_MACHINE<-setRefClass(
  Class="PROB_MACHINE",
  fields = list(
    name="character",
    fCall="character",
    options="vector"
  ),
  methods=list(
  getMethodLevel = function(){
    if(grepl("MC",fCall)){
      return(3)
    }else if(grepl("ORM",fCall)){
      return(2)
    }else if(strcmp("MVFOSM",fCall)){
      return(2)
    }else{
      return(1)
    }
  }
  )
)
