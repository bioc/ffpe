CATplot <-
function(vec1,vec2,maxrank=min(length(vec1),length(vec2)),make.plot=TRUE,...){ 
  ##----------------------------------------------------------------------
  ##----------------------------------------------------------------------
  ##Fuction for Concordance at the Top plots
  ##----------------------------------------------------------------------
  ## vec1 and vec2: ranked lists to use for CAT plot
  ## maxrank: only go up to maxrank
  ## make.plot: if TRUE, the plot will be made
  ## ...: arguments passed on to plot
  ##----------------------------------------------------------------------  
  if(class(vec1)=="numeric" & class(vec2)=="numeric" & !is.null(names(vec1)) & !is.null(names(vec1)))
    {
      vec1 <- sort(vec1)
      vec1 <- names(vec1)
      vec2 <- sort(vec2)
      vec2 <- names(vec2)
    }
  if(is.na(maxrank) | is.null(maxrank) | maxrank > min(length(vec1),length(vec2)))
    {
      maxrank <- min(length(vec1),length(vec2))
    }
  output <- data.frame(rank=1:maxrank,
                       concordance=NA)
  for (i in 1:nrow(output)){
    output[i,"concordance"] <- length(intersect(vec1[1:i],vec2[1:i]))/i
  }
  if(make.plot){
    plot(concordance~rank,data=output,type='l',...)
    abline(a=0,b=1/min(length(vec1),length(vec2)),lty=2)
  }
  return(output)
}

