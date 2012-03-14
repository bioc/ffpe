setGeneric("sampleQC", function(data.obj,
                                logtransform=TRUE,
                                goby=3,
                                xaxis="notindex",
                                QCmeasure="IQR",
                                cor.to="pseudochip",
                                pseudochip.samples=1:ncol(data.obj),
                                detectionTh=0.01,
                                manualcutoff=NULL,
                                mincor=0,
                                maxcor=0.8,
                                below.smoothed.threshold=1.5,
                                lowess.f=1/3,
                                labelnote=NULL,
                                pch=1,lw=4,
                                linecol="red",
                                make.legend=TRUE,
                                main.title=NA,
                                ...)
           standardGeneric("sampleQC"))

setMethod("sampleQC",signature(data.obj="matrix"), function(data.obj,
                                 logtransform, goby, xaxis, QCmeasure, cor.to,
                                 pseudochip.samples, detectionTh, manualcutoff, mincor,
                                 maxcor, below.smoothed.threshold, lowess.f, labelnote,
                                 pch,lw, linecol, make.legend, main.title,...){
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  ## This function sorts samples by either interquartile range or number of probes called as detected.
  ## It then computes Spearman correlation within groups of samples which are similar according to this
  ## measure.  Correlations are plotted against either this measure or by rank, with a Loess curve.
  ## The resulting plot can be helpful in deciding how many samples to discard on the basis of IQR or
  ## number of detected probes.
  ##-----------------------------------------------------------------------------------------------
  ## VARIABLE DEFINITIONS
  ##-----------------------------------------------------------------------------------------------
  ##logtransform: if TRUE, values will be log2-transformed before calculating IQR.  If using only ndetectedprobes as a QC measure, this is irrelevant
  ##goby: the number of chips above and below the center chip to include when calculating rank correlations.  Ignored if cor.to="pseudochip".
  ##xaxis: Plot the actual QC measure (IQR or number of detected probes), or relative rank of the chips (1 being lowest)
  ##QCmeasure: IQR - use interquartile range for ranking arrays
  ##           ndetectedprobes - use number of detected probes for ranking arrays
  ##           MAplot.var - use variance of M in MA plot
  ##           c("IQR","ndetectedprobes") - use both
  ##           numeric vector with length equal to the number of samples - use this for ranking arrays
  ##cor.to: "similar" for correlation within a sliding window, "pseudochip" for correlation to the median pseudochip.
  ##labelnote: optional label for the plots, used only if QCmeasure is a numeric vector
  ##detectionTh: nominal detection p-value to consider a probe as detected or not (0.01 by default)
  ##pch: plotting character to be used for points (see ?par)
  ##lw: line width for Loess curve
  ## linecol: line color for Loess curve
  ## ...: other arguments passed on to plot()
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  ##define the plotting function
  ##-----------------------------------------------------------------------------------------------
  ##-----------------------------------------------------------------------------------------------
  if(! (QCmeasure[1] == "IQR" | class(QCmeasure) == "numeric" | class(QCmeasure) == "integer")) stop("Invalid type for QCmeasure")
  makeplots <- function(expr.dat,goby,xaxis,QC.measure,manualcutoff,mincor,below.smoothed.threshold,labelnote="QC",pch,lw,linecol,cor.to,pseudochip.samples,make.legend){
    if(any(is.na(QC.measure)|is.nan(QC.measure))){
      navals <- which(is.na(QC.measure)|is.nan(QC.measure))
      warning(paste(length(navals),"NA and NaN values of QC.measure were removed, along with the associated samples"))
      expr.dat <- expr.dat[,-navals]
      QC.measure <- QC.measure[-navals]
    }
    expr.sort <- expr.dat[,order(QC.measure)]
    pseudochip.samples <- pseudochip.samples[order(QC.measure)]
    if(cor.to=="similar"){
      if(is.na(main.title)|is.null(main.title)){
        thismain <- paste("Rank correlation of groups of ",2*goby+1," arrays\n grouped by ",labelnote,", lowest to highest",sep="")
      }else{
        thismain <- main.title
      }
      print("Calculating correlation within sliding window...")
      cor.vector <- rep(NA,length(QC.measure))
      for (i in (1+goby):(length(QC.measure)-goby)){
        thiscor <- cor(expr.sort[,(i-goby):min((i+goby),ncol(expr.sort))],method="spearman",use="pairwise.complete.obs")
        cor.vector[i] <- median(thiscor[-(1+goby),1+goby],na.rm=TRUE)
      }
      ##For the first few samples, calculate the median correlation to the first goby samples
      thiscor <- cor(expr.sort[,1:(1+goby)],method="spearman",use="pairwise.complete.obs")
      for (i in 1:(1+goby)){
        cor.vector[i] <- median(thiscor[-i,i],na.rm=TRUE)
      }
      names(cor.vector) <- colnames(expr.sort)
      for (i in 1:length(cor.vector)){
        if (is.na(cor.vector[i])){  ## replace NA values
          if(i < (length(cor.vector)-goby)){  #for the low QC values, NA values should be replaced with 0
            cor.vector[i] <- 0
          }else{
            ##otherwise, replace the last few NA values with the mean correlation of 5 arrays at the top of the QC range.
            cor.vector[i] <- mean(cor.vector[(length(cor.vector)-5):length(cor.vector)],na.rm=TRUE) 
          }
        }
      }
    }else if(cor.to=="pseudochip"){
      if(is.na(main.title)|is.null(main.title)){
        thismain <- paste("Rank correlation to median pseudochip \n grouped by ",labelnote,", lowest to highest",sep="")
      }else{
        thismain <- main.title
      }
      print("Calculating Spearman correlation to median pseudochip using pairwise complete observations.")
      print(paste("Using samples:",paste(colnames(expr.sort)[pseudochip.samples],collapse=", ")))
      pseudochip <- apply(expr.sort[,pseudochip.samples],1,median,na.rm=TRUE)
      cor.vector <- cor(expr.sort,pseudochip,method="spearman",use="pairwise.complete.obs")[,1]
    }
    QC.measure.sort <- QC.measure[order(QC.measure)]
    if(xaxis[1]=="index"){
      cormat.2col <- data.frame(i=1:length(cor.vector),
                                spearman=cor.vector,
                                row.names=names(cor.vector))
    }else{
      cormat.2col <- data.frame(i=QC.measure.sort,
                                spearman=cor.vector,
                                row.names=names(cor.vector))
    }
    cormat.2col <- cormat.2col[!is.na(cormat.2col$spearman),]
    cor.vector <- cor.vector[match(rownames(cormat.2col),names(cor.vector))]
    rejectQC <- rep(FALSE,length(cor.vector));names(rejectQC) <- names(cor.vector)
    if(!is.na(lowess.f)&!is.null(lowess.f)){
      movingavg.n <- 2*goby+1
      myavg <- rev(SMA(rev(cormat.2col$spearman),n=movingavg.n))
      if(class(myavg)=="try-error"){
        warning(paste("error in loess fit for",labelnote))
        if(is.null(manualcutoff)) manualcutoff <- 1
      }else{
        was.na <- as.integer(attributes(na.omit(as.numeric(myavg)))$na.action)
        ##Use original values for those that couldn't be smoothed by SMA:
        myavg[was.na] <- cormat.2col[was.na,"spearman"]
        cormat.2col$movingaverage <- myavg
        cormat.loess <- loess(movingaverage~i,data=cormat.2col,span=lowess.f,)
        cormat.2col$interpolate.i <- seq(min(cormat.2col$i),max(cormat.2col$i),length.out=nrow(cormat.2col))
        cormat.2col$smoothed <- predict(cormat.loess,data.frame(i=cormat.2col$interpolate.i))
        cormat.2col$ddy <- D1D2(cormat.2col$interpolate.i,cormat.2col$smoothed,deriv=2,spar.offset=0.6)$D2
        cormat.2col$ddy[1] <- NA
        ##only consider samples with sufficiently low correlation as the cutoff:
        if(is.na(maxcor) | is.null(maxcor)) maxcor <- 1  ##interpret null or na maxcor as "no maximum" or 1.
        cormat.lowspearman <- cormat.2col[cormat.2col$spearman <= max(maxcor,min(cormat.2col$spearman)),]
        if(is.null(manualcutoff)){
          min.position <- which.min(cormat.lowspearman$ddy)
          if(length(min.position) == 0) min.position <- 1
          manualcutoff <- cormat.lowspearman$interpolate.i[min.position]
        }
        diff.from.smoothed <- cormat.2col$spearman-cormat.2col$smoothed
        reject.below.smoothed <- diff.from.smoothed < (-below.smoothed.threshold*IQR(diff.from.smoothed))
      }
      if(xaxis[1]!="index"){
        manualcutoff <- sum(cormat.2col$i < manualcutoff)
      }
      rejectQC[1:manualcutoff] <- TRUE
      if(exists("reject.below.smoothed")) rejectQC[reject.below.smoothed] <- TRUE
    } #end if(!is.na(lowess.f)&!is.null(lowess.f))
    rejectQC[cor.vector<mincor] <- TRUE
    if(identical(all.equal(names(rejectQC),rownames(cormat.2col)),TRUE)){
      cormat.2col$rejectQC <- rejectQC
    }else{
      warning("Something is wrong with the sample names.")
    }
    mycol <- ifelse(rejectQC,"red","black")
    ##correct label note
    if(xaxis[1]=="index"){
      xlab <- paste(labelnote,"rank, 1 is smallest")
    }else{
      xlab <- labelnote
    }
    ##make the plot:
    plot(spearman~i,data=cormat.2col,
         pch=pch,
         col=mycol,
         xlab=xlab,
         ylab="Spearman correlation",
         main=thismain
         )
    if(mincor>0) abline(h=mincor,col="red")
    abline(v=cormat.2col[manualcutoff,"i"],col="red")
    reject.summary <- summary(cormat.2col$rejectQC)
    if(make.legend){
      legend("bottomright",pch=1,lty=-1,legend=c(paste(reject.summary[3],"rejected"),paste(reject.summary[2],"not rejected")),bty='n',col=c("red","black"))
    }
    if("smoothed" %in% colnames(cormat.2col)){
      lines(smoothed~interpolate.i,data=cormat.2col,lw=lw,col="red")
    }
    return(cormat.2col)
  }
  if(logtransform){
    if(max(data.obj) < 25) warning("Log-transforming what looks like already log-transformed data")
    if(min(data.obj,na.rm=TRUE) <= 0){
      data.obj <- data.obj - min(data.obj,na.rm=TRUE) + 1
    }
    data.obj <- log2(data.obj)
  }
  if(identical(class(QCmeasure),"numeric") | identical(class(QCmeasure),"integer")){
    if(is.na(labelnote)|is.null(labelnote)) labelnote <- "custom QC"
    output <-
      makeplots(data.obj,
                goby,
                xaxis[1],
                QC.measure=QCmeasure,
                cor.to=cor.to[1],
                pseudochip.samples=pseudochip.samples,
                manualcutoff,
                mincor,
                below.smoothed.threshold=below.smoothed.threshold,
                labelnote,
                lw=lw,
                pch=pch,
                linecol=linecol,
                make.legend=make.legend)
  }
  if("IQR" %in% QCmeasure){
    Raw.IQR <- apply(data.obj,2,IQR,na.rm=TRUE)
    Raw.IQR[is.na(Raw.IQR)] <- 0
    thismethod <- paste(xaxis[1],cor.to,"IQR",sep="_")
    output <-
      makeplots(data.obj,
                goby,
                xaxis[1],
                QC.measure=Raw.IQR,
                cor.to=cor.to[1],
                pseudochip.samples=pseudochip.samples,
                manualcutoff,
                mincor,
                below.smoothed.threshold=below.smoothed.threshold,
                labelnote="IQR",
                lw=lw,
                pch=pch,
                linecol=linecol,
                make.legend=make.legend)
  }
  ##Sort the output to the original sample order
  mt <- na.omit(match(colnames(data.obj),rownames(output)))
  output <- output[mt,]
  return(output)
}
          )

setMethod("sampleQC",signature(data.obj="LumiBatch"),
          function(data.obj, logtransform, goby, xaxis, QCmeasure,
                   cor.to, pseudochip.samples, detectionTh, manualcutoff,
                   mincor, maxcor, below.smoothed.threshold, lowess.f,
                   labelnote, pch,lw, linecol, make.legend, main.title, ... ){
            ##
            expr.dat <- exprs(data.obj)
            if("ndetectedprobes" %in% QCmeasure){
              QCmeasure <- apply(detection(data.obj),2,function(x) sum(x<detectionTh)/length(x))
              labelnote <- "fraction detected"
            }
            sampleQC(data.obj=expr.dat, logtransform, goby, xaxis,
                     QCmeasure=QCmeasure,
                     cor.to, pseudochip.samples,
                     detectionTh, manualcutoff, mincor, maxcor,
                     below.smoothed.threshold, lowess.f,
                     labelnote=labelnote, pch, lw, linecol,
                     make.legend, ... )
          }
          )


setMethod("sampleQC",signature(data.obj="AffyBatch"),
          function(data.obj, logtransform, goby, xaxis, QCmeasure,
                   cor.to, pseudochip.samples, detectionTh, manualcutoff,
                   mincor, maxcor, below.smoothed.threshold, lowess.f,
                   labelnote, pch,lw, linecol, make.legend, main.title, ... ){
            ##
            expr.dat <- exprs(data.obj)
            if("ndetectedprobes" %in% QCmeasure)
              stop("QCmeasure = \"ndetectedprobes\" is only implemented for LumiBatch objects.")
            sampleQC(data.obj=expr.dat, logtransform, goby, xaxis,
                     QCmeasure, cor.to, pseudochip.samples,
                     detectionTh, manualcutoff, mincor, maxcor,
                     below.smoothed.threshold, lowess.f,
                     labelnote=labelnote, pch, lw, linecol,
                     make.legend, ... )
          }
          )
