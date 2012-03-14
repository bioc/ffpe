setGeneric("sortedIqrPlot", function(data.obj,
                                     batchvar=rep(1,ncol(data.obj)),
                                     dolog2=FALSE,...)
           standardGeneric("sortedIqrPlot"))

setMethod("sortedIqrPlot",signature(data.obj="matrix"), function(data.obj,
                                      batchvar=rep(1,ncol(data.obj)),
                                      dolog2=FALSE,...)
          {
### data.obj=matrix
            if(!length(batchvar)==ncol(data.obj))
              stop("length(batchvar) must be equal to ncol(data.obj)")
            if(dolog2)
              {
                if(min(data.obj) <= 0)
                  {
                    data.obj <- data.obj - min(data.obj) + 1
                  }
                data.obj <- log2(data.obj)
              }
            central50.range <- apply(data.obj,2,quantile,probs=c(0.25,0.75),na.rm=TRUE)
            central50.range <- central50.range[,order(batchvar,apply(central50.range,2,diff))]
            index <- 1:ncol(central50.range)
            expression.intensity <- seq(min(central50.range),max(central50.range),length.out=ncol(central50.range))
            plot(index,
                 expression.intensity,
                 type='n',
                 ...)
            for (i in 1:ncol(central50.range))
              {
                segments(x0=i,
                         y0=central50.range[1,i],
                         y1=central50.range[2,i])
              }
            persample.iqr <- apply(data.obj,2,IQR,na.rm=TRUE)
            invisible(persample.iqr)
          }
          )

setMethod("sortedIqrPlot",signature(data.obj="LumiBatch"), function(data.obj,
                                      batchvar=rep(1,ncol(data.obj)),
                                      dolog2=FALSE,...)
          {
### data.obj=LumiBatch
            expr.dat <- exprs(data.obj)
            sortedIqrPlot(expr.dat, batchvar, dolog2,...)
          }
          )

setMethod("sortedIqrPlot",signature(data.obj="AffyBatch"), function(data.obj,
                                      batchvar=rep(1,ncol(data.obj)),
                                      dolog2=FALSE,...)
          {
### data.obj=AffyBatch
            expr.dat <- exprs(data.obj)
            sortedIqrPlot(expr.dat, batchvar, dolog2,...)
          }
          )
