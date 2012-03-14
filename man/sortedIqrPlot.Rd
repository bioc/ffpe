\name{sortedIqrPlot}
\alias{sortedIqrPlot}
\alias{sortedIqrPlot,matrix-method}
\alias{sortedIqrPlot,AffyBatch-method}
\alias{sortedIqrPlot,LumiBatch-method}

\title{
  
  Minimal box plot of samples, which get ordered optionally by batch,
  then by IQR.
  
}

\description{
  
  This plots only the 25th to 75th percentile of expression intensities
  (Interquartile Range), sorted from smallest to largest IQR.  This
  modification is more readable than a regular boxplot for large sample
  sizes.  An optional batch variable may be specified, so that the sorting
  is done within each batch.
  
}

\usage{

sortedIqrPlot(data.obj, batchvar = rep(1, ncol(data.obj)), dolog2=FALSE, ...)
      \S4method{sortedIqrPlot}{LumiBatch}(data.obj, batchvar = rep(1, ncol(data.obj)), dolog2=FALSE, ...)
      \S4method{sortedIqrPlot}{AffyBatch}(data.obj, batchvar = rep(1, ncol(data.obj)), dolog2=FALSE, ...)
      \S4method{sortedIqrPlot}{matrix}(data.obj, batchvar = rep(1, ncol(data.obj)), dolog2=FALSE, ...)
}

\arguments{
  
  \item{data.obj}{
    An object of class LumiBatch, AffyBatch, AffyBatch, or matrix.  This
    should contain raw, unnormalized, expression intensities.
  }
  
  \item{batchvar}{
    Optional integer batch variable.  If specified, samples will be sorted
    within batches only.  Default is to assume a single batch for all data.
  }

  \item{dolog2}{
    If TRUE, data will be log2-transformed before plotting.  Default is FALSE.
  }

  \item{\dots}{
    Optional arguments passed on to the plot() function.
  }

}

\details{
  
  Function will be generally called for the side-effect of producing a
  plot of sorted IQRs, but if the output is redirected it also produces
  the IQRs.
  
}

\value{

  If the output is redirected, e.g.:

  output <- sortedIqrPlot()

  the function will return the IQR of each sample.
  
}

\references{
  
  Under review.
  
}

\author{
  
  Levi Waldron <lwaldron@hsph.harvard.edu>
  
}



\seealso{
  sampleQC
}
\examples{

library(ffpeExampleData)
data(lumibatch.GSE17565)
sortedIqrPlot(lumibatch.GSE17565,main="GSE17565")

}

\keyword{hplot}

