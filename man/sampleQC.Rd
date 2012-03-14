\name{sampleQC}
\alias{sampleQC}
\alias{sampleQC-method}
\alias{sampleQC,matrix-method}
\alias{sampleQC,AffyBatch-method}
\alias{sampleQC,LumiBatch-method}

\title{

  Sample quality control for FFPE expression data

}

\description{

  Expression data from FFPE tissues may contain a much larger range of
  quality than data from fresh-frozen tissues.  This function sorts
  samples by some specified measure of array quality, Interquartile Range
  (IQR) by default, and plots the correlation of each sample's expression
  values against a ``typical'' sample for the study, as a function of the
  quality measure.  A ``typical'' sample can either be the median
  pseudochip (default), or a specified number of samples with quality
  measure most similar to that sample.  An attempt is made to
  automatically select a threshold for rejection of low-quality samples,
  at the point of largest negative inflection of a Loess smoothing
  curve. for this plot.

}

\usage{

sampleQC(data.obj,logtransform = TRUE, goby = 3, xaxis = "notindex", QCmeasure = "IQR", cor.to = "pseudochip", pseudochip.samples = 1:ncol(data.obj), detectionTh = 0.01, manualcutoff = NULL, mincor = 0, maxcor = 0.8, below.smoothed.threshold = 1.5, lowess.f = 1/3, labelnote = NULL, pch = 1, lw = 4, linecol = "red", make.legend = TRUE, main.title = NA, ...)
       \S4method{sampleQC}{LumiBatch}(data.obj,logtransform = TRUE, goby = 3, xaxis = "notindex", QCmeasure = "IQR", cor.to = "pseudochip", pseudochip.samples = 1:ncol(data.obj), detectionTh = 0.01, manualcutoff = NULL, mincor = 0, maxcor = 0.8, below.smoothed.threshold = 1.5, lowess.f = 1/3, labelnote = NULL, pch = 1, lw = 4, linecol = "red", make.legend = TRUE, main.title = NA, ...)
       \S4method{sampleQC}{AffyBatch}(data.obj,logtransform = TRUE, goby = 3, xaxis = "notindex", QCmeasure = "IQR", cor.to = "pseudochip", pseudochip.samples = 1:ncol(data.obj), detectionTh = 0.01, manualcutoff = NULL, mincor = 0, maxcor = 0.8, below.smoothed.threshold = 1.5, lowess.f = 1/3, labelnote = NULL, pch = 1, lw = 4, linecol = "red", make.legend = TRUE, main.title = NA, ...)
       \S4method{sampleQC}{matrix}(data.obj,logtransform = TRUE, goby = 3, xaxis = "notindex", QCmeasure = "IQR", cor.to = "pseudochip", pseudochip.samples = 1:ncol(data.obj), detectionTh = 0.01, manualcutoff = NULL, mincor = 0, maxcor = 0.8, below.smoothed.threshold = 1.5, lowess.f = 1/3, labelnote = NULL, pch = 1, lw = 4, linecol = "red", make.legend = TRUE, main.title = NA, ...)

}

\arguments{
  
  \item{data.obj}{
    
    A data object of class LumiBatch, AffyBatch, or matrix.  If matrix,
    the columns should contain samples and the rows probes.  If using
    QCmeasure="IQR", it is critical that the data not be normalized.
    QCmeasure="ndetectedprobes" currently works only for LumiBatch
    objects.
    
  }
  
  \item{logtransform}{
    
    If TRUE, data will be log2-transformed before calculating IQR and
    correlation.
    
  }
  \item{goby}{
    
    This number of samples above and below each sample will be used to
    for calculating correlation.  Used only if cor.to="similar".
    
  }
  
  \item{xaxis}{
    If "index", the QC measure will be converted to ranks.  This can be
    useful for very discontinuous values of the QC measure, which
    interfere with generation of a smoothing line.  If "notindex", the
    QC measure is used as-is.
    
  }
  
  \item{QCmeasure}{
    
    Automated options include "IQR" and "ndetectedprobes".  These are
    the Interquartile Range and number of probes called present,
    respectively.  QCmeasure can also be a numeric vector of length
    equal to the number of samples, to manually specify some other
    quality metric.
    
  }

  \item{cor.to}{
    
    "similar" to calculate correlation of each chip to neighbors within
    a sliding window of size 2*goby+1, or "pseudochip" to calculate
    correlation to a study-wide pseudochip.  The former can be more
    sensitive, or more appropriate with a large number of failed chips,
    but does not work well with small sample size (<20).  The default is
    "pseudochip".
    
  }
  
  \item{pseudochip.samples}{
    
    An integer vector, specifying the column numbers of samples to use
    in calculation of the median pseudochip.  Default is to use all
    samples.
    
  }
  
  \item{detectionTh}{
    
    Nnominal detection p-value to consider a probe as detected or not
    (0.01 by default).  Used only if QCmeasure="ndetectedprobes" and
    class(data.obj)=="LumiBatch"
    
  }
  \item{manualcutoff}{
    
    Optional manual specification of a cutoff for good and bad samples.
    If xaxis="index", manualcutoff specifies the number of samples that
    will be rejected.  If xaxis="notindex", it is the value of the QC
    measure plotted on the x-axis, below which samples will be rejected.
    
  }

  \item{mincor}{
    
    Optional specification of a minimum correlation to the
    sliding window samples or median pseudochip, below which all samples
    will be rejected for QC.  This is drawn as a horizontal line on the
    output plot.  
    
  }
  \item{maxcor}{

    Optional specification of an upper limit of correlation to sliding
    window samples or median pseudochip, above which samples will not be
    considered in determining the QC cutoff.  This can be useful if some
    structure in the high-quality end of the plot has the maximum
    downward inflection, causing most samples to be incorrectly
    rejected.  In such case, specifying maxcor can force the otherwise
    automatically-determined threshold into a more reasonable region.
    Only values above at least 0.25, and probably below 0.8, make sense.
    
  }
  \item{below.smoothed.threshold}{
    
    Samples falling more than below.smoothed.threshold times the IQR of
    the residuals will be rejected for low QC.  Large negative residuals
    from the Loess best-fit line may indicate outlier samples even if
    that sample has a high IQR or other quality measure.
    
  }
  
  \item{lowess.f}{
    
    Degree of smoothing of the Loess best-fit line (see ?loess)
    
  }

  \item{labelnote}{
    
    An optional label for the plot, used only if QCmeasure is a numeric vector.
    
  }
  
  \item{pch}{
    
    Plotting character to be used for points (see ?par).
    
  }
  
  \item{lw}{
    
    Line width for Loess curve.
    
  }
  
  \item{linecol}{
    
    Line color for Loess curve.
    
  }
  
  \item{make.legend}{
    
    If TRUE, an automatic legend will be added to the plot.
    
  }
  
  \item{main.title}{
    
    If specified, this over-rides the automatically-generated title.
    
  }
  \item{\dots}{
    
    Other arguments passed on to plot().
    
  }
  
}

\details{
  
  These methods aid in the identification of low-quality samples from
  FFPE expression data, when technical replication is not available.
  
}
\value{

  If only one method is specified (one value each for xaxis, QCmeasure,
  and cor.to, the output is a dataframe with the following columns:

  \item{i}{index or QC measure of each sample}
  \item{spearman}{spearman correlation to sliding window samples or to
    median pseudochip}
  \item{movingaverage}{moving average smoothing of spearman correlation}
  \item{interpolate.i}{evenly spaced QC measure (index or actual QC
    measure) used for plotting Loess curve}
  \item{smoothed}{values of the Loess curve}
  \item{ddy}{second derivative of the Loess curve}
  \item{rejectQC}{was this sample rejected in the QC process?  logical TRUE or FALSE}

  If more than one method is specified, the output is a list, where each
  element contains a dataframe of the above description.
  
}

\references{
  
  Under review.
  
}
\author{
  
  Levi Waldron <lwaldron@hsph.harvard.edu>
  
}


\examples{

library(ffpeExampleData)
data(lumibatch.GSE17565)

QC <- sampleQC(lumibatch.GSE17565,xaxis="index",cor.to="pseudochip",QCmeasure="IQR")

##sort samples
QCvsRNA <- data.frame(inputRNA.ng=lumibatch.GSE17565$inputRNA.ng,rejectQC=QC$rejectQC)
QCvsRNA <- QCvsRNA[order(QCvsRNA$rejectQC,-QCvsRNA$inputRNA.ng),]

##QC rejects samples with lowest input RNA concentration\n
par(mgp=c(4,2,0))
dotchart(log10(QCvsRNA$inputRNA.ng),
         QCvsRNA$rejectQC,
         xlab="log10(RNA conc. in ng)",
         ylab="rejected?",
         col=ifelse(QCvsRNA$rejectQC,"red","black"))

}

\keyword{ hplot }% __ONLY ONE__ keyword per line
