\name{ffpe-package}
\alias{ffpe-package}
\alias{ffpe}
\docType{package}
\title{
  
  Quality assessment and control for FFPE microarray expression data

}

\description{

  Identify low-quality data using metrics developed for expression data
  derived from Formalin-Fixed, Paraffin-Embedded (FFPE) data.  Also a
  function for making Concordance at the Top plots (CAT-plots).
  
}

\details{
  
  \tabular{ll}{
    Package: \tab ffpe\cr
    Type: \tab Package\cr
    Version: \tab 1.0.0\cr
    Date: \tab 2011-11-17\cr
    License: \tab GPL (>=2)\cr
    LazyLoad: \tab yes\cr
    biocViews: \tab Microarray, GeneExpression, QualityControl, Bioinformatics\cr
  }

  Quality control of FFPE expression data for Illumina and Affymetrix
  microarrays.  The function sampleQC identifies low-quality
  expression data, using IQR or any other surrogate quality measure for
  expression data.  sortedIqrPlot provides a simplified, sorted boxplot of
  raw expression intensities as a quality summary for the experiment,
  suitable for large sample sizes and multiple batches.

}

\author{

  Levi Waldron

  Maintainer: Levi Waldron <lwaldron@hsph.harvard.edu>

}

\references{
  under review
}

\keyword{ package }

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
