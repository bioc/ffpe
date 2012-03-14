\name{CATplot}
\alias{CATplot}

\title{
  
  Make a Concordance at the Top plot.
  
}

\description{

  For the i top-ranked members of each list, concordance is defined as
  length(intersect(vec1[1:i],vec2[1:i]))/i.  This concordance is plotted
  as a function of i.
  
}

\usage{

CATplot(vec1, vec2, maxrank = min(length(vec1),length(vec2)), make.plot = TRUE, ...)

}

\arguments{

  \item{vec1, vec2}{
    Two numeric vectors, for computing concordance.  If these are
    numeric vectors with names, the numeric values will be used for
    sorting and the names will be used for calculating concordance.
    Otherwise, they are assumed to be already-ranked vectors, and the
    values themselves will be used for calculating concordance.
  }

  \item{maxrank}{

    Optionally specify the maximum size of top-ranked items that you
    want to plot.
    
  }
  
  \item{make.plot}{
  
    If TRUE, the plot will be made.  Set to FALSE if you just want the
    concordance calculations.
    
  }
  
  \item{\dots}{

    Optional arguments passed onto plot()
    
  }
  
}

\value{

  Returns a dataframe with two columns:

  \item{i}{length of top lists}
  \item{concordance}{fraction in common that the two provided lists have
    in the top i items}

}

\references{

  The CAT-plot was suggested by Irizarry et al.:
  
  Irizarry, R.A. et al. Multiple-laboratory comparison of microarray
  platforms. Nat Meth 2, 345-350 (2005).

  The CAT-boxplot for multiple splits of a single dataset was suggested
  by Waldron et al. (under review).
  
}

\author{

  Levi Waldron <lwaldron@hsph.harvard.edu>
  
}

\examples{

library(ffpeExampleData)
data(lumibatch.GSE17565)

##preprocessing, individually for rep1 and rep2
lumibatch.rep1 <- lumibatch.GSE17565[,lumibatch.GSE17565$replicate==1]
lumbiatch.rep1 <- lumiT(lumibatch.rep1,"log2")
lumbiatch.rep1 <- lumiN(lumibatch.rep1,"quantile")
probe.var <- apply(exprs(lumibatch.rep1),1,var)
lumibatch.rep1 <- lumibatch.rep1[probe.var > median(probe.var),]

lumibatch.rep2 <- lumibatch.GSE17565[,lumibatch.GSE17565$replicate==2]
lumibatch.rep2 <- lumiT(lumibatch.rep2,"log2")
lumibatch.rep2 <- lumiN(lumibatch.rep2,"quantile")
lumibatch.rep2 <- lumibatch.rep2[featureNames(lumibatch.rep1),]

##row t-tests for differential expression
library(genefilter)
ttests.rep1 <- rowttests(exprs(lumibatch.rep1),fac=factor(lumibatch.rep1$cell.type))
ttests.rep2 <- rowttests(exprs(lumibatch.rep2),fac=factor(lumibatch.rep2$cell.type))

pvals.rep1 <- ttests.rep1$p.value;names(pvals.rep1) <- rownames(ttests.rep1)
pvals.rep2 <- ttests.rep2$p.value;names(pvals.rep2) <- rownames(ttests.rep2)

## Very high concordance between top differentially expressed gene lists
## identified by different replicates
x <- CATplot(pvals.rep1,pvals.rep2,maxrank=1000,xlab="Size of top-ranked gene lists",ylab="Concordance")
legend("topleft",lty=1:2,legend=c("Actual concordance","Concordance expected by chance"), bty="n")

}

\keyword{ hplot }
