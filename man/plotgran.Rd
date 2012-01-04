\encoding{UTF8}
\name{granplot}
\alias{granplot}
\title{
Histogram with a cumulative percentage curve
}
\description{
This function provides a histogram of the grain-size distribution with a cumulative percentage curve
}
\usage{
granplot(x, xc = 1, hist = TRUE, cum = TRUE, main = "", 
			col.cum = "red", col.hist="gray")
}
\arguments{
  \item{x}{
A numeric matrix or data frame
}
  \item{xc}{
Define a column
}
  \item{hist}{
If TRUE, display a histogram; if FALSE, do not display a histogram
}
  \item{cum}{
If TRUE, display a cumulative percentage curve; if FALSE do not display a cumulative percentage curve
}
  \item{main}{
Add a title to the current plot
}
  \item{col.cum}{
Color in which cumulative percentage curve will be drawn
}
\item{col.hist}{
Color in which histogram will be drawn
}
}
\details{
The obtained graph is the most commonly used by Sedimentologists
}
\value{
A histogram with a cumulative percentage curve
}
\author{
Regis K. Gallon (MNHN) \email{regis.gallon@mnhn.fr}, 
Jerome Fournier (CNRS) \email{fournier@mnhn.fr}
}

\examples{
data(granulo)
granplot(granulo,xc=1,hist=TRUE,cum=TRUE,main="Grain-size Distribution",
			col.hist="grey",col.cum="black")
}