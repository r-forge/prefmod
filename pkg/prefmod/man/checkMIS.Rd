\name{checkMIS}
\alias{checkMIS}
\title{Function to check/report missing values in paired comparison studies}
\description{
For a given paired comparisons data set the function calculates and prints the
number of missing comparisons and the number of times objects are
missing.  It can also be used to avoid failure of nonresponse-parameter
for nonresponse models in
}
\usage{
checkMIS(obj, nitems, verbose = FALSE)
}
\arguments{
  \item{obj}{dataframe or datafile path/name (like \code{\link{pattPC.fit}}).}
  \item{nitems}{the number of compared objects, not the number of comparisons
       (like \code{\link{pattPC.fit}}).}
  \item{verbose}{if \code{TRUE} printed output, otherwise only invisible output
        to be used, e.g., in the specification of \code{MISalpha} in
        \code{\link{pattPC.fit}}.}
}
\value{
  a logical vector (returned invisibly) specifying for which object
  there are NA responses in the data (\code{obj)}.
}
\seealso{\code{\link{pattPC.fit}}}
\examples{
## no missing NAs in dataset dat4

data(dat4)
checkMIS(dat4,nitems=4,verbose=TRUE)

## generates data set with three items and some missing values
## in comparison (23), column 3, then there are no NAs for
## object 1
data3<-dat4[,1:3]
idx3<-sample(1:100,10)
data3[idx3,3]<-NA
checkMIS(data3,nitems=3,verbose=TRUE)

## estimate MCAR PC pattern model for data3 with NA indicators
## alpha1 cannot be estimated being accomodated by using checkMIS

pattPC.fit(data3, nitems=3, MISalpha=checkMIS(data3,nitems=3))
}
\keyword{models}
\keyword{multivariate}