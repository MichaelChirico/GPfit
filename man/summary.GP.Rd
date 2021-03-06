\name{summary.GP}
\alias{summary.GP}
\title{Summary of GP model fit}
\description{Prints the summary of a class \code{GP} object estimated by \code{GP_fit}}
\usage{
\method{summary}{GP}(object, \dots)}
\arguments{
  \item{object}{a class \code{GP} object estimated by \code{GP_fit}}
  \item{\dots}{for compatibility with generic method \code{\link{summary}}}
}
\details{
  prints the summary of the GP object (\code{object}), by calling 
  \code{\link{print.GP}}
}
\author{Blake MacDonald, Hugh Chipman, Pritam Ranjan}

\seealso{
\code{\link{print.GP}} for more description of the output; \cr
\code{\link{GP_fit}} for more information on estimating the model; \cr
\code{\link{summary}} for more description on the \code{summary} function.
}
\examples{
## 1D example
n = 5; d = 1;
computer_simulator <- function(x){
x = 2*x+0.5;
y = sin(10*pi*x)/(2*x) + (x-1)^4;
return(y)
}
set.seed(3);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
summary(GPmodel)


## 2D Example: GoldPrice Function
computer_simulator <- function(x) {
x1=4*x[,1] - 2; x2=4*x[,2] - 2;
t1 = 1 + (x1 + x2 + 1)^2*(19 - 14*x1 + 3*x1^2 - 14*x2 + 
	6*x1*x2 + 3*x2^2);
t2 = 30 + (2*x1 -3*x2)^2*(18 - 32*x1 + 12*x1^2 + 48*x2 - 
	36*x1*x2 + 27*x2^2);
y = t1*t2;
return(y)
}
n = 10; d = 2;
set.seed(1);
library(lhs);
x = maximinLHS(n,d); 
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
summary(GPmodel)
}
