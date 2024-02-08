\name{Katrina}
\alias{Katrina}
\alias{Katrina.raw}
\docType{data}
\title{New Orleans business recovery in the aftermath of Hurricane Katrina}
\description{
This dataset has been used in the LeSage et al. (2011) paper entitled 
"New Orleans business recovery in the aftermath of Hurricane Katrina"
to study the decisions of shop owners to reopen business after Hurricane Katrina. 
The dataset contains 673 observations on 3 streets in New Orleans 
and can be used to estimate the spatial probit models 
and to replicate the findings in the paper.
}
\usage{
 data(Katrina)
}

\format{
  \code{Katrina.raw} is a data frame with 673 observations on the following 15 variables.
  \describe{
    \item{\code{code}}{a numeric vector}
    \item{\code{long}}{longitude coordinate of store}
    \item{\code{lat}}{latitude coordinate of store}
    \item{\code{street1}}{a numeric vector}
    \item{\code{medinc}}{median income}
    \item{\code{perinc}}{a numeric vector}
    \item{\code{elevation}}{a numeric vector}
    \item{\code{flood}}{flood depth (measured in feet)}
    \item{\code{owntype}}{type of store ownership: "sole proprietorship" vs. "local chain" vs. "national chain"}
    \item{\code{sesstatus}}{socio-economic status of clientele (1-5): 1-2 = low status customers, 3 = middle, 4-5 = high status customers}
    \item{\code{sizeemp}}{"small size" vs. "medium size" vs. "large size" firms}
    \item{\code{openstatus1}}{a numeric vector}
    \item{\code{openstatus2}}{a numeric vector}
    \item{\code{days}}{days to reopen business}
    \item{\code{street}}{1=Magazine Street, 2=Carrollton Avenue, 3=St. Claude Avenue}
  }
  
  \code{Katrina} is a data frame with 673 observations on the following 13 variables.
  \describe{
    \item{\code{long}}{longitude coordinate of store}
    \item{\code{lat}}{latitude coordinate of store}
    \item{\code{flood_depth}}{flood depth (measured in feet)}
    \item{\code{log_medinc}}{log median income}
    \item{\code{small_size}}{binary variable for "small size" firms}
    \item{\code{large_size}}{binary variable for "large size" firms}
    \item{\code{low_status_customers}}{binary variable for low socio-economic status of clientele}
    \item{\code{high_status_customers}}{binary variable for high socio-economic status of clientele}
    \item{\code{owntype_sole_proprietor}}{a binary variable indicating "sole proprietor" ownership type}
    \item{\code{owntype_national_chain}}{a binary variable indicating "national_chain" ownership type}
    \item{\code{y1}}{reopening status in the very short period 0-3 months; 1=reopened, 0=not reopened}
    \item{\code{y2}}{reopening status in the period 0-6 months; 1=reopened, 0=not reopened}
    \item{\code{y3}}{reopening status in the period 0-12 months; 1=reopened, 0=not reopened}
  }
}

\details{
The \code{Katrina.raw} dataset contains the data found on the website
before some of the variables are recoded. For example, the 
socio-economic status of clientele is coded as 1-5 in the raw data, 
but only 3 levels will be used in estimation: 
1-2 = low status customers, 3 = middle, 4-5 = high status customers. Hence,
with "middle" as the reference category, 
\code{Katrina} contains 2 dummy variables for low status customers 
and high status customers.

The dataset \code{Katrina} is the result of these recoding operations and can be
directly used for model estimation.
}

\note{
When definining the reopening status variables \code{y1} (0-3 months), \code{y2} (0-6 months),
and \code{y3} (0-12 months) from the \code{days} variable, the Matlab code ignores the seven cases
where \code{days=90}. To be consistent with the number of cases in the paper,
we define \code{y1},\code{y2},\code{y3} in the same way: 
\code{y1=sum(days < 90)}, \code{y2=sum(days < 180 & days != 90)}, 
\code{y3=sum(days < 365 & days != 90)}.
So this is not a bug, its a feature.
}

\source{
The raw data was obtained from the Journal of the Royal Statistical Society dataset website
(was: \code{https://rss.onlinelibrary.wiley.com/pb-assets/hub-assets/rss/Datasets/})
and brought to RData format.
}

\references{
J. P. LeSage, R. K. Pace, N. Lam, R. Campanella and X. Liu  (2011),
New Orleans business recovery in the aftermath of Hurricane Katrina
\emph{Journal of the Royal Statistical Society A}, \bold{174}, 1007--1027
}

\examples{
data(Katrina)
attach(Katrina)
table(y1) # 300 of the 673 firms reopened during 0-3 months horizon, p.1016
table(y2) # 425 of the 673 firms reopened during 0-6 months horizon, p.1016
table(y3) # 478 of the 673 firms reopened during 0-12 months horizon, p.1016
detach(Katrina)

\dontrun{
# plot observations in New Orleans map; Google requires an API key; see `ggmap::register_google()`
if (require(ggmap)) {
  qmplot(long, lat, data = Katrina, maptype="roadmap", source="google")
}
}

\donttest{
# replicate LeSage et al. (2011), Table 3, p.1017
require(spatialreg)
 
# (a) 0-3 months time horizon
# LeSage et al. (2011) use k=11 nearest neighbors in this case
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11))
listw <- nb2listw(nb, style="W")
W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")

# Note: cannot replicate (a) 0-3 months time horizon model as of February 2024  
#fit1 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
#  low_status_customers +  high_status_customers + 
#  owntype_sole_proprietor + owntype_national_chain, 
#  W=W1, data=Katrina, ndraw=600, burn.in = 100, showProgress=FALSE)
#summary(fit1)
  
# (b) 0-6 months time horizon
# LeSage et al. (2011) use k=15 nearest neighbors
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=15))
listw <- nb2listw(nb, style="W")
W2 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")

fit2 <- sarprobit(y2 ~ flood_depth + log_medinc + small_size + large_size +
  low_status_customers + high_status_customers + 
  owntype_sole_proprietor + owntype_national_chain, 
  W=W2, data=Katrina, ndraw=600, burn.in = 100, showProgress=FALSE)
summary(fit2)  

# (c) 0-12 months time horizon
# LeSage et al. (2011) use k=15 nearest neighbors as in 0-6 months
W3 <- W2
fit3 <- sarprobit(y3 ~ flood_depth + log_medinc + small_size + large_size +
  low_status_customers + high_status_customers + 
  owntype_sole_proprietor + owntype_national_chain, 
  W=W3, data=Katrina, ndraw=600, burn.in = 100, showProgress=FALSE)
summary(fit3)

# replicate LeSage et al. (2011), Table 4, p.1018
# SAR probit model effects estimates for the 0-3-month time horizon
# impacts(fit1)  

# replicate LeSage et al. (2011), Table 5, p.1019
# SAR probit model effects estimates for the 0-6-month time horizon
impacts(fit2)

# replicate LeSage et al. (2011), Table 6, p.1020
# SAR probit model effects estimates for the 0-12-month time horizon
impacts(fit3)
}
}
\keyword{datasets}
