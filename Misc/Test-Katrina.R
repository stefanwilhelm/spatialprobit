library(spatialprobit)
library(McSpatial)
library(ggmap)
require(spdep)

# load data set
data(Katrina)

# plot data set
qmplot(long, lat, data = Katrina)
qmplot(long, lat, data = Katrina, maptype="roadmap", source="google")

# (a) 0-3 months time horizon
# LeSage et al. (2011) use k=11 nearest neighbors in this case
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11))
listw <- nb2listw(nb, style="W")
W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")

fit1 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size +
  low_status_customers +  high_status_customers +
  owntype_sole_proprietor + owntype_national_chain,
  W=W1, data=Katrina, ndraw=600, burn.in = 100, showProgress=TRUE)
summary(fit1)

# compute marginal effects
# LeSage et al. (2011), Table 4, p.1018
Rprof("mfx.out")
mfx <- marginal.effects(fit1)
Rprof(NULL)
summaryRprof("mfx.out")
# dauert nur noch 5 Sekunden! Yes!

# print impacts
impacts(fit1)

# (a) Direct effects
cbind(
 lower_005=apply(mfx$direct, 2, quantile, prob=0.05),
 posterior_mean=colMeans(mfx$direct),
 upper_095=apply(mfx$direct, 2, quantile, prob=0.95)
)

# (b) Indirect effects
cbind(
 lower_005=apply(mfx$indirect, 2, quantile, prob=0.05),
 posterior_mean=colMeans(mfx$indirect),
 upper_095=apply(mfx$indirect, 2, quantile, prob=0.95)
)

#                          lower_005 posterior_mean    upper_095
#flood_depth             -0.04742760    -0.02876699 -0.013035431
#log_medinc               0.04637806     0.12977771  0.232806634
#small_size              -0.10524148    -0.04753350 -0.003754389
#large_size              -0.17945257    -0.06090693  0.038052944
#low_status_customers    -0.13313983    -0.06175697 -0.011374714
#high_status_customers   -0.02443048     0.01520672  0.065258738
#owntype_sole_proprietor  0.02439708     0.10198152  0.208338850
#owntype_national_chain  -0.10537807     0.01815482  0.161831983

# (c) Total effects
cbind(
 lower_005=apply(mfx$total, 2, quantile, prob=0.05),
 posterior_mean=colMeans(mfx$total),
 upper_095=apply(mfx$total, 2, quantile, prob=0.95)
)

#                          lower_005 posterior_mean   upper_095
#flood_depth             -0.10554640    -0.07565948 -0.04711513
#log_medinc               0.15784622     0.34112048  0.52018421
#small_size              -0.23801467    -0.12264812 -0.01381990
#large_size              -0.41901052    -0.15464251  0.10096781
#low_status_customers    -0.28680190    -0.15916203 -0.03426724
#high_status_customers   -0.06530528     0.03849546  0.14675022
#owntype_sole_proprietor  0.08618174     0.26302559  0.43712570
#owntype_national_chain  -0.25576245     0.04428228  0.35195825





