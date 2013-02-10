library(spdep)

################################################################################
#
# Load some Ireland data set with n=25
#
################################################################################
fn <- system.file("etc/shapes/eire.shp", package="spdep")[1]
prj <- CRS("+proj=utm +zone=30 +units=km")
eire <- readShapeSpatial(fn, ID="names", proj4string=prj)
fn <- system.file("etc/misc/geary_eire.txt", package="spdep")[1]
ge <- read.table(fn, header=TRUE)
row.names(ge) <- as.character(ge$county)
all.equal(row.names(ge), row.names(eire))
eire_ge <- spCbind(eire, ge)
eire_ge1 <- eire_ge[!(row.names(eire_ge) %in% "Dublin"),]
length(row.names(eire_ge1))
nb <- poly2nb(eire_ge1)
nb



################################################################################
#
# Compute Impacts for SAR LM
#
################################################################################
n   <- 25
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
listw <- nb2listw(nb, style="W")
W   <- as_dgRMatrix_listw(listw)
rho <- 0.75
S <- (I_n - rho * W)
beta <- c(-0.5, 0.7)
X <- cbind(1, rnorm(n))
eps <- rnorm(n, sd=0.3)

z   <- as.double(solve(S)%*%(X %*% beta + eps))

# estimate spatial lag model
fit1 <- lagsarlm(z ~ X - 1, listw=listw)
summary(fit1)

impacts(fit1, listw=listw)

#
#Impact measures (lag, exact):
#       Direct  Indirect     Total
#X1 -0.8636143 -1.413445 -2.277059
#X2  0.8880665  1.453465  2.341531


# Theoretische Impacts
# Average Direct Impact =  n^(-1) diag(S^(-1)) beta_r
direct <- mean(diag(solve(S))) * beta
# Average Total Impact =  n^(-1) 1' S^(-1) 1 beta_r
total <- mean(rowSums(solve(S))) * beta
# Average Indirect Impacts =  Average Total Impact - Average Direct Impact
indirect <- total - direct
cbind(direct, indirect, total)

#       direct  indirect total
#X1 -0.6301404 -1.369860  -2.0
#X2  0.8821966  1.917803   2.8

# Impacts mit den Schätzern:
# X1=-0.728322, X2=0.748944
beta_est <- c(-0.728322, 0.748944)
rho_est  <- 0.68015
S_est    <- I_n - rho_est * W

direct_est <- mean(diag(solve(S_est))) * beta_est
total_est <- mean(rowSums(solve(S_est))) * beta_est
indirect_est <- total_est - direct_est
cbind(direct_est, indirect_est, total_est)

#   direct_est indirect_est total_est
#X1 -0.8636156    -1.413458 -2.277074
#X2  0.8880683     1.453479  2.341548

