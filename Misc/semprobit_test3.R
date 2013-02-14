par(mfrow=c(2,2))

plot(density(as.double(yfull)[ind0], to=0))
lines(density(as.double(z)[ind0], to=0), col="red")

plot(density(as.double(yfull)[ind1], from=0))
lines(density(as.double(z)[ind1], from=0), col="red")

