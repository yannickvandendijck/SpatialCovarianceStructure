##
## Analysis of Meuse data set with clustering
##
## Author: Yannick Vandendijck
## Paper: Estimating the spatial covariance structure using 
## the geoadditive model. Environ Ecol Stat (2017) 24:341–361.
## Last update: 20/01/2017
##


rm(list=ls())
source(".../Rsource_code.R")
read.in.libraries()


#--- READ IN DATA ---#
#--------------------#
data(meuse)
str(meuse)
x.coord = meuse$x/1000
y.coord = meuse$y/1000
resp = log(meuse$lead)
dist = (meuse$dist)
sqrtdist = sqrt(meuse$dist)
cadmium = meuse$cadmium
cluster.var = as.numeric(meuse$landuse)[!is.na(meuse$landuse)]
application.data = data.frame(cbind(resp,x.coord,y.coord,dist,sqrtdist,cadmium))

plot(x.coord, y.coord, pch=16, cex=resp/4, col=3,
	xlim=c(178,182), ylim=c(329,334))
plot(dist,resp)
plot(cadmium,resp)
plot(sqrtdist,resp)

hist(meuse$zinc)
hist(meuse$cadmium)
hist(meuse$copper)
hist(meuse$lead)


#--- PLOT LANDUSE
#----------------
fit.data = application.data[!is.na(meuse$landuse),]
plot.landuse = meuse$landuse[!is.na(meuse$landuse)]
levels(plot.landuse)

xx = fit.data$x.coord[plot.landuse==(levels(plot.landuse)[1])]
yy = fit.data$y.coord[plot.landuse==(levels(plot.landuse)[1])]
plot(xx, yy, pch=1, xlim=c(178,182), ylim=c(329,334))
for (i in 2:15){
  xx = fit.data$x.coord[plot.landuse==(levels(plot.landuse)[i])]
  yy = fit.data$y.coord[plot.landuse==(levels(plot.landuse)[i])]
  points(xx, yy, pch=i)
}


#--- FIT MODEL TO DATA: WITH DISTANCE TO RIVER ---#
#-------------------------------------------------#
length(cluster.var)
fit.data = application.data[!is.na(meuse$landuse),]

# ignoring landuse
spatial.fit1 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit1) 

spatial.fit2 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="matern", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit2) 

spatial.fit3 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="spherical", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit3) 

ptm = proc.time()
spatial.fit4 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="circular", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
proc.time() - ptm
extract.covariogram.parameters(spatial.fit4) 

spatial.fit5 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="exponential", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit5) 


# including landuse as random effect
spatial.fit6 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5, cluster.input=cluster.var)
extract.covariogram.parameters(spatial.fit6) 			

spatial.fit7 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="matern", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5, cluster.input=cluster.var)
extract.covariogram.parameters(spatial.fit7) 		

spatial.fit8 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="spherical", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5, cluster.input=cluster.var)
extract.covariogram.parameters(spatial.fit8) 			

ptm = proc.time()
spatial.fit9 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="circular", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5, cluster.input=cluster.var)
proc.time() - ptm
extract.covariogram.parameters(spatial.fit9)
		

spatial.fit10 = fit.spatial.process(data=fit.data, method=6, 
	ini.parms=c(1,1,0.20), covariogram.model="exponential", number.of.knots=154,
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5, cluster.input=cluster.var)
extract.covariogram.parameters(spatial.fit10)			




length(spatial.fit4$df.adjusted)
sum(spatial.fit4$df.adjusted[3:156])
sum(spatial.fit4$df.adjusted[157:184])

length(spatial.fit9$df.adjusted)
sum(spatial.fit9$df.adjusted[3:156])
sum(spatial.fit9$df.adjusted[157:184])


#--- SPATIAL COMPONENT WITHOUT USING LANDUSE
#-------------------------------------------------#

extract.covariogram.parameters(spatial.fit4)

data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame( rep(meuse.grid$dist[1000],3103))
cov.pred = cbind(dist.pred)
names(cov.pred) = c("dist")
pred.loc = data.frame(cbind(x.pred,y.pred))
pred1 = spatial.effect(fit=spatial.fit4, data=fit.data, model="circular", 
		locations=pred.loc )

plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.mean1 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.mean1[i,j] <- pred1[count]
		count <- count+1
	}
}}
image.plot(plot.x, plot.y, plot.M.mean1, col=terrain.colors(300),xlab="x",ylab="y")


	### variance
input.data = fit.data[,c("resp","x.coord","y.coord","dist")]

boot.var.spatial1 = c()
for (i in 1:32){
	if (i!=32) pred.loc2 = pred.loc[ ((i-1)*100+1):(i*100), ]
	if (i==32) pred.loc2 = pred.loc[ ((i-1)*100+1):nrow(pred.loc), ]
	if (i!=32) cov.pred2 = cov.pred[ ((i-1)*100+1):(i*100), ]
	if (i==32) cov.pred2 = cov.pred[ ((i-1)*100+1):nrow(cov.pred), ]

	int.result = bootstrap.variance.spatial.effect(spatial.fit4, input.data, model="circular", method=6,
		trend.input="~dist", nonlinear.trend.input="~dist",	
		locations=pred.loc2, cov.locations=cov.pred2,
		boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)
	boot.var.spatial1 = c(boot.var.spatial1,int.result)
print(i)
}
pred.spatial1.all = cbind(pred1, boot.var.spatial1)
head(pred.spatial1.all)


plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.var1 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.var1[i,j] <- pred.spatial1.all[count,2]
		count <- count+1
	}
}}
image.plot(plot.x, plot.y, plot.M.var1, col=terrain.colors(300),xlab="x",ylab="y")



#--- PREDICT USING LANDUSE
#-------------------------------------------------#
extract.covariogram.parameters(spatial.fit9)

data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame( rep(meuse.grid$dist[1000],3103))
cov.pred = cbind(dist.pred)
names(cov.pred) = c("sqrtdist")
cluster.pred = rep(7,nrow(meuse.grid))
pred.loc = data.frame(cbind(x.pred,y.pred))
pred2 = spatial.effect(fit=spatial.fit9, data=fit.data, model="circular", 
		locations=pred.loc )

plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.mean2 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.mean2[i,j] <- pred2[count]
		count <- count+1
	}
}}
image.plot(plot.x, plot.y, plot.M.mean2, col=terrain.colors(300),xlab="x",ylab="y")


	### variance
input.data = fit.data[,c("resp","x.coord","y.coord","dist")]

boot.var.spatial2 = c()
for (i in 1:32){
	if (i!=32) pred.loc2 = pred.loc[ ((i-1)*100+1):(i*100), ]
	if (i==32) pred.loc2 = pred.loc[ ((i-1)*100+1):nrow(pred.loc), ]
	if (i!=32) cov.pred2 = cov.pred[ ((i-1)*100+1):(i*100), ]
	if (i==32) cov.pred2 = cov.pred[ ((i-1)*100+1):nrow(cov.pred), ]
	if (i!=32) cluster.pred2 = cluster.pred[ ((i-1)*100+1):(i*100) ]
	if (i==32) cluster.pred2 = cluster.pred[ ((i-1)*100+1):length(cluster.pred) ]

	int.result = bootstrap.variance.spatial.effect(spatial.fit9, input.data, model="circular", method=6,
		trend.input="~dist", nonlinear.trend.input="~dist",	
		locations=pred.loc2, cov.locations=cov.pred2,
		cluster.input = cluster.var, cluster.locations= cluster.pred2,
		boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)
	boot.var.spatial2 = c(boot.var.spatial2,int.result)
print(i)
}
pred.spatial2.all = cbind(pred2, boot.var.spatial2)
head(pred.spatial2.all)


plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.var2 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.var2[i,j] <- pred.spatial2.all[count,2]
		count <- count+1
	}
}}
image.plot(plot.x, plot.y, plot.M.var2, col=terrain.colors(300),xlab="x",ylab="y")




#--- PLOT NONLINEAR COVARIATE
#----------------------------
data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame((meuse.grid$dist))
cov.pred = cbind(dist.pred)
names(cov.pred) = c("dist")
pred.loc = data.frame(cbind(x.pred,y.pred))
cluster.pred = rep(1,nrow(meuse.grid))
data.boot = fit.data[,c("resp","x.coord","y.coord","dist")]

# calculate variance
one.result = bootstrap.effect(spatial.fit9, data.boot, model="circular", method=6,
	trend.input="~dist", nonlinear.trend.input="~dist",	
	locations=pred.loc, cov.locations=cov.pred,
	cluster.input = cluster.var, cluster.locations= cluster.pred,
	boot.samples=1, approximate=TRUE, seeding.bootstrap=987654)

C.fit.new = one.result$C.pred
dim(C.fit.new)
Cr = C.fit.new[,2:199]
tilde.C = cbind( C.fit.new [,1] ,  (diag(3103) - (1/3103))%*%Cr)

boot.effect.results = bootstrap.effect(spatial.fit9, data.boot, model="circular", method=6,
	trend.input="~dist", nonlinear.trend.input="~dist",	
	locations=pred.loc[1000,], cov.locations=cov.pred[1000,],
	cluster.input = cluster.var, cluster.locations= cluster.pred[1000],
	boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)

# calculate mean effect
plot.y.unsorted = spatial.fit9$coef$fixed[2]*tilde.C[,2] + tilde.C[,157:189] %*% unlist(spatial.fit9$coef$rand)[155:187]
plot.dist = sort(cov.pred$dist)
plot.y.sorted = plot.y.unsorted[order(cov.pred$dist)]

all.dist.lines = matrix(NA,nrow=3103, ncol=500 )
for (i in 1:500){
  all.dist.lines[,i] = unlist( boot.effect.results$pars[i,2]*tilde.C[,2] + tilde.C[,157:184] %*% boot.effect.results$pars[i,155:182])
}

var.vec = apply(all.dist.lines, 1, var)
ll = plot.y.sorted-1.96*sqrt(var.vec[order(cov.pred$dist)])
ul = plot.y.sorted+1.96*sqrt(var.vec[order(cov.pred$dist)])

par(oma=c(0,1,0,1),mar=c(6,6,2,0.5))
plot(plot.dist, plot.y.sorted,	type="l", lwd=3, ylim=c(-2.0,1.4), xaxt="n", yaxt="n", 
	xlab="Distance to river (km)", ylab="Effect on mean log(lead)", cex.lab=1.8)
axis(1,at=seq(0,1,0.25), seq(0,1,0.25), cex.axis=1.8)
axis(2,at=seq(-2, 1, 1), seq(-2, 1, 1), cex.axis=1.8)
polygon(c(plot.dist, rev(plot.dist)), 
	c(ul, rev(ll) ),
     	col = "grey80", border = NA)
lines(plot.dist, plot.y.sorted, lwd=3)
points(dist, rep(-2.08,length(dist)), pch="|")

plot.y.unsorted2 = spatial.fit4$coef$fixed[2]*tilde.C[,2] + tilde.C[,157:184] %*% unlist(spatial.fit4$coef$rand)[155:182]
plot.y.sorted2 = plot.y.unsorted2[order(cov.pred$dist)]
lines(plot.dist, plot.y.sorted2, lwd=3, col=2, lty=2)
legend(0.205,1.25, c("no clustering","with clustering"), col=2:1, lwd=3, lty=2:1, bty="n", cex=2)


# COVARIOGRAM PLOT
#-----------------
parms1 = extract.covariogram.parameters(spatial.fit4)
parms2 = extract.covariogram.parameters(spatial.fit9)

x.dist = seq(0.0,1.0,0.01)
theta = unlist( lapply(x.dist, function(t){min(1,t/parms1[3])}) )
cov.without = parms1[1]*(1-2*( theta*sqrt(1-theta^2) +  asin(theta)) / pi)
theta = unlist( lapply(x.dist, function(t){min(1,t/parms2[3])}) )
cov.with = parms2[1]*(1-2*( theta*sqrt(1-theta^2) +  asin(theta)) / pi)

plot(x.dist, cov.without, type="l", lwd=4.5, ylim=c(0,0.25),
	xlab="Distance (km)",ylab="Covariance",
	xaxt="n", yaxt="n", main="", cex.lab=1.8, cex.main=1.8)
points(0, parms1[1] + parms1[2], pch=16, cex=1.2)
lines(x.dist, cov.with, lwd=4.5, col="gray40")
points(0, parms2[1] + parms2[2], pch=16, cex=1.2, col="gray40")
axis(1,at=seq(0.0,1.0,0.2),seq(0.0,1.0,0.2), cex.axis=1.8)
axis(2,at=seq(0,0.25,0.05),seq(0,0.25,0.05), cex.axis=1.8)
legend(0.4,0.20,c("No clustering","With clustering"), lwd=3, 
	col=c(1,"gray40"), bty="n", cex=1.75)


dev.off()

----------------------------------------------------------------------------------


data(meuse.riv)
meuse.sr = Polygon(meuse.riv)
x.meuse = meuse.sr@coords[,1]/1000
y.meuse = meuse.sr@coords[,2]/1000

# 2 x 2 plot
par(mfrow=c(3,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))

	# plot 1
min(resp);quantile(resp,p=0.25);quantile(resp,p=0.50);quantile(resp,p=0.75);max(resp)
q.index = rep(NA,length(resp))
q.index[3.6<=resp & resp<4.3] = 1
q.index[4.3<=resp & resp<4.8] = 2
q.index[4.8<=resp & resp<5.3] = 3
q.index[5.3<=resp & resp<6.5] = 4

plot(x.coord, y.coord, xlim=c(178,182),ylim=c(329.5,334), pch=16,col=1,cex=0.01,
 xlab="Coord X (km)",ylab="Coord Y (km)", xaxt="n", yaxt="n", main="(a)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
points(x.coord[q.index==1], y.coord[q.index==1], pch=16,col="gray60" ,cex=0.75)
points(x.coord[q.index==2], y.coord[q.index==2], pch=16,col="gray40" ,cex=1.00)
points(x.coord[q.index==3], y.coord[q.index==3], pch=16,col="gray20" ,cex=1.50)
points(x.coord[q.index==4], y.coord[q.index==4], pch=16,col=1 ,cex=2.00)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)
legend(178, 334.35, c("[3.6-4.3[","[4.3-4.8[","[4.8-2.3[","[5.3-6.5["), pch=16, 
   pt.cex=c(0.75,1.00,1.50,2.00), col=c("gray60","gray40","gray20",1),
   ncol=1, bty="n", cex=1.7)

	# plot 2
plot(x.coord, y.coord, xlim=c(178,182),ylim=c(329.5,334), pch=16, cex=0.01, col="white",
 xlab="Coord X (km)",ylab="Coord Y (km)", xaxt="n", yaxt="n", main="(b)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
for (i in 1:15){
  xx = fit.data$x.coord[plot.landuse==(levels(plot.landuse)[i])]
  yy = fit.data$y.coord[plot.landuse==(levels(plot.landuse)[i])]
  points(xx, yy, pch=i, col=1)
}
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 3
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, plot.M.mean1, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(-0.75, 0.75),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(c)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 4
plot(x.coord, y.coord, xlim=c(178,182),ylim=c(329.5,334), pch=16, cex=0.01, col="white",
 xlab="",ylab="", xaxt="n", yaxt="n", main="", cex.lab=1.8, cex.main=1.8)
legend(178.5,334, levels(plot.landuse), pch=1:16, col=1, cex=1.7, bty="n", ncol=2, x.intersp=2, text.width=1)

	# plot 5
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, plot.M.mean2, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(-0.75, 0.75),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="Estimated spatial component", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 6
plot.M.diff = plot.M.mean2-plot.M.mean1
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, plot.M.diff, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(-0.20, 0.10),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(e)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)
i=7
  xx = fit.data$x.coord[plot.landuse==(levels(plot.landuse)[i])]
  yy = fit.data$y.coord[plot.landuse==(levels(plot.landuse)[i])]
  points(xx, yy, pch=i, col=1)



unlist(spatial.fit9$coef$random)[188:202]
fit.data[cluster.var==7,]




#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

data(meuse.riv)
meuse.sr = Polygon(meuse.riv)
x.meuse = meuse.sr@coords[,1]/1000
y.meuse = meuse.sr@coords[,2]/1000
par(mfrow=c(1,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))

	# plot 3
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, plot.M.var1, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(0.0, 0.20),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="no clustering", cex.lab=1.8, cex.main=1.2)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 5
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, plot.M.var2, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(0.0, 0.20),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="Estimated variance spatial component", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------


spatial.kriging = fit.spatial.process(data=application.data, method=1, 
	ini.parms=c(0.10,0.80,0.10), covariogram.model="circular",
	trend.input="~sqrtdist", phi.upper.limit=5)
extract.covariogram.parameters(spatial.kriging) 

data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame( rep(sqrt(meuse.grid$dist)[2500], 3103)  )
cov.pred = cbind(dist.pred)
names(cov.pred) = c("sqrtdist")
pred.loc = data.frame(cbind(x.pred,y.pred))

pred.krig = predict.results(fit=spatial.kriging, data=application.data, model="circular", 
		trend.input="~sqrtdist", locations=pred.loc, cov.locations=cov.pred )

plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.mean <- matrix(nrow=length(plot.x),ncol=length(plot.y))
plot.M.var <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.mean[i,j] <- pred.krig$pred[count]
		plot.M.var[i,j] <- pred.krig[count,4]
		count <- count+1
	}
}}

image.plot(plot.x, plot.y, plot.M.mean-5.6454, 
	col=heat.colors(200),xlab="x",ylab="y")

image.plot(plot.x, plot.y, plot.M.var,
	col=terrain.colors(300),xlab="x",ylab="y")

















data.boot[107,]

bootstrap.variance(spatial.fit9, data.boot, model="circular", method=6,
	trend.input="~dist", nonlinear.trend.input="~dist",	
	locations=data.boot[107,c("x.coord","y.coord")], cov.locations=data.boot[107,c("dist")],
	cluster.input = cluster.var, cluster.locations= cluster.var[107],
	boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)






