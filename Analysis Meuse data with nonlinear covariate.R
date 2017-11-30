##
## Analysis of Meuse data set
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
resp = log(meuse$zinc)
dist = (meuse$dist)
sqrtdist = sqrt(meuse$dist)
cadmium = meuse$cadmium
application.data = data.frame(cbind(resp,x.coord,y.coord,dist,sqrtdist,cadmium))

plot(x.coord, y.coord, pch=16, cex=resp/4, col=3,
	xlim=c(178,182), ylim=c(329,334))
plot(dist,resp)
plot(cadmium,resp)
plot(sqrtdist,resp)


#--- FIT MODEL TO DATA: WITH DISTANCE TO RIVER ---#
#-------------------------------------------------#
# with direct likelihood estimation
spatial.fit1 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(1,1,0.20), covariogram.model="exponential", trend.input="~sqrtdist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit1) # AIC=159.84

spatial.fit2 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(1,1,0.20), covariogram.model="circular", trend.input="~sqrtdist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit2) # AIC=158.12

spatial.fit3 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", trend.input="~sqrtdist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit3) # AIC=157.44
summary(spatial.fit3)

spatial.fit4 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(1,1,0.20), covariogram.model="matern", trend.input="~sqrtdist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit4) # AIC=158.44

spatial.fit5 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(1,1,0.20), covariogram.model="spherical", trend.input="~sqrtdist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit5) # AIC=158.21


# with direct likelihood doubly-iterative estimation
spatial.fit6 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="exponential", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit6) # AIC=161.71

# with direct likelihood doubly-iterative estimation
spatial.fit7 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="circular", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit7) # AIC=161.47

# with direct likelihood doubly-iterative estimation
spatial.fit8 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit8) # AIC=162.96

# with direct likelihood doubly-iterative estimation
spatial.fit9 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="matern", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit9) # AIC=161.76

spatial.fit10 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="spherical", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit10) # AIC=161.44
summary(spatial.fit10)



spatial.fit11 = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="thin.plate", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
extract.covariogram.parameters(spatial.fit11) # AIC=198.92

# final.fit
ptm = proc.time()
final.fit = fit.spatial.process(data=application.data, method=6, number.of.knots=100,
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", 
	trend.input="~dist",nonlinear.trend.input="~dist",
	phi.upper.limit=5)
proc.time() - ptm
extract.covariogram.parameters(final.fit) # AIC=162.96
summary(final.fit)

sum(final.fit$df.adjusted[3:102])
sum(final.fit$df.adjusted[103:130])



fit.spatial.process(data=application.data, method=6, number.of.knots=155,
	ini.parms=c(1,1,0.20), covariogram.model="gaussian", 
	trend.input="~sqrtdist",nonlinear.trend.input=NULL,
	phi.upper.limit=5)




spatial.fit3$model.components$residuals
par(mfrow=c(1,2))
plot(1:155, predict(final.fit) - application.data[,1])
plot(1:155, -spatial.fit3$model.components$residuals) 



#------- PREDICT KRIGING: WITH DISTANCE TO RIVER -----------#
#-----------------------------------------------------------#
data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame(sqrt(meuse.grid$dist))
cov.pred = cbind(dist.pred)
names(cov.pred) = c("sqrtdist")
pred.loc = data.frame(cbind(x.pred,y.pred))

pred.krig = predict.results(fit=spatial.fit3, data=application.data, model="gaussian", 
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

image.plot(plot.x, plot.y, plot.M.mean, zlim=c(4.4, 7.5),
	col=terrain.colors(300),xlab="x",ylab="y")

image.plot(plot.x, plot.y, plot.M.var, zlim=c(0.0, 0.15),
	col=terrain.colors(300),xlab="x",ylab="y")


#-------- PREDICT SPLINES: WITH DISTANCE TO RIVER ----------#
#-----------------------------------------------------------#
data(meuse.grid)
meuse.grid <- meuse.grid[order(meuse.grid$x, meuse.grid$y), ]
x.pred = meuse.grid$x/1000
y.pred = meuse.grid$y/1000
dist.pred = data.frame((meuse.grid$dist))
cov.pred = cbind(dist.pred)
names(cov.pred) = c("dist")
pred.loc = data.frame(cbind(x.pred,y.pred))
data.boot = application.data[,c("resp","x.coord","y.coord","dist")]


pred.spline = predict.results(fit=final.fit, data=application.data, model="gaussian", 
		trend.input="~dist", nonlinear.trend.input="~dist",
		locations=pred.loc, cov.locations=cov.pred )

boot.var.spline = bootstrap.variance(fit=final.fit, data=data.boot, 
		model="gaussian", method=6, 
		trend.input="~dist", nonlinear.trend.input="~dist",	
		locations=pred.loc, cov.locations=cov.pred,
		boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)

pred.spline.all = cbind(pred.spline, boot.var.spline)
head(pred.spline.all)
str(pred.spline.all)

plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
plot.M.mean2 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
plot.M.var2 <- matrix(nrow=length(plot.x),ncol=length(plot.y))
count <- 1
for (i in 1:length(plot.x)){
for (j in 1:length(plot.y)){
	if (count==3104) stop("it worked")
	if ( (plot.x[i]==meuse.grid$x[count]) &  (plot.y[j]==meuse.grid$y[count]) ){
		plot.M.mean2[i,j] <- pred.spline.all$pred[count]
		plot.M.var2[i,j] <- pred.spline.all[count,5]
		count <- count+1
	}
}}

image.plot(plot.x, plot.y, plot.M.mean2, zlim=c(4.4, 7.5),
	col=terrain.colors(300),xlab="x",ylab="y")

image.plot(plot.x, plot.y, plot.M.var2, zlim=c(0.0, 0.150),
	col=terrain.colors(300),xlab="x",ylab="y")



#--- PLOT NONLINEAR COVARIATE
#----------------------------

# calculate variance
boot.effect.one = bootstrap.effect(final.fit, data.boot, model="gaussian", method=6,
	trend.input="~dist", nonlinear.trend.input="~dist",	
	locations=pred.loc, cov.locations=cov.pred,
	boot.samples=1, approximate=TRUE, seeding.bootstrap=987654)
C.fit.new = boot.effect.one$C.pred
dim(C.fit.new)
Cr = C.fit.new[,2:130]
tilde.C = cbind( C.fit.new [,1] ,  (diag(3103) - (1/3103))%*%Cr)

boot.effect.results = bootstrap.effect(final.fit, data.boot, model="gaussian", method=6,
	trend.input="~dist", nonlinear.trend.input="~dist",	
	locations=pred.loc[1000,], cov.locations=cov.pred[1000,],
	boot.samples=500, approximate=TRUE, seeding.bootstrap=987654)

# calculate mean effect
plot.y.unsorted = final.fit$coef$fixed[2]*tilde.C[,2] + tilde.C[,103:130] %*% unlist(final.fit$coef$rand)[101:128]
plot.dist = sort(cov.pred$dist)
plot.y.sorted = plot.y.unsorted[order(cov.pred$dist)]

all.dist.lines = matrix(NA, nrow=3103, ncol=500 )
for (i in 1:500){
  all.dist.lines[,i] = unlist( boot.effect.results$pars[i,2]*tilde.C[,2] + tilde.C[,103:130] %*% boot.effect.results$pars[i,103:130])
}

var.vec = apply(all.dist.lines, 1, var)
ll = plot.y.sorted-1.96*sqrt(var.vec[order(cov.pred$dist)])
ul = plot.y.sorted+1.96*sqrt(var.vec[order(cov.pred$dist)])

plot(plot.dist, plot.y.sorted,	type="l", lwd=3, ylim=c(-2.0,1.4), 
	xlab="Distance to river", ylab="effect on mean log(cadmium) level")
polygon(c(plot.dist, rev(plot.dist)), 
	c(ul, rev(ll) ),
     	col = "grey80", border = NA)
lines(plot.dist, plot.y.sorted, lwd=3)
points(dist, rep(-2.08,length(dist)), pch="|")
sqrt.trend = as.vector((diag(3103) - (1/3103))%*%sqrt(plot.dist))*spatial.fit3$beta[2]
lines(plot.dist, sqrt.trend, col=2, lwd=3)




# COVARIOGRAM PLOT
#-----------------
x.dist = seq(0.0,1.0,0.01)
cov.krig = 0.102*exp(-(x.dist/0.218)^2)
cov.k.spline = 0.105*exp(-(x.dist/0.189)^2)
plot(x.dist, cov.krig, type="l", lwd=4.5, ylim=c(0,0.20),
	xlab="Distance (km)",ylab="Covariance",
	xaxt="n", yaxt="n", main="", cex.lab=1.8, cex.main=1.8)
points(0, 0.102 + 0.086, pch=16, cex=1.2)
lines(x.dist, cov.k.spline, lwd=4.5, col="gray40")
points(0, 0.105+0.078, pch=16, cex=1.2, col="gray40")
axis(1,at=seq(0.0,1.0,0.2),seq(0.0,1.0,0.2), cex.axis=1.8)
axis(2,at=seq(0,0.2,0.05),seq(0,0.2,0.05), cex.axis=1.8)
legend(0.4,0.20,c("D-ML covariogram","KS-ML: covariogram"), lwd=3, 
	col=c(1,"gray40"), bty="n", cex=1.35)







-----------------------------------

# SPLINE MODELS WITH LESS KNOTS THAN OBSERVATIONS
#------------------------------------------------

knots.vector = c(10,15,20,30,40,50,60,80,100,125,150)
results = matrix(nrow=length(knots.vector), ncol=5)

for (i in 1:length(knots.vector)){
fit = fit.spatial.process(data=application.data, method=6,
	number.of.knots = knots.vector[i], ini.parms=c(1,1,0.20),
 	covariogram.model="gaussian", trend.input="~dist",
	nonlinear.trend.input="~dist",, phi.upper.limit=2)
parms = extract.covariogram.parameters(fit)
results[i,1] = parms[1]
results[i,2] = parms[2]
results[i,3] = parms[3]
results[i,4] = parms[4]
results[i,5] = parms[5]
print(i)
}

fit = fit.spatial.process(data=application.data, method=6,
	number.of.knots = 155, ini.parms=c(1,1,0.20),
 	covariogram.model="gaussian", trend.input="~dist",
	nonlinear.trend.input="~dist",, phi.upper.limit=2)
parms = extract.covariogram.parameters(fit)
results2 = as.matrix(rbind(results, parms))


par(mfrow=c(2,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))
plot( c(knots.vector,155), results2[,1] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",c[s])),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,155), results2[,1], pch=16, cex=1.5)

plot( c(knots.vector,155), results2[,2] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",sigma[epsilon]^2)),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,155), results2[,2], pch=16, cex=1.5)

plot( c(knots.vector,155), results2[,3] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",tau)),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,155), results2[,3], pch=16, cex=1.5)

plot( c(knots.vector,155), results2[,4] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab="AIC value",
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,155), results2[,4], pch=16, cex=1.5)



-----------------------------------------------------------------------------------------------------



data(meuse.riv)
meuse.sr = Polygon(meuse.riv)
x.meuse = meuse.sr@coords[,1]/1000
y.meuse = meuse.sr@coords[,2]/1000

del.plot.M.mean = plot.M.mean
del.plot.M.mean2 = plot.M.mean2
del.plot.M.var = plot.M.var
del.plot.M.var2 = plot.M.var2
del.plot.M.mean[plot.M.var2>0.15]=NA
del.plot.M.mean2[plot.M.var2>0.15]=NA
del.plot.M.var[plot.M.var2>0.15]=NA
del.plot.M.var2[plot.M.var2>0.15]=NA

# 2 x 2 plot
par(mfrow=c(3,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))

	# plot 1
min(resp);quantile(resp,p=0.25);quantile(resp,p=0.50);quantile(resp,p=0.75);max(resp)
q.index = rep(NA,length(resp))
q.index[4.7<=resp & resp<5.3] = 1
q.index[5.3<=resp & resp<5.8] = 2
q.index[5.8<=resp & resp<6.5] = 3
q.index[6.5<=resp & resp<7.6] = 4

plot(x.coord, y.coord, xlim=c(178,182),ylim=c(329.5,334), pch=16,col=1,cex=0.01,
 xlab="Coord X (km)",ylab="Coord Y (km)", xaxt="n", yaxt="n", main="(a)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
points(x.coord[q.index==1], y.coord[q.index==1], pch=16,col="gray60" ,cex=0.75)
points(x.coord[q.index==2], y.coord[q.index==2], pch=16,col="gray40" ,cex=1.00)
points(x.coord[q.index==3], y.coord[q.index==3], pch=16,col="gray20" ,cex=1.50)
points(x.coord[q.index==4], y.coord[q.index==4], pch=16,col=1 ,cex=2.00)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)
legend(178, 334.35, c("[4.7-5.3[","[5.3-5.8[","[5.8-6.5[","[6.5-7.6["), pch=16, 
   pt.cex=c(0.75,1.00,1.50,2.00), col=c("gray60","gray40","gray20",1),
   ncol=1, bty="n", cex=1.7)

	# plot 2
plot(plot.dist, plot.y.sorted, type="l", lwd=3, ylim=c(-2.0,1.4), xaxt="n", yaxt="n",
	xlab="Distance to river (km)", ylab="Effect on mean log(zinc)",
	main="(b)", cex.lab=1.8, cex.main=1.8)
polygon(c(plot.dist, rev(plot.dist)), 
	c(ul, rev(ll) ),
     	col = "grey80", border = NA)
lines(plot.dist, plot.y.sorted, lwd=3)
points(dist, rep(-2.08,length(dist)), pch="|")
sqrt.trend = as.vector((diag(3103) - (1/3103))%*%sqrt(plot.dist))*spatial.fit3$beta[2]
lines(plot.dist, sqrt.trend, col=1, lwd=3, lty=2)
axis(1,at=seq(0,1,0.25),seq(0,1,0.25), cex.axis=1.8)
axis(2,at=seq(-2,1,1),seq(-2,1,1), cex.axis=1.8)
legend(0.30,1.45, c("non-linear effect","square root effect"), col=1, lty=c(1,2), lwd=3,
	 bty="n", cex=1.7)

	# plot 3
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, del.plot.M.mean, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(4.4, 7.5),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(c)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 4
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, del.plot.M.var, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(0.0, 0.15),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(d)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 5
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, del.plot.M.mean2, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(4.4, 7.5),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(e)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)

	# plot 6
plot.x = seq(min(meuse.grid$x),max(meuse.grid$x),40)
plot.y = seq(min(meuse.grid$y),max(meuse.grid$y),40)
image.plot(plot.x/1000, plot.y/1000, del.plot.M.var2, col=heat.colors(200), xlim=c(177.90,182.10),
	ylim=c(329.5,334.10), zlim=c(0.0, 0.15),  xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(f)", cex.lab=1.8, cex.main=1.8)
polygon(c(x.meuse), c(y.meuse), col = "royalblue", border = NA)
axis(1,at=c(178,179,180,181,182),c(178,179,180,181,182), cex.axis=1.8)
axis(2,at=c(330,331,332,333,334),c(330,NA,332,NA,334), cex.axis=1.8)















image.plot(plot.x, plot.y, plot.M.var, zlim=c(0.0, 0.15),
	col=terrain.colors(300),xlab="x",ylab="y")
