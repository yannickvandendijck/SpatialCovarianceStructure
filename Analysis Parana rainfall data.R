##
## Analysis of Parana rainfall data set	
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
data(parana)
plot(parana, bor=borders)

x.coord = parana$coords[,1]
y.coord = parana$coords[,2]
resp = parana$data
application.data = data.frame(cbind(resp,x.coord,y.coord))

plot(x.coord,y.coord, xlim=c(100,800),ylim=c(0,550), pch=16,col=5 ,cex=1*(resp/300),
xlab="Coord X (km)",ylab="Coord Y (km)", xaxt="n")
polygon(parana$borders, lwd=5)
points(x.coord,y.coord, pch=16,col=5 ,cex=1*(resp/200))
axis(1,at=seq(100,800,100),seq(100,800,100))


#--- FIT MODEL TO DATA ---#
#-------------------------#
# Diggle & Ribeiro (2002) - kriging - exponential
spatial.fit1 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(100,100,100), covariogram.model="exponential", trend.input="1st",
	phi.upper.limit=5000)
extract.covariogram.parameters(spatial.fit1)

# kriging - gaussian
spatial.fit1.1 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(800,100,100), covariogram.model="gaussian", trend.input="1st",
	phi.upper.limit=5000)
extract.covariogram.parameters(spatial.fit1.1)

# kriging - spherical
spatial.fit1.2 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(100,100,100), covariogram.model="spherical", trend.input="1st",
	phi.upper.limit=5000)
extract.covariogram.parameters(spatial.fit1.2)

# kriging - circular
spatial.fit1.3 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(100,100,100), covariogram.model="circular", trend.input="1st",
	phi.upper.limit=5000)
extract.covariogram.parameters(spatial.fit1.3)

# kriging - matern
spatial.fit1.4 = fit.spatial.process(data=application.data, method=1,
	ini.parms=c(100,100,100), covariogram.model="matern", trend.input="1st",
	phi.upper.limit=5000)
extract.covariogram.parameters(spatial.fit1.4)



# Kammann and Wand (2003) - exponential
spatial.fit2 = fit.spatial.process(data=application.data, method=4, number.of.knots=143,
 	covariogram.model="exponential", trend.input="1st")
extract.covariogram.parameters(spatial.fit2)

# Kammann and Wand (2003) - circular
spatial.fit3 = fit.spatial.process(data=application.data, method=4, number.of.knots=143,
 	covariogram.model="circular", trend.input="1st")
extract.covariogram.parameters(spatial.fit3)

# Kammann and Wand (2003) - gaussian
spatial.fit10 = fit.spatial.process(data=application.data, method=4, number.of.knots=143,
 	covariogram.model="gaussian", trend.input="1st")
extract.covariogram.parameters(spatial.fit10)

# Kammann and Wand (2003) - spherical
spatial.fit11 = fit.spatial.process(data=application.data, method=4, number.of.knots=143,
 	covariogram.model="circular", trend.input="1st")
extract.covariogram.parameters(spatial.fit11)

# Kammann and Wand (2003) - matern
spatial.fit12 = fit.spatial.process(data=application.data, method=4, number.of.knots=143,
 	covariogram.model="matern", trend.input="1st")
extract.covariogram.parameters(spatial.fit12)



# K-splines - exponential
spatial.fit5 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="exponential", trend.input="1st")
extract.covariogram.parameters(spatial.fit5)

# K-splines - Gaussian
spatial.fit6 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="gaussian", trend.input="1st")
extract.covariogram.parameters(spatial.fit6)

# K-splines - spherical
spatial.fit7 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="spherical", trend.input="1st")
extract.covariogram.parameters(spatial.fit7)

# K-splines - circular
ptm = proc.time()
spatial.fit8 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="circular", trend.input="1st")
proc.time() - ptm

extract.covariogram.parameters(spatial.fit8)
spatial.fit8$df.adjusted
sum(spatial.fit8$df.adjusted)-3


# K-splines - matern (nu=3/2)
spatial.fit9 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="matern", trend.input="1st")
extract.covariogram.parameters(spatial.fit9)


# Thin plate splines
spatial.fit4 = fit.spatial.process(data=application.data, method=6, number.of.knots=143,
 	covariogram.model="thin.plate", trend.input="1st")
extract.covariogram.parameters(spatial.fit4)



# plot covariogram of best model
cov.parms = extract.covariogram.parameters(spatial.fit8)
plot.covariogram(cov.parms,"circular",application.data)



#--- PREDICT SOME FITS ---#
#-------------------------#
#debug(utils:::unpackPkgZip)
#install.packages("SDMTools")

library(SDMTools)
x.pred.vec = seq(100,800,10)
y.pred.vec = seq(0,520,10)
x.pred = expand.grid(x.pred.vec, y.pred.vec)[,1]
y.pred = expand.grid(x.pred.vec, y.pred.vec)[,2]
pred.loc = data.frame(cbind(x.pred,y.pred))

### using K-splines
pred.spline = predict.results(fit=spatial.fit8, data=application.data,
	model="circular", trend.input="1st", locations=pred.loc)
pred.spline.all = cbind(pred.spline)

boot.var.spline = c()
for (i in 1:38){
	if (i!=38) pred.loc2 = pred.loc[ ((i-1)*100+1):(i*100), ]
	if (i==38) pred.loc2 = pred.loc[ ((i-1)*100+1):nrow(pred.loc), ]

	int.result = bootstrap.variance(fit=spatial.fit8, data=application.data, 
		model="circular", method=6, trend.input="1st", locations=pred.loc2,
		boot.samples=500, approximate=TRUE, seeding.bootstrap=98765)
	boot.var.spline = c(boot.var.spline,int.result)
print(i)
}
pred.spline2 = cbind(pred.spline, boot.var.spline)
head(pred.spline2 )

point.in.polygon = pnt.in.poly(pred.loc, parana$borders)
plot.index = matrix(point.in.polygon[,3], nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.mean = matrix(pred.spline.all$pred, nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.mean[plot.index==0] = NA
plot.var = matrix(pred.spline2[,5], nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.var[plot.index==0] = NA

image.plot(x.pred.vec, y.pred.vec, plot.mean, col=topo.colors(200), zlim=c(150,400),
	xlab="Coord X (km)",ylab="Coord Y (km)", main="(b)")
contour(x.pred.vec, y.pred.vec, plot.mean, add = TRUE, drawlabels = TRUE)
polygon(parana$borders,lwd=5)

image.plot(x.pred.vec, y.pred.vec, plot.var, col=topo.colors(200), zlim=c(50, 500), 
	xlab="Coord X (km)",ylab="Coord Y (km)", main="(c)")
contour(x.pred.vec, y.pred.vec, plot.var, add = TRUE, drawlabels = TRUE)
polygon(parana$borders,lwd=5)



# 2 x 2 plot
par(mfrow=c(2,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))

	# plot 1
min(resp);quantile(resp,p=0.25);quantile(resp,p=0.50);quantile(resp,p=0.75);max(resp)
q.index = rep(NA,length(resp))
q.index[162<=resp & resp<234] = 1
q.index[234<=resp & resp<270] = 2
q.index[270<=resp & resp<318] = 3
q.index[318<=resp & resp<414] = 4

plot(x.coord,y.coord, xlim=c(110,800),ylim=c(-80,575), pch=16,col=1,cex=0.01,
xlab="Coord X (km)",ylab="Coord Y (km)", xaxt="n", yaxt="n", main="(a)", cex.lab=1.8, cex.main=1.8)
polygon(parana$borders, lwd=2)
points(x.coord[q.index==1], y.coord[q.index==1], pch=16,col="gray60" ,cex=0.75)
points(x.coord[q.index==2], y.coord[q.index==2], pch=16,col="gray40" ,cex=1.00)
points(x.coord[q.index==3], y.coord[q.index==3], pch=16,col="gray20" ,cex=1.50)
points(x.coord[q.index==4], y.coord[q.index==4], pch=16,col=1 ,cex=2.00)
axis(1,at=seq(200,800,100),c("200","","400","","600","","800"), cex.axis=1.8)
axis(2,at=seq(0,500,100),c("","100","","300","","500"), cex.axis=1.8)
legend(200, 50, c("[162-234[","[234-270[","[270-318[","[318-414["), pch=16, 
   pt.cex=c(0.75,1.00,1.50,2.00), col=c("gray60","gray40","gray20",1),
   ncol=2, bty="n", cex=1.5)

	# plot 2
image.plot(x.pred.vec, y.pred.vec, plot.mean, col=topo.colors(200), xlim=c(110,800),ylim=c(-80,575),
	zlim=c(150,400), xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(b)", cex.lab=1.8, cex.main=1.8)
contour(x.pred.vec, y.pred.vec, plot.mean, add = TRUE, drawlabels = TRUE, labcex=1.0)
polygon(parana$borders,lwd=2)
axis(1,at=seq(200,800,100),c("200","","400","","600","","800"), cex.axis=1.8)
axis(2,at=seq(0,500,100),c("","100","","300","","500"), cex.axis=1.8)

	# plot 3
image.plot(x.pred.vec, y.pred.vec, plot.var, col=topo.colors(200), xlim=c(110,800),ylim=c(-80,575),
	zlim=c(60, 400), xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(c)", cex.lab=1.8, cex.main=1.8)
polygon(parana$borders,lwd=2)
axis(1,at=seq(200,800,100),c("200","","400","","600","","800"), cex.axis=1.8)
axis(2,at=seq(0,500,100),c("","100","","300","","500"), cex.axis=1.8)

	# plot 4
x.dist = seq(0,800,0.10)
cov.exp = 785.6*exp(-x.dist/184.4)
theta = unlist( lapply(x.dist, function(t){min(1,t/330.4)}) )
cov.circ = 706.4*(1 - 2*( theta*sqrt(1-theta^2) +  asin(theta)) / pi)
plot(x.dist, cov.exp, type="l", lwd=4.5, ylim=c(0,1200),
	xlab="Distance (km)",ylab="Covariance",
	xaxt="n", yaxt="n", main="(d)", cex.lab=1.8, cex.main=1.8)
points(0, 785.6 + 385.5, pch=16, cex=1.2)
lines(x.dist, cov.circ, lwd=4.5, col="gray40")
points(0, 706.4 + 415.0, pch=16, cex=1.2, col="gray40")
axis(1,at=seq(100,700,200),c("100","300","500","700"), cex.axis=1.8)
axis(2,at=seq(0,1200,300),seq(0,1200,300), cex.axis=1.8)
legend(80,1200,c("Diggle and Ribeiro Jr (2002):","exponential",
	"Geoadditive model: circular"), lwd=3, 
	col=c(1,NA,"gray40"), bty="n", cex=1.35)


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

knots.vector = c(5,10,15,20,30,40,50,60,80,100,120,140)
results = matrix(nrow=length(knots.vector), ncol=4)

for (i in 1:length(knots.vector)){
fit = fit.spatial.process(data=application.data, method=6,
	number.of.knots = knots.vector[i],
 	covariogram.model="circular", trend.input="1st")
parms = extract.covariogram.parameters(fit)
results[i,1] = parms[1]
results[i,2] = parms[2]
results[i,3] = parms[3]
results[i,4] = parms[4]
}


fit = fit.spatial.process(data=application.data, method=6,
	number.of.knots = NULL,
 	covariogram.model="circular", trend.input="1st")
parms = extract.covariogram.parameters(fit)
results2 = as.matrix(rbind(results, parms))

par(mfrow=c(2,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))
plot( c(knots.vector,143), results2[,1] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",c[s])),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,143), results2[,1], pch=16, cex=1.5)

plot( c(knots.vector,143), results2[,2] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",sigma[epsilon]^2)),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,143), results2[,2], pch=16, cex=1.5)

plot( c(knots.vector,143), results2[,3] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab=expression(paste("parameter  ",tau)),
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,143), results2[,3], pch=16, cex=1.5)

plot( c(knots.vector,143), results2[,4] ,type="l", lwd=2.5,
	xlab="Number of kntos", ylab="AIC value",
	cex.axis=1.8, cex.lab=1.8 )
points(c(knots.vector,143), results2[,4], pch=16, cex=1.5)




















point.in.polygon = pnt.in.poly(pred.loc, parana$borders)
plot.index = matrix(point.in.polygon[,3], nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.var = matrix(pred.spline2[,4], nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.var[plot.index==0] = NA
plot.var2 = matrix(pred.spline2[,5], nrow=length(x.pred.vec), ncol=length(y.pred.vec))
plot.var2[plot.index==0] = NA

par(mfrow=c(1,2), oma=c(0,1,0,1),mar=c(6,6,2,0.5))
image.plot(x.pred.vec, y.pred.vec, plot.var, col=topo.colors(200), xlim=c(110,800),ylim=c(-80,575),
	zlim=c(0, 400), xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(c)", cex.lab=1.8, cex.main=1.8)
polygon(parana$borders,lwd=2)
axis(1,at=seq(200,800,100),c("200","","400","","600","","800"), cex.axis=1.8)
axis(2,at=seq(0,500,100),c("","100","","300","","500"), cex.axis=1.8)

image.plot(x.pred.vec, y.pred.vec, plot.var2, col=topo.colors(200), xlim=c(110,800),ylim=c(-80,575),
	zlim=c(0, 400), xlab="Coord X (km)",ylab="Coord Y (km)", 
	xaxt="n", yaxt="n", main="(c)", cex.lab=1.8, cex.main=1.8)
polygon(parana$borders,lwd=2)
axis(1,at=seq(200,800,100),c("200","","400","","600","","800"), cex.axis=1.8)
axis(2,at=seq(0,500,100),c("","100","","300","","500"), cex.axis=1.8)


