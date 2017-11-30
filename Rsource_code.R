##
## R-functions to run a geostatistical analysis.
##
## Author: Yannick Vandendijck
## Paper: Estimating the spatial covariance structure using 
## the geoadditive model. Environ Ecol Stat (2017) 24:341–361.
## Last update: 20/01/2017
##


read.in.libraries = function(){
	library(geoR)
	library(gstat)
	library(nlme)
	library(fields)
	library(RandomFields)
}


### Mean function for fitting the spatial process
###----------------------------------------------

fit.spatial.process = function(data, method, ini.parms=c(1,1,1), 
				covariogram.model="exponential", number.of.knots=NULL,
				trend.input="cte", nonlinear.trend.input=NULL, cluster.input=NULL,
				tol1=1e-06, tol2=1e-06, max.iter=500, phi.upper.limit=1000){

	# check if input is correct
	if ( !(is.null(nonlinear.trend.input)) & method %in% c(1,2,3) ) print("WARNING: Nonlinear covariates are not used")
	if (method %in% c(1:5)) print("The inputs tol1, tol2, max.iter and phi.upper.limit are not used for this method.")
	if (ini.parms[1]==1 & ini.parms[2]==1 & ini.parms[3]==1) print("The default starting values are used.")
	if ( !(is.null(cluster.input)) & method %in% c(1,2,3) ) print("WARNING: Cluster variable is not used")

	# Fit covariogram with ML
	if (method==1){

		if (ncol(data)<4){
		geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3))
		}
		if (ncol(data)>3){
		geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3), covar.col=c(4:ncol(data)))
		}

		if (trend.input %in% c("cte","1st","2nd")){
		covariogram.results = likfit(geodata.fit,
			ini.cov.pars=c(ini.parms[1],ini.parms[3]), 
			cov.model=covariogram.model, kappa=3/2,
			trend=trend.input,
			nugget=ini.parms[2], 
			fix.nugget=F, nospatial=F, lik.method = "ML")
		}

		if (!(trend.input %in% c("cte","1st","2nd"))){
		covariogram.results = likfit(geodata.fit,
			ini.cov.pars=c(ini.parms[1],ini.parms[3]), 
			cov.model=covariogram.model, kappa=3/2,
			trend=trend.spatial(as.formula(trend.input) , geodata.fit),
			nugget=ini.parms[2], 
			fix.nugget=F, nospatial=F, lik.method = "ML")
		}
	
	print("Estimation method is based on Maximum Likelihood (ML) estimation of the covariogram parameters")

	}

	# Fit covariogram with REML
	if (method==2){

		if (ncol(data)<4){
		geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3))
		}
		if (ncol(data)>3){
		geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3), covar.col=c(4:ncol(data)))
		}

		if (trend.input %in% c("cte","1st","2nd")){
		covariogram.results = likfit(geodata.fit,
			ini.cov.pars=c(ini.parms[1],ini.parms[3]), 
			cov.model=covariogram.model, kappa=3/2,
			trend=trend.input,
			nugget=ini.parms[2], 
			fix.nugget=F, nospatial=F, lik.method = "REML")
		}

		if (!(trend.input %in% c("cte","1st","2nd"))){
		covariogram.results = likfit(geodata.fit,
			ini.cov.pars=c(ini.parms[1],ini.parms[3]), 
			cov.model=covariogram.model, kappa=3/2, 
			trend=trend.spatial(as.formula(trend.input) , geodata.fit),
			nugget=ini.parms[2], 
			fix.nugget=F, nospatial=F, lik.method = "REML")
		}

	print("Estimation method is based on Restricted Maximum Likelihood (REML) estimation of the covariogram parameters")
	
	}

	# Fit covariogram with Weighted Least squares
	if (method==3){
		
		geodata.fit = data
		colnames(geodata.fit)[1] = "Z"

		if (covariogram.model=="exponential") model.input="Exp"
		if (covariogram.model=="gaussian") model.input="Gau"
		if (covariogram.model=="spherical") model.input="Sph"
		if (covariogram.model=="circular") model.input="Cir"
		if (covariogram.model=="matern") model.input="Mat"

		if (trend.input=="cte") formula.input = Z~1
		if (trend.input=="1st"){
			colnames(geodata.fit)[2] = "x"
			colnames(geodata.fit)[3] = "y"
			formula.input = Z~x+y
		}
		if (trend.input=="2nd"){
			colnames(geodata.fit)[2] = "x"
			colnames(geodata.fit)[3] = "y"
			formula.input = Z~x+y+x^2+y^2+x*y
		}
		if (!(trend.input %in% c("cte","1st","2nd"))){
			formula.input = as.formula(paste("Z",trend.input,sep=""))
		}

		coordinates(geodata.fit) = c(2,3)
		vt = variogram(formula.input, data=geodata.fit)
		if (model.input!="Mat"){
		covariogram.results = fit.variogram(vt, vgm(ini.parms[1],model.input,ini.parms[3],ini.parms[2]),
			fit.sills = T, fit.ranges = T, fit.method=2)
		}
		if (model.input=="Mat"){
		covariogram.results = fit.variogram(vt, vgm(ini.parms[1],model.input,ini.parms[3],ini.parms[2],kappa=3/2),
			fit.sills = T, fit.ranges = T, fit.method=2)
		}
		print("Estimation method is based on weighted least squares estimation of the variogram parameters")

	}

	# Fit covariogram parameters with Kammann and Wand (2002) using ML
	if (method==4){

		geodata.fit = data
		if (is.null(number.of.knots)){
			index.knots = 1:nrow(geodata.fit)
			
		}
		if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(geodata.fit[,2:3]))) ){
			set.seed(12345)
			space.filling = cover.design(R=unique(geodata.fit[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
			index.knots = space.filling$best.id
		}
		if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(geodata.fit[,2:3]))) ){
			index.knots = 1:nrow(unique(geodata.fit[,2:3]))
		}

		x.dist = abs(outer(geodata.fit[,2], geodata.fit[,2],"-"))
		y.dist = abs(outer(geodata.fit[,3], geodata.fit[,3],"-"))
		dist.locations = sqrt(x.dist^2 + y.dist^2)
		phi = max(dist.locations)

		x.dist1 = abs(outer(geodata.fit[,2], unique(geodata.fit[,2:3])[index.knots,1],"-"))
		y.dist1 = abs(outer(geodata.fit[,3], unique(geodata.fit[,2:3])[index.knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(geodata.fit[,2:3])[index.knots,1], unique(geodata.fit[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(unique(geodata.fit[,2:3])[index.knots,2], unique(geodata.fit[,2:3])[index.knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)
		
		X = create.X.matrix(geodata.fit, trend.input)
		if (is.null(nonlinear.trend.input)){
			if (is.null(cluster.input)){
				Z = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				covariogram.results = fit.spline.model.ML(geodata.fit, X, Z)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, ncol(Z))
			}
			if (!is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.cluster.matrix(input = cluster.input)
				Z = cbind(Z1,Z2)
				re.block.inds = c(ncol(Z1), ncol(Z2))
				covariogram.results = fit.spline.model.ML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
		}

		if (!(is.null(nonlinear.trend.input))){
			if (is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.matrix.nonlinear(geodata.fit, nonlinear.trend.input)
				Z = cbind(Z1,Z2$matrix)
				re.block.inds = c(ncol(Z1) , Z2$number.of.knots)
				covariogram.results = fit.spline.model.ML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
			if (!is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.matrix.nonlinear(geodata.fit, nonlinear.trend.input)
				Z3 = create.Z.cluster.matrix(input = cluster.input)
				Z = cbind(Z1,Z2$matrix,Z3)
				re.block.inds = c( ncol(Z1) , Z2$number.of.knots, ncol(Z3) )
				covariogram.results = fit.spline.model.ML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
		}

		print("Spline construction is based on Kammann and Wand (2003). Estimated using ML")
	}

	# Fit covariogram parameters with Kammann and Wand (2002) using REML
	if (method==5){

		geodata.fit = data
		if (is.null(number.of.knots)){
			index.knots = 1:nrow(geodata.fit)
			
		}
		if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(geodata.fit[,2:3]))) ){
			set.seed(12345)
			space.filling = cover.design(R=unique(geodata.fit[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
			index.knots = space.filling$best.id
		}
		if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(geodata.fit[,2:3]))) ){
			index.knots = 1:nrow(unique(geodata.fit[,2:3]))
		}

		x.dist = abs(outer(geodata.fit[,2], geodata.fit[,2],"-"))
		y.dist = abs(outer(geodata.fit[,3], geodata.fit[,3],"-"))
		dist.locations = sqrt(x.dist^2 + y.dist^2)
		phi = max(dist.locations)

		x.dist1 = abs(outer(geodata.fit[,2], unique(geodata.fit[,2:3])[index.knots,1],"-"))
		y.dist1 = abs(outer(geodata.fit[,3], unique(geodata.fit[,2:3])[index.knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(geodata.fit[,2:3])[index.knots,1], unique(geodata.fit[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(unique(geodata.fit[,2:3])[index.knots,2], unique(geodata.fit[,2:3])[index.knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)
		
		X = create.X.matrix(geodata.fit, trend.input)
		if (is.null(nonlinear.trend.input)){
			if (is.null(cluster.input)){
				Z = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				covariogram.results = fit.spline.model.REML(geodata.fit, X, Z)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, ncol(Z))
			}
			if (!is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.cluster.matrix(input = cluster.input)
				Z = cbind(Z1,Z2)
				re.block.inds = c(ncol(Z1), ncol(Z2))
				covariogram.results = fit.spline.model.REML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
		}

		if (!(is.null(nonlinear.trend.input))){
			if (is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.matrix.nonlinear(geodata.fit, nonlinear.trend.input)
				Z = cbind(Z1,Z2$matrix)
				re.block.inds = c(ncol(Z1) , Z2$number.of.knots)
				covariogram.results = fit.spline.model.REML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
			if (!is.null(cluster.input)){
				Z1 = create.Z.matrix(geodata.fit, dist.data, dist.knots, covariogram.model, phi)
				Z2 = create.Z.matrix.nonlinear(geodata.fit, nonlinear.trend.input)
				Z3 = create.Z.cluster.matrix(input = cluster.input)
				Z = cbind(Z1,Z2$matrix,Z3)
				re.block.inds = c( ncol(Z1) , Z2$number.of.knots, ncol(Z3) )
				covariogram.results = fit.spline.model.REML.nonlinear.covariates(geodata.fit, X, Z, re.block.inds)
				covariogram.results$phi = phi
				covariogram.results$stop = 1
				covariogram.results$knots = index.knots
				covariogram.results$df.adjusted = NA # degrees.of.freedom(covariogram.results, X, Z, re.block.inds)
			}
		}

		print("Spline construction is based on Kammann and Wand (2003). Estimated using REML")
	}

	# Fit covariogram parameters with our method using ML
	if (method==6){

		geodata.fit = data
		phi.null = ini.parms[3]
		if (is.null(nonlinear.trend.input)){
			covariogram.results = run.doubly.iterative.method.ML(geodata.fit, phi.null, covariogram.model, trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit)
		}
		if (!(is.null(nonlinear.trend.input))){		
			covariogram.results = run.doubly.iterative.method.ML.nonlinear.covariates(geodata.fit, phi.null, covariogram.model, trend.input, nonlinear.trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit)
		}
		if (covariogram.results$stop==0) {print("Model did not converge. You can either increase your maximum number of iterations, decrease the tolerance values or construct another model")}
		else {print("Spline construction is based on our doubly iterative method. Estimated using ML")}
	}

	# Fit covariogram parameters with our method using REML
	if (method==7){

		geodata.fit = data
		phi.null = ini.parms[3]
		if (is.null(nonlinear.trend.input)){
			covariogram.results = run.doubly.iterative.method.REML(geodata.fit, phi.null, covariogram.model, trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit)
		}
		if (!(is.null(nonlinear.trend.input))){	
			covariogram.results = run.doubly.iterative.method.REML.nonlinear.covariates(geodata.fit, phi.null, covariogram.model, trend.input, nonlinear.trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit)
		}
		if (covariogram.results$stop==0) {print("Model did not converge. You can either increase your maximum number of iterations, decrease the tolerance values or construct another model")}
		else {print("Spline construction is based on our doubly iterative method. Estimated using REML")}
	}

return(covariogram.results)

}


### Functions to create important matrices
###---------------------------------------

create.X.matrix <- function(data, input){
	if (input=="cte") X = rep(1,length(data[,1]))
	if (input=="1st") X = cbind( rep(1,length(data[,1])) ,data[,2] ,data[,3] )
	if (input=="2nd") X = cbind( rep(1,length(data[,1])) ,data[,2] ,data[,3], data[,2]^2 ,data[,3]^2 , data[,2]*data[,3])
	if (!(input %in% c("cte","1st","2nd"))){
		t1 = strsplit(input, "[~]")
		t2 = strsplit( t1[[1]][2] , "[+]" )
		X = rep(1,length(data[,1]))
		for (i in 1:length(t2[[1]])){
 			X = cbind( X , eval(parse(text=paste("data$",t2[[1]][i],sep=""))) )
		}
	}
return(X)
}


create.Z.matrix <- function(data, dist.data, dist.knots, model, phi){
	if (model=="exponential"){
		Z = exp(-dist.data/phi)
		K = exp(-dist.knots/phi)
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="gaussian"){
		Z = exp(-(dist.data/phi)^2)
		K = exp(-(dist.knots/phi)^2)
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="spherical"){
		Z = ( 1 - 3/2*(dist.data/phi) + 1/2*(dist.data/phi)^3 ) * (dist.data < phi)
		K = ( 1 - 3/2*(dist.knots/phi) + 1/2*(dist.knots/phi)^3 ) * (dist.knots < phi)
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="circular"){
		theta1v = unlist( lapply(dist.data, function(t){min(1,t/phi)}) )
		theta2v = unlist( lapply(dist.knots, function(t){min(1,t/phi)}) )
		theta1 = matrix(as.vector(theta1v), nrow=nrow(dist.data), ncol=ncol(dist.data))
		theta2 = matrix(as.vector(theta2v), nrow=nrow(dist.knots), ncol=ncol(dist.knots))
		g1 = 2*( theta1*sqrt(1-theta1^2) +  asin(theta1)) / pi
		g2 = 2*( theta2*sqrt(1-theta2^2) +  asin(theta2)) / pi
		Z = 1-g1
		K = 1-g2
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="matern"){
		Z = exp(-dist.data/phi)*(1 + dist.data/phi)
		K = exp(-dist.knots/phi)*(1 + dist.knots/phi)
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="inverse.multiquadratic"){
		Z = 1/sqrt(1 + dist.data^2/phi)
		K = 1/sqrt(1 + dist.knots^2/phi)
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
	if (model=="thin.plate"){
		Z = matrix(0, nrow(dist.data), ncol(dist.data))
		K = matrix(0, nrow(dist.knots), ncol(dist.knots))
		Z[dist.data==0] = (dist.data^2)[dist.data==0]
		Z[dist.data!=0] = (dist.data^2*log(dist.data))[dist.data!=0]
		K[dist.knots==0] = (dist.knots^2)[dist.knots==0]
		K[dist.knots!=0] = (dist.knots^2*log(dist.knots))[dist.knots!=0]
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		Z.tilde = t(solve(sqrt.K,t(Z)))
	}
return(Z.tilde)
}


create.Zpred.matrix <- function(data, pred.locations, model, phi, index.knots){
		x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)

		x.dist2 = abs(outer(pred.locations[,1], unique(data[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(pred.locations[,2], unique(data[,2:3])[index.knots,2],"-"))
		dist.pred.loc = sqrt(x.dist2^2 + y.dist2^2)

		Z.tilde = create.Z.matrix(data=data, dist.data=dist.pred.loc, dist.knots=dist.knots, model=model, phi=phi)
return(Z.tilde)
}


create.Z.matrix.nonlinear <- function(data, nonlinear.trend){
	t1 = strsplit(nonlinear.trend, "[~]")
	t2 = strsplit( t1[[1]][2] , "[+]" )
	Z.tilde = matrix(nrow=nrow(data), ncol=0)
	v = vector()
	for (i in 1:length(t2[[1]])){
		data.points = eval(parse(text=paste("data$",t2[[1]][i],sep="")))
		#knots = as.vector( quantile(unique(data.points) , seq(0,1,length=round(length(data.points)/4)+2) )[-c(1 , round(length(data.points)/4)+2)] )
		knots = as.vector( quantile(unique(data.points) , seq(0,1,length=30) )[-c(1 , 30)] )
		knots.dist = abs(outer(knots,knots,"-"))
		K = knots.dist^3
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))
		points.dist = abs(outer(data.points,knots,"-"))
		Z = points.dist^3
		Z.tilde.append = t(solve(sqrt.K,t(Z)))
		Z.tilde = cbind(Z.tilde,Z.tilde.append)
		v[i] = length(knots) 
	}
results = list(matrix=Z.tilde , number.of.knots = v)
return(results)
}


create.Zpred.matrix.nonlinear <- function(data, prediction.data, nonlinear.trend){
	t1 = strsplit(nonlinear.trend, "[~]")
	t2 = strsplit( t1[[1]][2] , "[+]" )
	Z.tilde = matrix(nrow=nrow(prediction.data), ncol=0)
	v = vector()
	for (i in 1:length(t2[[1]])){
		data.points = eval(parse(text=paste("data$",t2[[1]][i],sep="")))
		#knots = as.vector( quantile(unique(data.points) , seq(0,1,length=round(length(data.points)/4)+2) )[-c(1 , round(length(data.points)/4)+2)] )
		knots = as.vector( quantile(unique(data.points) , seq(0,1,length=30) )[-c(1 , 30)] )
		knots.dist = abs(outer(knots,knots,"-"))
		K = knots.dist^3
		svd.K = svd(K)
		sqrt.K = t(svd.K$v %*% (t(svd.K$u)*sqrt(svd.K$d)))

		pred.data.points = eval(parse(text=paste("prediction.data$",t2[[1]][i],sep="")))
		points.dist = abs(outer(pred.data.points,knots,"-"))
		Z = points.dist^3
		Z.tilde.append = t(solve(sqrt.K,t(Z)))
		Z.tilde = cbind(Z.tilde,Z.tilde.append)
		v[i] = length(knots) 
	}
return(Z.tilde)
}

create.Z.cluster.matrix <- function(input){
	factor.input = as.factor(input)
	cluster.Zmatrix = matrix(0, nrow=length(factor.input), ncol=length(levels(factor.input)))
	for (i in 1:length(factor.input)){
		cluster.Zmatrix[i,] = as.numeric( factor.input[i] == levels(factor.input) )
	}
return(cluster.Zmatrix)
}


create.Zpred.cluster.matrix <- function(input, cluster.locations){
	factor.input = as.factor(input)
	factor.cluster.locations = as.factor(cluster.locations)
	cluster.Zmatrix = matrix(0, nrow=length(cluster.locations), ncol=length(levels(factor.input)))
	for (i in 1:length(cluster.locations)){
		cluster.Zmatrix[i,] = as.numeric( factor.cluster.locations[i] == levels(factor.input) )
	}
return(cluster.Zmatrix)
}


### Functions to fit likelihood models
###-----------------------------------

fit.spline.model.ML <- function(data, fixed.matrix, random.matrix){
	resp = data[,1]
	d = list(z=resp , X=fixed.matrix , Z=random.matrix , g=rep(1,length(resp)))
	d$g = as.factor(d$g)
	fit = lme(z ~ X-1,random=list(g=pdIdent(~-1+Z)),data=d,
		na.action=na.omit,method="ML",
		control=lmeControl(returnObject=TRUE))
return(fit)
}



fit.spline.model.ML.nonlinear.covariates <- function(data, fixed.matrix, random.matrix, block.ind){
	resp = data[,1]
	.GlobalEnv$group =  rep(1,length(resp))
	.GlobalEnv$X = fixed.matrix
	.GlobalEnv$Z = random.matrix
	Z.block <- list()
	v = cumsum(block.ind)
	for (i in 1:length(v)){
		if (i==1){
			Z.block[[i]] = as.formula(paste("~Z[,1:",v[1],"]-1",sep=""))
		}
		if (i==length(v)){
			Z.block[[i]] = as.formula(paste("~Z[,",v[length(v)-1]+1,":",v[length(v)],"]-1",sep=""))
		}
		if (i!=1 & i!=length(v)){
			Z.block[[i]] = as.formula(paste("~Z[,",v[i-1]+1,":",v[i],"]-1",sep=""))
		}
	}

	fitData = groupedData(resp~1|group, data=data.frame(cbind(resp)))
	fit = lme(resp~-1+X , random=pdBlocked(Z.block,pdClass="pdIdent"), method="ML", data=fitData,
		control=lmeControl(returnObject=TRUE))

	remove("group", envir = .GlobalEnv)
	remove("X", envir = .GlobalEnv)
	remove("Z", envir = .GlobalEnv)
return(fit)
}


fit.spline.model.REML <- function(data, fixed.matrix, random.matrix){
	resp = data[,1]
	d = list(z=resp , X=fixed.matrix , Z=random.matrix , g=rep(1,length(resp)))
	d$g = as.factor(d$g)	
	fit = lme(z ~ X-1,random=list(g=pdIdent(~-1+Z)),data=d,
		na.action=na.omit,method="REML",control=lmeControl(returnObject=TRUE))
return(fit)
}


fit.spline.model.REML.nonlinear.covariates <- function(data, fixed.matrix, random.matrix, block.ind){
	resp = data[,1]
	.GlobalEnv$group =  rep(1,length(resp))
	.GlobalEnv$X = fixed.matrix
	.GlobalEnv$Z = random.matrix
	Z.block <- list()
	v = cumsum(block.ind)
	for (i in 1:length(v)){
		if (i==1){
			Z.block[[i]] = as.formula(paste("~Z[,1:",v[1],"]-1",sep=""))
		}
		if (i==length(v)){
			Z.block[[i]] = as.formula(paste("~Z[,",v[length(v)-1]+1,":",v[length(v)],"]-1",sep=""))
		}
		if (i!=1 & i!=length(v)){
			Z.block[[i]] = as.formula(paste("~Z[,",v[i-1]+1,":",v[i],"]-1",sep=""))
		}
	}

	fitData = groupedData(resp~1|group, data=data.frame(cbind(resp)))
	fit = lme(resp~-1+X , random=pdBlocked(Z.block,pdClass="pdIdent"), method="REML", data=fitData,
		control=lmeControl(returnObject=TRUE))
	remove("group", envir = .GlobalEnv)
	remove("X", envir = .GlobalEnv)
	remove("Z", envir = .GlobalEnv)
return(fit)
}


### Functions to fit doubly iterative likelihood models
###----------------------------------------------------

run.doubly.iterative.method.ML <- function(data, phi.null, model, trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit){
	phi.old = phi.null
	stop = 0
	iteration = 0
	phi.opt.vec = vector()
	phi.opt.vec[1] = phi.old
	k = 1
	
	X = create.X.matrix(data, trend.input)

	if (is.null(number.of.knots)){
		index.knots = 1:nrow(data)
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(data[,2:3]))) ){
		set.seed(12345)
		space.filling = cover.design(R=unique(data[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
		index.knots = space.filling$best.id
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(data[,2:3]))) ){
		index.knots = 1:nrow(unique(data[,2:3]))
	}

	x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
	y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
	dist.data = sqrt(x.dist1^2 + y.dist1^2)

	x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
	y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
	dist.knots = sqrt(x.dist2^2 + y.dist2^2)

	while (stop!=1 & iteration < max.iter){
		if (is.null(cluster.input)){
			Z = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			re.block.inds = c(ncol(Z))
			fit = fit.spline.model.ML(data, X, Z)
			maximize.phi = optimize(loglik.function.ML, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit = fit, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
		if (!is.null(cluster.input)){
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.cluster.matrix(input = cluster.input)
			Z = cbind(Z1,Z2)
			re.block.inds = c(ncol(Z1), ncol(Z2))
			fit = fit.spline.model.ML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.ML, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
			
		phi.new = maximize.phi$maximum
		phi.opt.vec[k+1] = phi.new
		if (abs(phi.old - phi.new) < tol2) stop=1
		print(c(iteration, phi.old))
		phi.old = phi.new
		k = k+1
		iteration = iteration + 1
	}
fit$stop.criterion = stop
fit$phi = tail(phi.opt.vec,1)
fit$knots = index.knots
fit$df.adjusted = degrees.of.freedom(fit=fit, fixed.matrix=X, random.matrix=Z, block.inds=re.block.inds) 
return(fit)
}


run.doubly.iterative.method.ML.nonlinear.covariates <- function(data, phi.null, model, trend.input, nonlinear.trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit){
	phi.old = phi.null
	stop = 0
	iteration = 0
	phi.opt.vec = vector()
	phi.opt.vec[1] = phi.old
	k = 1

	X = create.X.matrix(data, trend.input)

	if (is.null(number.of.knots)){
		index.knots = 1:nrow(data)
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(data[,2:3]))) ){
		set.seed(12345)
		space.filling = cover.design(R=unique(data[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
		index.knots = space.filling$best.id
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(data[,2:3]))) ){
		index.knots = 1:nrow(unique(data[,2:3]))
	}
	
	x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
	y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
	dist.data = sqrt(x.dist1^2 + y.dist1^2)

	x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
	y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
	dist.knots = sqrt(x.dist2^2 + y.dist2^2)
	
	while (stop!=1 & iteration < max.iter){
		if (is.null(cluster.input)){	
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.matrix.nonlinear(data, nonlinear.trend.input)
			Z = cbind(Z1,Z2$matrix)
			re.block.inds <- c(ncol(Z1) , Z2$number.of.knots)
			fit = fit.spline.model.ML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.ML.nonlinear.covariates, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, nonlinear.trend=nonlinear.trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
		if (!is.null(cluster.input)){
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.matrix.nonlinear(data, nonlinear.trend.input)
			Z3 = create.Z.cluster.matrix(input = cluster.input)
			Z = cbind(Z1, Z2$matrix, Z3)
			re.block.inds <- c(ncol(Z1), Z2$number.of.knots, ncol(Z3))
			fit = fit.spline.model.ML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.ML.nonlinear.covariates, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, nonlinear.trend=nonlinear.trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}		
		
		phi.new = maximize.phi$maximum
		phi.opt.vec[k+1] = phi.new
		if (abs(phi.old - phi.new) < tol2) stop=1
		print(c(iteration, phi.old))
		phi.old = phi.new
		k = k+1
		iteration = iteration + 1
	}
fit$stop.criterion = stop
fit$phi = tail(phi.opt.vec,1)
fit$knots = index.knots
fit$df.adjusted = degrees.of.freedom(fit=fit, fixed.matrix=X, random.matrix=Z, block.inds=re.block.inds) 
return(fit)
}


run.doubly.iterative.method.REML <- function(data, phi.null, model, trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit){
	phi.old = phi.null
	stop = 0
	iteration = 0
	phi.opt.vec = vector()
	phi.opt.vec[1] = phi.old
	k = 1

	X = create.X.matrix(data, trend.input)

	if (is.null(number.of.knots)){
		index.knots = 1:nrow(data)
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(data[,2:3]))) ){
		set.seed(12345)
		space.filling = cover.design(R=unique(data[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
		index.knots = space.filling$best.id
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(data[,2:3]))) ){
		index.knots = 1:nrow(unique(data[,2:3]))
	}
	
	x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
	y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
	dist.data = sqrt(x.dist1^2 + y.dist1^2)

	x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
	y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
	dist.knots = sqrt(x.dist2^2 + y.dist2^2)

	while (stop!=1 & iteration < max.iter){
		if (is.null(cluster.input)){
			Z = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			re.block.inds = c(ncol(Z))
			fit = fit.spline.model.REML(data, X, Z)
			maximize.phi = optimize(loglik.function.REML, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit = fit, trend.input=trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
		if (!is.null(cluster.input)){
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.cluster.matrix(input = cluster.input)
			Z = cbind(Z1,Z2)
			re.block.inds = c(ncol(Z1), ncol(Z2))
			fit = fit.spline.model.REML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.REML, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, trend.input=trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
			
		phi.new = maximize.phi$maximum
		phi.opt.vec[k+1] = phi.new
		if (abs(phi.old - phi.new) < tol2) stop=1
		print(c(iteration, phi.old))
		phi.old = phi.new
		k = k+1
		iteration = iteration + 1
	}
fit$stop.criterion = stop
fit$phi = tail(phi.opt.vec,1)
fit$knots = index.knots
fit$df.adjusted = NA # degrees.of.freedom(fit=fit, fixed.matrix=X, random.matrix=Z, block.inds=re.block.inds)
return(fit)
}


run.doubly.iterative.method.REML.nonlinear.covariates <- function(data, phi.null, model, trend.input, nonlinear.trend.input, cluster.input, number.of.knots, tol1, tol2, max.iter, phi.upper.limit){
	phi.old = phi.null
	stop = 0
	iteration = 0
	phi.opt.vec = vector()
	phi.opt.vec[1] = phi.old
	k = 1

	X = create.X.matrix(data, trend.input)

	if (is.null(number.of.knots)){
		index.knots = 1:nrow(data)
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots < nrow(unique(data[,2:3]))) ){
		set.seed(12345)
		space.filling = cover.design(R=unique(data[,2:3]), nd=number.of.knots, nn=TRUE, nruns=1)
		index.knots = space.filling$best.id
	}
	if ( (!is.null(number.of.knots)) & (number.of.knots == nrow(unique(data[,2:3]))) ){
		index.knots = 1:nrow(unique(data[,2:3]))
	}
	
	x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
	y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
	dist.data = sqrt(x.dist1^2 + y.dist1^2)

	x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
	y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
	dist.knots = sqrt(x.dist2^2 + y.dist2^2)
	
	while (stop!=1 & iteration < max.iter){
		if (is.null(cluster.input)){	
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.matrix.nonlinear(data, nonlinear.trend.input)
			Z = cbind(Z1,Z2$matrix)
			re.block.inds <- c(ncol(Z1) , Z2$number.of.knots)
			fit = fit.spline.model.REML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.REML.nonlinear.covariates, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, trend.input=trend.input, nonlinear.trend=nonlinear.trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}
		if (!is.null(cluster.input)){
			Z1 = create.Z.matrix(data, dist.data, dist.knots, model, phi.old)
			Z2 = create.Z.matrix.nonlinear(data, nonlinear.trend.input)
			Z3 = create.Z.cluster.matrix(input = cluster.input)
			Z = cbind(Z1, Z2$matrix, Z3)
			re.block.inds <- c(ncol(Z1), Z2$number.of.knots, ncol(Z3))
			fit = fit.spline.model.REML.nonlinear.covariates(data, X, Z, re.block.inds)
			maximize.phi = optimize(loglik.function.REML.nonlinear.covariates, c(0.001,phi.upper.limit), data=data, dist.data=dist.data, dist.knots=dist.knots, model=model, fit=fit, trend.input=trend.input, nonlinear.trend=nonlinear.trend.input, cluster.input=cluster.input, tol=tol1, maximum=TRUE)
		}		
		
		phi.new = maximize.phi$maximum
		phi.opt.vec[k+1] = phi.new
		if (abs(phi.old - phi.new) < tol2) stop=1
		print(c(iteration, phi.old))
		phi.old = phi.new
		k = k+1
		iteration = iteration + 1
	}
fit$stop.criterion = stop
fit$phi = tail(phi.opt.vec,1)
fit$knots = index.knots
fit$df.adjusted = NA # degrees.of.freedom(fit=fit, fixed.matrix=X, random.matrix=Z, block.inds=re.block.inds)
return(fit)
}


loglik.function.ML <- function(phi.opt, data, dist.data, dist.knots, model, fit, cluster.input){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))
	M1 = sigmasq.epsilon*diag(nrow(data))

	if (is.null(cluster.input)){
		Z.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		M2 = sigmasq.u*diag(ncol(Z.opt))
	}
	if (!is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.cluster.matrix(input = cluster.input)
		Z.opt = cbind(Z1.opt, Z2.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt))) )
	}
		
	Sigma.matrix = Z.opt%*%M2%*%t(Z.opt) + M1
	ll1 = -1/2*determinant(Sigma.matrix, logarithm=TRUE)$modulus[1]
	ll2 = -nrow(data)/2 * log(2*pi)
	ll3 = -1/2*(fit$residuals[,1]%*%solve(Sigma.matrix,fit$residuals[,1]))
return(ll1 + ll2 + ll3)
}


loglik.function.ML.nonlinear.covariates <- function(phi.opt, data, dist.data, dist.knots, model, fit, nonlinear.trend, cluster.input){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))
	M1 = sigmasq.epsilon*diag(nrow(data))

	if (is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.matrix.nonlinear(data, nonlinear.trend)$matrix
		Z.opt = cbind(Z1.opt,Z2.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt))) )
	}
	if (!is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.matrix.nonlinear(data, nonlinear.trend)$matrix
		Z3.opt = create.Z.cluster.matrix(input = cluster.input)
		Z.opt = cbind(Z1.opt, Z2.opt, Z3.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt),ncol(Z3.opt))) )
	}	

	Sigma.matrix = Z.opt%*%M2%*%t(Z.opt) + M1
	ll1 = -1/2*determinant(Sigma.matrix, logarithm=TRUE)$modulus[1]
	ll2 = -nrow(data)/2 * log(2*pi)
	ll3 = -1/2*(fit$residuals[,1]%*%solve(Sigma.matrix,fit$residuals[,1]))
return(ll1 + ll2 + ll3)
}


loglik.function.REML <- function(phi.opt, data, dist.data, dist.knots, model, fit, trend.input, cluster.input){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))
	M1 = sigmasq.epsilon*diag(nrow(data))
	X.opt = create.X.matrix(data, trend.input)
	
	if (is.null(cluster.input)){
		Z.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		M2 = sigmasq.u*diag(ncol(Z.opt))
	}
	if (!is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.cluster.matrix(input = cluster.input)
		Z.opt = cbind(Z1.opt, Z2.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt))) )
	}
		
	Sigma.matrix = Z.opt%*%M2%*%t(Z.opt) + M1
	ll1 = -1/2*determinant(Sigma.matrix, logarithm=TRUE)$modulus[1]
	ll2 = -1/2*determinant( (t(X.opt)%*%solve(Sigma.matrix,X.opt)) , logarithm=TRUE)$modulus[1]
	if (class(X.opt)=="numeric") ll3 = -(nrow(data) - 1)/2 * log(2*pi)
	if (class(X.opt)!="numeric") ll3 = -(nrow(data) - ncol(X.opt))/2 * log(2*pi)
	ll4 = -1/2*(fit$residuals[,1]%*%solve(Sigma.matrix,fit$residuals[,1]))
return(ll1 + ll2 + ll3 + ll4)
}


loglik.function.REML.nonlinear.covariates <- function(phi.opt, data, dist.data, dist.knots, model, fit, trend.input, nonlinear.trend, cluster.input){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))
	M1 = sigmasq.epsilon*diag(nrow(data))
	X.opt = create.X.matrix(data, trend.input)

	if (is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.matrix.nonlinear(data, nonlinear.trend)$matrix
		Z.opt = cbind(Z1.opt,Z2.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt))) )
	}
	if (!is.null(cluster.input)){
		Z1.opt = create.Z.matrix(data, dist.data, dist.knots, model, phi.opt)
		Z2.opt = create.Z.matrix.nonlinear(data, nonlinear.trend)$matrix
		Z3.opt = create.Z.cluster.matrix(input = cluster.input)
		Z.opt = cbind(Z1.opt, Z2.opt, Z3.opt)
		M2 = diag( rep(sigmasq.u, c(ncol(Z1.opt),ncol(Z2.opt),ncol(Z3.opt))) )
	}	

	Sigma.matrix = Z.opt%*%M2%*%t(Z.opt) + M1
	ll1 = -1/2*determinant(Sigma.matrix, logarithm=TRUE)$modulus[1]
	ll2 = -1/2*determinant( (t(X.opt)%*%solve(Sigma.matrix,X.opt)) , logarithm=TRUE)$modulus[1]
	ll3 = -( nrow(data) - ncol(X.opt))/2 * log(2*pi)
	ll4 = -1/2*(fit$residuals[,1]%*%solve(Sigma.matrix,fit$residuals[,1]))
return(ll1 + ll2 + ll3 + ll4)
}


### Functions to extract covarigram parameters and plot it
###-------------------------------------------------------

extract.covariogram.parameters <- function(input.results){
covariogram.parms = vector()
	if (class(input.results)[1]=="likGRF"){
		covariogram.parms[1] = input.results$sigmasq
		covariogram.parms[2] = input.results$tausq
		covariogram.parms[3] = input.results$phi
		covariogram.parms[4] = AIC(input.results)
	}

	if (class(input.results)[1]=="lme"){
		covariogram.parms[1] = input.results$sigma^2*exp(2*unlist(input.results$modelStruct))[1]
		covariogram.parms[2] = input.results$sigma^2
		covariogram.parms[3] = input.results$phi
		covariogram.parms[4] = AIC(input.results)+2
		v = input.results$sigma^2*exp(2*unlist(input.results$modelStruct))[-1]
		if (length(v)!=0){
			covariogram.parms[5:(4+length(v))] = v
		}
	}

	if (class(input.results)[1]=="variogramModel"){
		covariogram.parms[1] = input.results[2,2]
		covariogram.parms[2] = input.results[1,2]
		covariogram.parms[3] = input.results[2,3]
		covariogram.parms[4] = NA
	}

	print(paste("Estimated sill is", covariogram.parms[1]))		
	print(paste("Estimated error term is", covariogram.parms[2]))
	print(paste("Estimated phi parameter is", covariogram.parms[3]))
	print(paste("The AIC value is", covariogram.parms[4]))
	if (!(is.na(covariogram.parms[5]))){
	print(paste("Variance parameters of non-linear covariates", covariogram.parms[-c(1:4)])) }

return(covariogram.parms)
}


plot.covariogram <- function(parameters, model, data){
	x.dist = abs(outer(data[,2],data[,2],"-"))
	y.dist = abs(outer(data[,3],data[,3],"-"))
	max.dist = max(sqrt(x.dist^2 + y.dist^2))
	x = seq(0,max.dist,max.dist/1000)
	if (model=="exponential"){
		y = parameters[1]*exp( - x / parameters[3])
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Exponential covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
	if (model=="gaussian"){
		y = parameters[1]*exp( - (x / parameters[3])^2 )
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Gaussian covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
	if (model=="spherical"){
		y = parameters[1]*(1 - 3/2*(x/parameters[3]) + 1/2*(x/parameters[3])^3 ) * (x < parameters[3])
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Spherical covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
	if (model=="circular"){
		theta1 = unlist( lapply(x, function(t){min(1,t/parameters[3])}) )
		g1 = 2*( theta1*sqrt(1-theta1^2) +  asin(theta1)) / pi
		y = parameters[1]*(1-g1)
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Circular covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
	if (model=="matern"){
		y = parameters[1]*exp( - x / parameters[3]) * (1 + x / parameters[3])
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Matern (nu=3/2) covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
	if (model=="inverse.multiquadratic"){
		y = parameters[1] * (1/sqrt(1 + x^2/parameters[3]))
		y[0] = parameters[1] + parameters[2]
		plot(x,y,type="l", ylim=c(0,(parameters[1] + parameters[2])*1.05),lwd=3,
		ylab="Covariogram", xlab="Distance", main="Inverse multiquadratic covariogram function")
		points(0,parameters[1] + parameters[2],pch=16)
	}
}


### Functions to predict results
###-----------------------------

predict.results <- function(fit, data, model, trend.input, nonlinear.trend.input=NULL, locations, cov.locations=NULL, cluster.input=NULL, cluster.locations=NULL){
	parameters = extract.covariogram.parameters(fit)

	if (class(fit)[1]!="lme"){
		if (ncol(data)<4){
			geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3))
		}
		if (ncol(data)>3){
			geodata.fit = as.geodata(data, data.col=1, coords.col=c(2,3), covar.col=c(4,ncol(data)))
		}

		if (is.null(cov.locations)){
			prediction.data = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data = cbind(locations,cov.locations)
		}

		if (trend.input %in% c("cte","1st","2nd")){
		kriging = krige.conv(geodata.fit, loc=locations,
			krige=krige.control(type.krige="ok", cov.pars=
			c(parameters[1],parameters[3]),cov.model=model, lambda=1, kappa=3/2,
			nugget=parameters[2], trend.d=trend.input, trend.l=trend.input, 
			micro.scale=0), output=output.control(signal=T))
		}

		if (!(trend.input %in% c("cte","1st","2nd"))){
		t1 = strsplit(trend.input,"[~]")
		t2 = strsplit(t1[[1]][2],"[+]")
		v1 = vector()
		v2 = vector()
		for ( i in 1:length(t2[[1]]) ){
			v1[i] = paste("data",t2[[1]][i],sep="$")
			v2[i] = paste("prediction.data",t2[[1]][i],sep="$")
		}

		kriging = krige.conv(geodata.fit, loc=locations,
			krige=krige.control(type.krige="ok", cov.pars=
			c(parameters[1],parameters[3]),cov.model=model, lambda=1, kappa=3/2,
			nugget=parameters[2],
			trend.d=as.formula(paste("~",paste(v1,collapse="+"),sep="")),
			trend.l=as.formula(paste("~",paste(v2,collapse="+"),sep="")),
			micro.scale=0), output=output.control(signal=T))
		}
	pred.matrix = data.frame(cbind(locations , kriging$predict , kriging$krige.var))
	names(pred.matrix) = c(names(locations),"pred","var")
	}

	if (class(fit)[1]=="lme"){
		beta.hat = fit$coef$fixed
		prediction.data1 = rep(0,nrow(locations))

		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations, cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		X.pred = create.X.matrix(prediction.data, trend.input)

		if (is.null(cluster.locations)){
			if (is.null(nonlinear.trend.input)){
				u.hat = unlist(fit$coef$rand)
				Z.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				pred.spline = as.matrix(X.pred)%*%as.matrix(beta.hat,length(beta.hat),1) + Z.pred%*%u.hat
				var.splines = calculate.variance.spline(fit, data, model, trend.input, nonlinear.trend.input, locations, cov.locations, parameters[3], fit$knots, cluster.input, cluster.locations)
			}
			if (!(is.null(nonlinear.trend.input))){
				u.hat = as.vector(unlist(fit$coef$rand))
				Z1 = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2 = create.Zpred.matrix.nonlinear(data, cov.locations, nonlinear.trend.input)
				Z.pred = cbind(Z1,Z2)
				pred.spline = as.matrix(X.pred)%*%as.matrix(beta.hat,length(beta.hat),1) + Z.pred%*%u.hat
				var.splines = calculate.variance.spline(fit, data, model, trend.input, nonlinear.trend.input, locations, cov.locations, parameters[3], fit$knots, cluster.input, cluster.locations)	
			}
		}

		if (!is.null(cluster.locations)){
			if (is.null(nonlinear.trend.input)){
				u.hat = unlist(fit$coef$rand)
				Z1 = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2 = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1,Z2)
				pred.spline = as.matrix(X.pred)%*%as.matrix(beta.hat,length(beta.hat),1) + Z.pred%*%u.hat
				var.splines = calculate.variance.spline(fit, data, model, trend.input, nonlinear.trend.input, locations, cov.locations, parameters[3], fit$knots, cluster.input, cluster.locations)
			}
			if (!(is.null(nonlinear.trend.input))){
				u.hat = unlist(fit$coef$rand)
				Z1 = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2 = create.Zpred.matrix.nonlinear(data, cov.locations, nonlinear.trend.input)
				Z3 = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1,Z2,Z3)
				pred.spline = as.matrix(X.pred)%*%as.matrix(beta.hat,length(beta.hat),1) + Z.pred%*%u.hat
				var.splines = calculate.variance.spline(fit, data, model, trend.input, nonlinear.trend.input, locations, cov.locations, parameters[3], fit$knots, cluster.input, cluster.locations)
			}
		}

	pred.matrix = data.frame(cbind(locations , pred.spline , var.splines ))
	names(pred.matrix) = c(names(locations),"pred","var")
	}

return(pred.matrix)
}



calculate.variance.spline <- function(fit, data, model, trend.input, nonlinear.trend.input, locations, cov.locations, optimal.phi, index.knots, cluster.input, cluster.locations){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))

	if (is.null(nonlinear.trend.input)){
		X.fit = create.X.matrix(data, trend.input)

		x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
		y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)
	
		Sigma.inv = solve( sigmasq.epsilon*diag(nrow(data)) )
		
		if (is.null(cluster.input)){
			Z.fit = create.Z.matrix(data, dist.data, dist.knots, model, optimal.phi)
			C.fit = cbind(X.fit, Z.fit)
			B = matrix(0, ncol(Z.fit)+length(fit$coef$fixed) , ncol(Z.fit)+length(fit$coef$fixed) )
			i1 = length(fit$coef$fixed) + 1; i2 = ncol(Z.fit) + length(fit$coef$fixed)
			B[i1:i2,i1:i2] = solve( sigmasq.u*diag(ncol(Z.fit)) ) 
		}

		if (!(is.null(cluster.input))){
			Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, optimal.phi)
			Z2.fit = create.Z.cluster.matrix(cluster.input)
			C.fit = cbind(X.fit, Z1.fit, Z2.fit)
			B = matrix(0, ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed) , ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed) )
			i1 = length(fit$coef$fixed)+1; i2 = ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed)
			B[i1:i2,i1:i2] = solve( diag(rep(sigmasq.u,c(ncol(Z1.fit),ncol(Z2.fit)))) ) 
		}
		
		M1 = solve( t(C.fit)%*%Sigma.inv%*%C.fit + B )
		covariance.matrix = M1 					###%*% t(C.fit)%*%Sigma.inv%*%C.fit %*% M1

		prediction.data1 = rep(0,nrow(locations))
		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations,cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		X.pred = create.X.matrix(prediction.data, trend.input)

		if (is.null(cluster.locations)){
			Z.pred = create.Zpred.matrix(data, locations, model, optimal.phi, fit$knots)
		}
		if (!(is.null(cluster.locations))){
			Z1 = create.Zpred.matrix(data, locations, model, optimal.phi, index.knots)
			Z2 = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
			Z.pred = cbind(Z1,Z2)
		}

		C.pred = cbind(X.pred, Z.pred)
		result= diag(C.pred %*% covariance.matrix %*% t(C.pred))

	}

	
	if (!(is.null(nonlinear.trend.input))){
		X.fit = create.X.matrix(data, trend.input)
		

		x.dist1 = abs(outer(data[,2], unique(data[,2:3])[index.knots,1],"-"))
		y.dist1 = abs(outer(data[,3], unique(data[,2:3])[index.knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(data[,2:3])[index.knots,1], unique(data[,2:3])[index.knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[index.knots,2], unique(data[,2:3])[index.knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)

		Sigma.inv = solve( sigmasq.epsilon*diag(nrow(data)) )
		
		if (is.null(cluster.input)){
			Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, optimal.phi)
			Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
			C.fit = cbind(X.fit, Z1.fit, Z2.fit)
			B = matrix(0, ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed) , ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed) )
			i1 = length(fit$coef$fixed)+1; i2 = ncol(Z1.fit)+ncol(Z2.fit)+length(fit$coef$fixed)
			B[i1:i2,i1:i2] = solve( diag(rep(sigmasq.u,c(ncol(Z1.fit),ncol(Z2.fit)))) ) 
		}

		if (!(is.null(cluster.input))){
			Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, optimal.phi)
			Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
			Z3.fit = create.Z.cluster.matrix(cluster.input)
			C.fit = cbind(X.fit, Z1.fit, Z2.fit, Z3.fit)
			B = matrix(0, ncol(Z1.fit)+ncol(Z2.fit)+ncol(Z3.fit)+length(fit$coef$fixed) , ncol(Z1.fit)+ncol(Z2.fit)+ncol(Z3.fit)+length(fit$coef$fixed) )
			i1 = length(fit$coef$fixed)+1; i2 = ncol(Z1.fit)+ncol(Z2.fit)+ncol(Z3.fit)+length(fit$coef$fixed)
			B[i1:i2,i1:i2] = solve( diag(rep(sigmasq.u,c(ncol(Z1.fit),ncol(Z2.fit),ncol(Z3.fit)))) ) 
		}

		M1 = solve( t(C.fit)%*%Sigma.inv%*%C.fit + B )
		covariance.matrix = M1 					###%*% t(C.fit)%*%Sigma.inv%*%C.fit %*% M1

		prediction.data1 = rep(0,nrow(locations))
		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations,cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		X.pred = create.X.matrix(prediction.data, trend.input)
		
		if (is.null(cluster.locations)){
			Z1 = create.Zpred.matrix(data, locations, model, optimal.phi, fit$knots)
			Z2 = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
			Z.pred = cbind(Z1,Z2)
		}
		if (!(is.null(cluster.locations))){
			Z1 = create.Zpred.matrix(data, locations, model, optimal.phi, index.knots)
			Z2 = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
			Z3 = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
			Z.pred = cbind(Z1,Z2,Z3)
		}

		C.pred = cbind(X.pred, Z.pred)
		result= diag(C.pred %*% covariance.matrix %*% t(C.pred))
	}
return(result)
}


### Functions to perform bootstrap calculations of variance
###--------------------------------------------------------

bootstrap.variance <- function(fit, data, model, method, trend.input, nonlinear.trend.input=NULL, locations, cov.locations=NULL, cluster.input=NULL, cluster.locations=NULL, boot.samples=100, approximate=TRUE ,seeding.bootstrap=12345){
	if (!(method %in% c(1,2,6,7))) print("Bootstrap variance function only implemented for methods 1, 2, 6 and 7.")

	else{
		parameters = extract.covariogram.parameters(fit)
		xboot = c(data[,2], locations[,1])
		yboot = c(data[,3], locations[,2])

		prediction.data1 = rep(0,nrow(locations))
		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations,cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		names(prediction.data) = names(data)
		all.data = rbind(data, prediction.data)
		Xmatrix = create.X.matrix(all.data, trend.input)
		result.matrix = matrix(nrow=nrow(locations), ncol=boot.samples)
	
	if (class(fit)[1]!="lme"){
		set.seed(seeding.bootstrap)
		b = 1
		k = 1
		while(b < (boot.samples+1) & k<(2*boot.samples) ){
			fixed.par = fit$beta
			if (length(fixed.par)==1) mu = Xmatrix * fixed.par
			if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par

			Sboot = simulate.Gaussian.Field(model, parameters, xboot, yboot)
			rand.numb.boot = rnorm(length(xboot), mean=0, sd=sqrt(parameters[2]))
			boot.resp = Sboot + mu + rand.numb.boot

			boot.sample = data.frame(cbind(boot.resp, all.data[,-1]))
			boot.data = boot.sample[1:nrow(data),]
			true.value = (Sboot + mu)[(nrow(data)+1):nrow(all.data)]

			if (approximate!="TRUE"){
				fit.boot = try( fit.spatial.process(data=boot.data, method=method, covariogram.model=model,					
					ini.parms=round(parameters[1:3],3), trend.input=trend.input, 
					phi.upper.limit=parameters[3]*10) )
			} 
			if (approximate=="TRUE"){
				fit.boot = fit
			}

			if (class(fit.boot)!="try-error"){
				pred.boot = predict.results(fit.boot, data=boot.data, model=model,
					trend.input=trend.input, nonlinear.trend.input=nonlinear.trend.input,
					locations=locations, cov.locations=cov.locations)

				result.matrix[,b] = (pred.boot[,3] - true.value)^2
				print(paste("Bootstrap run ",b,sep=""))
				b=b+1
			}
			k=k+1		
		}
	} # end of bootstrap for kriging methods

	if (class(fit)[1]=="lme"){
		fixed.par = fit$coefficients$fixed
		X.fit = create.X.matrix(data, trend.input)
		X.pred = create.X.matrix(prediction.data, trend.input)

		x.dist1 = abs(outer(data[,2], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist1 = abs(outer(data[,3], unique(data[,2:3])[fit$knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(data[,2:3])[fit$knots,1], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[fit$knots,2], unique(data[,2:3])[fit$knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)

		random.par = as.vector(unlist(fit$coefficients$random))
		if (is.null(nonlinear.trend.input)){
			if (is.null(cluster.input)){
				Z.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				block.inds = NULL				
				number.of.knots = ncol(Z.fit)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1.pred, Z2.pred)
				parm.cluster = parameters[5]
			}
		}

		if (!(is.null(nonlinear.trend.input))){
			if (is.null(cluster.input)){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				Z.pred = cbind(Z1.pred, Z2.pred)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z3.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit, Z3.fit)
				block.inds = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))
				number.of.knots = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				Z3.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1.pred, Z2.pred, Z3.pred)
				parm.cluster = parameters[5]
			}
		}

		C.pred = cbind(X.pred,Z.pred)
		C.fit = cbind(X.fit, Z.fit)
		Z.all = rbind(Z.fit,Z.pred)

		set.seed(seeding.bootstrap)
		b = 1
		k = 1
		while( b<(boot.samples+1) & k<(2*boot.samples) ){

			if (is.null(nonlinear.trend.input)){
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par
				}
				if (!(is.null(cluster.input))){
					Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
					random.par[(ncol(Z1.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z2.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}						
			}
			if (!(is.null(nonlinear.trend.input))){
				Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
				if (!(is.null(cluster.input))){
					random.par[(ncol(Z1.fit)+ncol(Z2.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z3.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
			}

			Sboot = simulate.Gaussian.Field(model, parameters, xboot, yboot)
			rand.numb.boot = rnorm(length(xboot), mean=0, sd=sqrt(parameters[2]))
			boot.resp = Sboot + mu + rand.numb.boot

			boot.sample = data.frame(cbind(boot.resp, all.data[,-1]))
			boot.data = boot.sample[1:nrow(data),]
			true.value = (Sboot + mu)[(nrow(data)+1):nrow(all.data)]

			if (approximate != "TRUE"){
				fit.boot = try( fit.bootstrap.lme.model(method, boot.data, X.fit, Z.fit, block.inds) )
			} 
			if (approximate == "TRUE"){
				fit.boot = try( bootstrap.lme.approx.model(fit, boot.data, X.fit, Z.fit, number.of.knots) )
			}

			if (class(fit.boot)!="try-error"){
				par.boot = c(fit.boot$coef$fixed,unlist(fit.boot$coef$random))
				result.matrix[,b] = (C.pred%*%par.boot - true.value)^2
				print(paste("Bootstrap run ",b,sep=""))
				b=b+1
			}
			k=k+1		
		}
	} # end of bootstrap for spline method
	
output = apply(result.matrix,1,mean)
return = output
return(return)
	} # end of else statement
}


fit.bootstrap.lme.model <- function(method, data, fixed.matrix, random.matrix, block.ind){
	if (is.null(block.ind)){
		if(method==6){
			fit.boot.model = fit.spline.model.ML(data, fixed.matrix, random.matrix)
		}
		if(method==7){
			fit.boot.model = fit.spline.model.REML(data, fixed.matrix, random.matrix)
		}	
	}

	if (!(is.null(block.ind))){
		if(method==6){
			fit.boot.model = fit.spline.model.ML.nonlinear.covariates(data, fixed.matrix, random.matrix, block.ind)
		}
		if(method==7){
			fit.boot.model = fit.spline.model.REML.nonlinear.covariates(data, fixed.matrix, random.matrix, block.ind)
		}	
	}
	return(fit.boot.model)
}


bootstrap.lme.approx.model <- function(fit, data, fixed.matrix, random.matrix, number.of.knots){
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))

	R = sigmasq.epsilon*diag(nrow(data))
	G = diag( rep(sigmasq.u,number.of.knots))
	V = random.matrix %*% G %*% t(random.matrix) + R
		
	opt.beta = solve( t(fixed.matrix)%*%solve(V,fixed.matrix) , t(fixed.matrix) ) %*% solve(V, data[,1])
	opt.u = G%*%t(random.matrix)%*%solve(V,(data[,1] - fixed.matrix%*%opt.beta))

	fit$coef$fixed = opt.beta
	fit$coef$random = opt.u
return(fit)
}


simulate.Gaussian.Field <- function(model.input, parameters, xboot, yboot){

	if (model.input=="exponential"){
		model = RMexp(var=parameters[1], scale=parameters[3]) + RMnugget(var=0) + RMtrend(mean=0)
		S.vec = RFsimulate(model, x=xboot, y=yboot, grid=FALSE)
		Sboot = S.vec@data$variable1
	}

	if (model.input=="matern"){
		model = RMwhittle(nu=1.5, var=parameters[1], scale=parameters[3]) + RMnugget(var=0) + RMtrend(mean=0)
		S.vec = RFsimulate(model, x=xboot, y=yboot, grid=FALSE)
		Sboot = S.vec@data$variable1
	}

	if (model.input=="gaussian"){
		model = RMgauss(var=parameters[1], scale=parameters[3]) + RMnugget(var=0) + RMtrend(mean=0)
		S.vec = RFsimulate(model, x=xboot, y=yboot, grid=FALSE)
		Sboot = S.vec@data$variable1
	}

	if (model.input=="circular"){
		model = RMcircular(var=parameters[1], scale=parameters[3]) + RMnugget(var=0) + RMtrend(mean=0)
		S.vec = RFsimulate(model, x=xboot, y=yboot, grid=FALSE)
		Sboot = S.vec@data$variable1
	}

	if (model.input=="spherical"){
		model = RMspheric(var=parameters[1], scale=parameters[3]) + RMnugget(var=0) + RMtrend(mean=0)
		S.vec = RFsimulate(model, x=xboot, y=yboot, grid=FALSE)
		Sboot = S.vec@data$variable1
	}
return(Sboot)
}




### Function to obtain the spatial effect
###--------------------------------------

spatial.effect <- function(fit, model, data, locations){
	Zspatial = create.Zpred.matrix(data, locations, model, fit$phi, fit$knots)
	random.par = as.vector(unlist(fit$coefficients$random))
	results = Zspatial %*% random.par[1:length(fit$knots)]
return(results)
}




### Function to obtain the degrees of freedom of the fit
###-----------------------------------------------------

degrees.of.freedom <- function(fit, fixed.matrix, random.matrix, block.inds){
	C = cbind(fixed.matrix, random.matrix)
	sigmasq.epsilon = fit$sigma^2
	sigmasq.u = sigmasq.epsilon*exp(2*unlist(fit$modelStruct))

	Lambda = matrix(0, ncol(C), ncol(C))
	i1 = length(fit$coef$fixed)+1; i2 = ncol(C)
	Lambda[i1:i2,i1:i2] = sigmasq.epsilon * solve( diag(rep(sigmasq.u, block.inds)) ) 

	df = (diag( solve(t(C)%*%C + Lambda) %*% t(C)%*%C ))

return(df)
}




### Functions to perform bootstrap calculations of variance
###--------------------------------------------------------

bootstrap.effect <- function(fit, data, model, method, trend.input, nonlinear.trend.input=NULL, locations, cov.locations=NULL, cluster.input=NULL, cluster.locations=NULL, boot.samples=100, approximate=TRUE ,seeding.bootstrap=12345){
	if (!(method %in% c(6,7))) print("Bootstrap.effect function only implemented for methods 6 and 7.")

	else{
		parameters = extract.covariogram.parameters(fit)
		xboot = c(data[,2], locations[,1])
		yboot = c(data[,3], locations[,2])

		prediction.data1 = rep(0,nrow(locations))
		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations,cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		names(prediction.data) = names(data)
		all.data = rbind(data, prediction.data)
		Xmatrix = create.X.matrix(all.data, trend.input)
		
		output = list()
	
	if (class(fit)[1]=="lme"){
		fixed.par = fit$coefficients$fixed
		X.fit = create.X.matrix(data, trend.input)
		X.pred = create.X.matrix(prediction.data, trend.input)

		x.dist1 = abs(outer(data[,2], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist1 = abs(outer(data[,3], unique(data[,2:3])[fit$knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(data[,2:3])[fit$knots,1], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[fit$knots,2], unique(data[,2:3])[fit$knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)

		random.par = as.vector(unlist(fit$coefficients$random))
		if (is.null(nonlinear.trend.input)){
			if (is.null(cluster.input)){
				Z.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				block.inds = NULL				
				number.of.knots = ncol(Z.fit)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1.pred, Z2.pred)
				parm.cluster = parameters[5]
			}
		}

		if (!(is.null(nonlinear.trend.input))){
			if (is.null(cluster.input)){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				Z.pred = cbind(Z1.pred, Z2.pred)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z3.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit, Z3.fit)
				block.inds = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))
				number.of.knots = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				Z3.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				Z.pred = cbind(Z1.pred, Z2.pred, Z3.pred)
				parm.cluster = parameters[5]
			}
		}

		C.pred = cbind(X.pred,Z.pred)
		C.fit = cbind(X.fit, Z.fit)
		Z.all = rbind(Z.fit,Z.pred)

		output$Z.fit = Z.fit
		output$Z.pred = Z.pred
		output$C.fit = C.fit
		output$C.pred = C.pred

		set.seed(seeding.bootstrap)
		b = 1
		k = 1
		result.matrix = matrix(nrow=boot.samples, ncol=ncol(C.fit))
		while( b<(boot.samples+1) & k<(2*boot.samples) ){

			if (is.null(nonlinear.trend.input)){
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par
				}
				if (!(is.null(cluster.input))){
					Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
					random.par[(ncol(Z1.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z2.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}						
			}
			if (!(is.null(nonlinear.trend.input))){
				Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
				if (!(is.null(cluster.input))){
					random.par[(ncol(Z1.fit)+ncol(Z2.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z3.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
			}

			Sboot = simulate.Gaussian.Field(model, parameters, xboot, yboot)
			rand.numb.boot = rnorm(length(xboot), mean=0, sd=sqrt(parameters[2]))
			boot.resp = Sboot + mu + rand.numb.boot

			boot.sample = data.frame(cbind(boot.resp, all.data[,-1]))
			boot.data = boot.sample[1:nrow(data),]
			true.value = (Sboot + mu)[(nrow(data)+1):nrow(all.data)]

			if (approximate != "TRUE"){
				fit.boot = try( fit.bootstrap.lme.model(method, boot.data, X.fit, Z.fit, block.inds) )
			} 
			if (approximate == "TRUE"){
				fit.boot = try( bootstrap.lme.approx.model(fit, boot.data, X.fit, Z.fit, number.of.knots) )
			}

			if (class(fit.boot)!="try-error"){
				par.boot = c(fit.boot$coef$fixed,unlist(fit.boot$coef$random))
				result.matrix[b,] = par.boot
				print(paste("Bootstrap run ",b,sep=""))
				b=b+1
			}
			k=k+1		
		}
	} # end of bootstrap for spline method
	
output$pars = result.matrix
return(output)
	} # end of else statement
}



bootstrap.variance.spatial.effect <- function(fit, data, model, method, trend.input, nonlinear.trend.input=NULL, locations, cov.locations=NULL, cluster.input=NULL, cluster.locations=NULL, boot.samples=100, approximate=TRUE ,seeding.bootstrap=12345){
	if (!(method %in% c(6,7))) print("Bootstrap variance function only implemented for methods 1, 2, 6 and 7.")

	else{
		parameters = extract.covariogram.parameters(fit)
		xboot = c(data[,2], locations[,1])
		yboot = c(data[,3], locations[,2])

		prediction.data1 = rep(0,nrow(locations))
		if (is.null(cov.locations)){
			prediction.data2 = locations
		}
		if (!(is.null(cov.locations))){
			prediction.data2 = cbind(locations,cov.locations)
		}
		prediction.data = cbind(prediction.data1, prediction.data2)
		names(prediction.data) = names(data)
		all.data = rbind(data, prediction.data)
		Xmatrix = create.X.matrix(all.data, trend.input)
		Xmatrix[(nrow(data)+1):nrow(Xmatrix),]  = 0
		result.matrix = matrix(nrow=nrow(locations), ncol=boot.samples)

	if (class(fit)[1]=="lme"){
		fixed.par = fit$coefficients$fixed
		X.fit = create.X.matrix(data, trend.input)
		X.pred = create.X.matrix(prediction.data, trend.input)
		 X.pred = matrix(0, nrow(X.pred), ncol(X.pred))

		x.dist1 = abs(outer(data[,2], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist1 = abs(outer(data[,3], unique(data[,2:3])[fit$knots,2],"-"))
		dist.data = sqrt(x.dist1^2 + y.dist1^2)

		x.dist2 = abs(outer(unique(data[,2:3])[fit$knots,1], unique(data[,2:3])[fit$knots,1],"-"))
		y.dist2 = abs(outer(unique(data[,2:3])[fit$knots,2], unique(data[,2:3])[fit$knots,2],"-"))
		dist.knots = sqrt(x.dist2^2 + y.dist2^2)

		random.par = as.vector(unlist(fit$coefficients$random))
		if (is.null(nonlinear.trend.input)){
			if (is.null(cluster.input)){
				Z.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				block.inds = NULL				
				number.of.knots = ncol(Z.fit)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				 Z2.pred = matrix(0, nrow(Z2.pred), ncol(Z2.pred) ) 
				Z.pred = cbind(Z1.pred, Z2.pred)
				parm.cluster = parameters[5]
			}
		}

		if (!(is.null(nonlinear.trend.input))){
			if (is.null(cluster.input)){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z.fit = cbind(Z1.fit, Z2.fit)
				block.inds = c(ncol(Z1.fit) , ncol(Z2.fit))
				number.of.knots = c(ncol(Z1.fit) , ncol(Z2.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				 Z2.pred = matrix(0, nrow(Z2.pred), ncol(Z2.pred) ) 
				Z.pred = cbind(Z1.pred, Z2.pred)
			}
			if (!(is.null(cluster.input))){
				Z1.fit = create.Z.matrix(data, dist.data, dist.knots, model, parameters[3])
				Z2.fit = create.Z.matrix.nonlinear(data, nonlinear.trend.input)$matrix
				Z3.fit = create.Z.cluster.matrix(cluster.input)
				Z.fit = cbind(Z1.fit, Z2.fit, Z3.fit)
				block.inds = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))
				number.of.knots = c(ncol(Z1.fit), ncol(Z2.fit), ncol(Z3.fit))

				Z1.pred = create.Zpred.matrix(data, locations, model, parameters[3], fit$knots)
				Z2.pred = create.Zpred.matrix.nonlinear(data, prediction.data, nonlinear.trend.input)
				 Z2.pred = matrix(0, nrow(Z2.pred), ncol(Z2.pred) ) 
				Z3.pred = create.Zpred.cluster.matrix(cluster.input, cluster.locations)
				 Z3.pred = matrix(0, nrow(Z3.pred), ncol(Z3.pred) ) 
				Z.pred = cbind(Z1.pred, Z2.pred, Z3.pred)
				parm.cluster = parameters[5]
			}
		}

		C.pred = cbind(X.pred,Z.pred)
		C.fit = cbind(X.fit, Z.fit)
		Z.all = rbind(Z.fit,Z.pred)

		set.seed(seeding.bootstrap)
		b = 1
		k = 1
		while( b<(boot.samples+1) & k<(2*boot.samples) ){

			if (is.null(nonlinear.trend.input)){
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par
				}
				if (!(is.null(cluster.input))){
					Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
					random.par[(ncol(Z1.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z2.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}						
			}
			if (!(is.null(nonlinear.trend.input))){
				Zmatrix = Z.all[, (ncol(Z1.fit)+1) : ncol(Z.all) ]
				if (is.null(cluster.input)){
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
				if (!(is.null(cluster.input))){
					random.par[(ncol(Z1.fit)+ncol(Z2.fit)+1) : ncol(Z.all)] = rnorm(ncol(Z3.fit),0,sqrt(parameters[5]))
					if (length(fixed.par)==1) mu = Xmatrix * fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
					if (length(fixed.par)!=1) mu = Xmatrix %*% fixed.par + Zmatrix %*% random.par[(ncol(Z1.fit)+1) : ncol(Z.all)]
				}
			}

			Sboot = simulate.Gaussian.Field(model, parameters, xboot, yboot)
			rand.numb.boot = rnorm(length(xboot), mean=0, sd=sqrt(parameters[2]))
			boot.resp = Sboot + mu + rand.numb.boot

			boot.sample = data.frame(cbind(boot.resp, all.data[,-1]))
			boot.data = boot.sample[1:nrow(data),]
			true.value = (Sboot)[(nrow(data)+1):nrow(all.data)]

			if (approximate != "TRUE"){
				fit.boot = try( fit.bootstrap.lme.model(method, boot.data, X.fit, Z.fit, block.inds) )
			} 
			if (approximate == "TRUE"){
				fit.boot = try( bootstrap.lme.approx.model(fit, boot.data, X.fit, Z.fit, number.of.knots) )
			}

			if (class(fit.boot)!="try-error"){
				par.boot = c(fit.boot$coef$fixed,unlist(fit.boot$coef$random))
				result.matrix[,b] = (C.pred%*%par.boot - true.value)^2
				print(paste("Bootstrap run ",b,sep=""))
				b=b+1
			}
			k=k+1		
		}
	} # end of bootstrap for spline method
	
output = apply(result.matrix,1,mean)
return = output
return(return)
	} # end of else statement
}

