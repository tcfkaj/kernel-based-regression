library(kernlab, quietly=TRUE, warn.conflicts=FALSE)


## Kernel-based ridge regression

klm <- function(G, y, lambda){
		z <- solve(G + lambda * diag(length(y)), y)
		b <- -mean(G %*% z);
		pred <- G %*% z + b;
		err <- pred - y;
		mse <- mean(err^2);
		output <- list(z=z, b=b, pred=pred,
				err=err, mse=mse, G=G);
		return(output)
}


## Predict test set with:
## G - full kernel matrix
## test_samp - indices of test subset
## train_model - model trained on train subset

predict.klm <- function(G, test_samp, train_model){
			G_test <- G[test_samp,-test_samp];
			return(G_test %*% train_model$z + train_model$b)
}


## Regular n-fold cross validation
## cv="LOOCV" for Leave-one-out

cv.klm <- function(G, y, lambda, cv=5){
		train_mse = test_mse = preds =  c();
		samp_size <- ifelse(cv=="LOOCV",1,floor(length(y)/cv));
		main_samp <- sample(c(1:length(y)));
		for (i in 1:ifelse(cv=="LOOCV",length(y), cv)){
			test_samp <- main_samp[((i-1)*samp_size+1):
						ifelse(i==cv,length(y),
						i*samp_size)];
			G_train <- G[-test_samp, -test_samp];
			y_train <- y[-test_samp];
			model <- klm(G_train, y_train, lambda);
			pred <- predict.klm(G, test_samp, model);
			err <- pred - y[test_samp];
			mse <- mean(err^2);
			train_mse <- c(train_mse, model$mse);
			test_mse <- c(test_mse, mse);
			preds <- c(preds, pred);
			}
		train_rse <- sqrt(train_mse)
		test_rse <- sqrt(test_mse)
		test_sd <- sd(test_rse);
		output <- list(train_mse=train_mse, test_mse=test_mse,
				test_rse=test_rse, preds=preds,
				train_rse=train_rse, test_sd=test_sd);
		return(output)
}


## Monte carlo cross-validation

mccv.klm <- function(G, y, lambda, cv=5){
		train_mse = test_mse = preds =  c();
		samp_size <- floor(length(y)/cv);
		for (i in 1:cv){
			test_samp <- sample(seq(1:length(y)),
						 size=samp_size);
			G_train <- G[-test_samp, -test_samp];
			y_train <- y[-test_samp];
			model <- klm(G_train, y_train, lambda);
			pred <- predict.klm(G, test_samp, model);
			err <- pred - y[test_samp];
			mse <- mean(err^2);
			train_mse <- c(train_mse, model$mse);
			test_mse <- c(test_mse, mse);
			preds <- c(preds, pred);
			}
		train_rse <- sqrt(train_mse)
		test_rse <- sqrt(test_mse)
		test_sd <- sd(test_rse);
		output <- list(train_mse=train_mse, test_mse=test_mse,
				test_rse=test_rse, test_sd=test_sd,
				train_rse=train_rse, preds=preds);
		return(output)
}



## Find best gamma with fixed lambda

find_best_gamma <- function(X, y, gammas, cv=5, lambda=1 ,
				scale=TRUE, scale_lam=0.1, cv_type="cv"){
		train = test = test_sds = train_rse = test_rse = c();
		for (i in gammas){
			print(paste("Gamma: ",i));
			G <- kernelMatrix(X, kernel=rbfdot(sigma=i));
			G_eig <- eigen(G)$values
			lambda <- ifelse(!scale, lambda,
					scale_lam * G_eig[1]);
			if(cv_type=="mccv"){
				model <- mccv.klm(G, y, lambda, cv=cv)}
			else{
				model <- cv.klm(G, y, lambda, cv=cv)}
			train <- c(train, mean(model$train_mse))
			test <- c(test, mean(model$test_mse));
			m.test_rse <- mean(model$test_rse)
			m.train_rse <- mean(model$train_rse)
			test_rse <- c(test_rse, m.test_rse)
			train_rse <- c(train_rse, m.train_rse)
			test_sds <- c(test_sds, model$test_sd)
			}
		best_mse <- min(test)
		best_rse <- min(test_rse)
		best_sd <- test_sds[which(test==best_mse)]
		best_gamma <- gammas[which(test==best_mse)]
		output <- list(train=train, test=test, best_mse=best_mse,
				best_sd=best_sd, best_gamma=best_gamma,
				best_rse=best_rse, train_rse=train_rse,
				test_rse=test_rse, test_sd=test_sds);
		return(output)
}

find_best_lambda <- function(X, y, lambdas, cv=5, gamma=1 , cv_type="cv"){
		train = test = test_sds = train_rse = test_rse = c();
		G <- kernelMatrix(X, kernel=rbfdot(sigma=gamma));
		for (i in lambdas){
			print(paste("Lambda: ",i));
			if(cv_type=="mccv"){
				model <- mccv.klm(G, y, lambda, cv=cv)}
			else{
				model <- cv.klm(G, y, lambda, cv=cv)}
			train <- c(train, mean(model$train_mse))
			test <- c(test, mean(model$test_mse));
			m.test_rse <- mean(model$test_rse)
			m.train_rse <- mean(model$train_rse)
			test_rse <- c(test_rse, m.test_rse)
			train_rse <- c(train_rse, m.train_rse)
			test_sds <- c(test_sds, model$test_sd)
			}
		best_mse <- min(test)
		best_rse <- min(test_rse)
		best_sd <- test_sds[which(test==best_mse)]
		best_lambda <- lambdas[which(test==best_mse)]
		output <- list(train=train, test=test, best_mse=best_mse,
				best_sd=best_sd, best_lambda=best_lambda,
				best_rse=best_rse, train_rse=train_rse,
				test_rse=test_rse, test_sd=test_sds);
		return(output)
}

## Example code
## Leave commented
## Use this file with source("filename")
## X is scaled predictor data
## y is shifted response


## source("ridge.R")

## gammas <- c(0.00001, 0.0001, 0.001,
## 		0.005, 0.01, 0.05,
## 		0.1, 0.5, 2,
## 		5, 10, 50, 100, 1000)
## lambda=500

## ptm <- proc.time()
## best <- find_best_gamma(X, y, gammas, cv="LOOCV", lambda=lambda)

## Vector containing mean training mse's for each gamma
## best$train


## Vector containing mean training mse's for each gamma
## best$test

## Vector containing mean test mse's for each gamma
## best$best_mse


## Vector containing mean test rse=sqrt(mse) for each gamma
## best$rse


## To get list of test performances
## best$rse/mean(UnshiftedResponseData)


## Vector containing mean test rse for best gamma
## best$best_rse

## To get performance of best gamma
## best$best_rse/mean(UnshiftedResponseData)


## Vector containing mean sd for each gamma
## best$test_sd

## Sd of best gamma
## best$best_sd

## Sd of best gamma
## best$best_gamma

## proc.time() - ptm
##
## save(gammas, lambda, best, file="bestgamma-L500.rda")
