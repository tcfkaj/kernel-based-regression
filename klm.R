source(dependencies.R)

## Kernel-based regression
klm <- function(X, y, lamba=0.1, kern=rbf(sigma=0.01)){
		G <- kernelMatrix(X, kernel=kern)
		z <- solve(G + lambda * diag(length(y)), y)
		b <- -mean(G %*% z);
		pred <- G %*% z + b;
		err <- pred - y;
		mse <- mean(err^2);
		output <- list(z=z, b=b,pred=pred,
			       kern=kern, err=err,
			       mse=mse, y=y, X=X,
			       G=G, lambda=lambda);
		return(output)
}

## Predict
predict.klm <- function(model, new_data){
		X <- rbind(model$X, new_data)
		end_here = nrow(model$X)
		G <- kernelMatrix(X, kernel=model$kern)
		G_pred <- G[(end_here+1):ncol(X),1:end_here]
		return(G_pred %*% model$z + model$b)
}

