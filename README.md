## Introduction and Dependencies

We define several functions to do a kernel-based ridge regression, and subsequent cross-validation. This is a technique for general non-linear regression using kernels. A complete treatment of the theory behind this can be found in Hastie, Tibshirani and Friedman (2016), specifically section 5.8. A linear kernel can be used for a generalized kernel-based linear regression. However, even non-linear kernels map features to a linear function in Hilbert space, hence the name of the function `klm`. The only **dependency** is :

* `kernlab` available in
[CRAN](https://cran.r-project.org/web/packages/kernlab/index.html).

Please consult the `kernlab`
[manual](https://cran.r-project.org/web/packages/kernlab/kernlab.pdf) to become
familiar with this package, and available kernels it provides.

## Usage

Suppose we have a dataset with n observations, k descriptors and 1 response
variable. Then we have an n x 1 response vector, `Y` and an n x k descriptor
matrix, `X`. It is recommended that `X` be normalized by column (by descriptor) and `Y` be shifted by its mean.

For a given kernel, we refer to the n x n kernel matrix as `G`.

### `klm`

`klm(G, y, lambda)`

* `G` - the n x n kernel matrix
* `y` - the n x 1 response vector
* `lambda` - a meta-parameter, scale factor. Values too large will result in
robust models at the expense of predictive performance. Values too small will
result in overfitting (or error in computation if extremely small).

##### Output of `klm`

Outputs a named list, which can be referenced with usual `$`, containing:

* `$z` - the solution matrix
* `$b` - a constant added to each variable, often referred to to b<sub>0</sub> in standard linear models. It is calculated as the negative mean of `G` x `z`.
* `$pred` - the vector of predicted response values on the training data.
* `$err` - the vector of errors of the predicted values.
* `$mse` - the mean squared error.
* `$G` - returns original `G`.

### `predict.klm`

`predict.klm(G, test_samp, train_model)`

* `G` - the n x n kernel matrix
* `test_samp` - the indices of the columns of `G` that will act as the test
set.
* `train_model` - trained `klm`, trained on the training set.

##### Output of `predict.klm`

Outputs a vector consisting of the predicted values of the test set. One random
partition is used and all partitions are exhausted.

### `cv.klm`

Perform cross-validation on `klm`.

`cv.klm(G, y, lambda, cv=5)`


* `G` - the n x n kernel matrix
* `y` - the n x 1 response vector
* `lambda` - a meta-parameter, scale factor. Values too large will result in
* `cv` - number of cross validation "folds". For "leave one out" set to `'LOOCV'`
or the size of data set.

#### Output of `cv.klm`


Outputs a named list, which can be referenced with usual `$`, containing:

* `$train_mse` - list of training set mean squared errors, for each "fold".
* `$test_mse` - list of test set mean squared errors, for each "fold".
* `$train_rse` - square root of `train_mse`.
* `$test_rse` - square root of `test_mse`.
* `$preds` - a list containing `cv` number of lists, each the predicted values
for each fold.
* `$test_sd` - standard deviation of `test_rse`.


### `mccv.klm`

Perform monte-carlo cross-validation on `klm`. A new random partition is
generated with each fold. Observations can be used multiple times, and some
observations may not be used.

`mccv.klm(G, y, lambda, cv=5)`


* `G` - the n x n kernel matrix
* `y` - the n x 1 response vector
* `lambda` - a meta-parameter, scale factor. Values too large will result in
* `cv` - number of cross validation "folds". For "leave one out" set to `'LOOCV'`
or the size of data set.

#### Output of `mccv.klm`

Outputs a named list, which can be referenced with usual `$`, containing:

* `$train_mse` - list of training set mean squared errors, for each "fold".
* `$test_mse` - list of test set mean squared errors, for each "fold".
* `$train_rse` - square root of `train_mse`.
* `$test_rse` - square root of `test_mse`.
* `$preds` - a list containing `cv` number of lists, each the predicted values
for each fold.
* `$test_sd` - standard deviation of `test_rse`.

### `find_best_gamma`

Iterates `klm` method over possible values of `gamma` to find the best
performing one, based on test-set performance with `cv` or `mccv`.

`find_best_gamma(X, y, gammas, cv=5, lambda=1 , scale=FALSE, scale_lam=0.1, cv_type="cv")`

* `X` - the original predictor matrix (not kernel matrix).
* `y` - the n x 1 response vector.
* `gammas` - list of values of gamma to be iterated over.
* `cv` - number of cross validation "folds" for each iteration. For "leave one out" set to 'LOOCV'
or the size of data set.
* `lambda` - global value to be used for all iterations, if not `scale`.
* `scale` - boolean, if TRUE (default), a new `lambda` is choosen during each
iteration, as a percentage of the largest eigenvalue of the kernel matrix.
* `scale_lam` - percentage of largest eigenvalue to be used as `lambda` on each
iteration.
* `cv_type` - `'cv'` (default) or `'mccv'`.

#### Output of `find_best_gamma`

Outputs a named list, which can be referenced with usual `$`, containing:

* `$train` - list of the mean training mse for each `gamma`.
* `$test` - list of the mean test mse for each `gamma`.
* `$best_mse` - test set mse for best `gamma`.
* `$best_sd` - mean standard deviation of best `gamma`.
* `$best_gamma` - best `gamma`, with minimum test set prediction error.
* `$best_rse` - mean test set `rse` for best `gamma`.
* `$train_rse` - list of mean training set `rse` for each `gamma`.
* `$test_rse` - list of mean test set `rse` for each `gamma`.
* `$test_sds` - list of mean test set `sd` for each `gamma`.


### `find_best_lambda`

Iterates `klm` method over possible values of `lambda`, while fixing `gamma`, to find the best
performing one, based on test-set performance with `cv` or `mccv`.

`find_best_lambda(X, y, lambdas, cv=5, gamma=1 , cv_type="cv")`

* `X` - the original predictor matrix (not kernel matrix).
* `y` - the n x 1 response vector.
* `lambdas` - list of values of `lambda` to be iterated over.
* `cv` - number of cross validation "folds" for each iteration. For "leave one out" set to 'LOOCV'
or the size of data set.
* `gamma` - global value to fix `gamma` at.
* `cv_type` - `'cv'` (default) or `'mccv'`.

#### Output of `find_best_lambda`

Outputs a named list, which can be referenced with usual `$`, containing:


* `$train` - list of the mean training mse for each `lambda`.
* `$test` - list of the mean test mse for each `lambda`.
* `$best_mse` - test set mse for best `lambda`.
* `$best_sd` - mean standard deviation of best `lambda`.
* `$best_lambda` - best `lambda`, with minimum test set prediction error.
* `$best_rse` - mean test set `rse` for best `lambda`.
* `$train_rse` - list of mean training set `rse` for each `lambda`.
* `$test_rse` - list of mean test set `rse` for each `lambda`.
* `$test_sds` - list of mean test set `sd` for each `lambda`.

### Some example code

```R
library(ISwR)

X <- scale(bp.obese$obese)
y <- bp.obese$bp - mean(bp.obese$bp)

## Save ridge.R in pwd and include with source("filename")
source("ridge.R")

gammas <- c(0.00001, 0.0001, 0.001,
	    0.005, 0.01, 0.05,
	    0.1, 0.5, 2,
	    5, 10, 50, 100, 1000)
lambda=500
ptm <- proc.time()
best <- find_best_gamma(X, y, gammas, cv=20, lambda=lambda)

## Vector containing mean training-set mse's for each gamma
print("$train: ")
best$train


## Vector containing mean test-set mse's for each gamma
print("$test: ")
best$test

## Vector containing mean test mse's for each gamma
print("$best_mse: ")
best$best_mse


## Vector containing mean test rse=sqrt(mse) for each gamma
print("$test_rse: ")
best$rse


## To get list of normalized test performance errors
print("Normalized performance errors: ")
best$test_rse/mean(bp.obese$bp)


## Mean test-set rse for best gamma
print("$best_rse: ")
best$best_rse

## To get performance of best gamma
print("Best normalized performance error: ")
best$best_rse/mean(bp.obese$bp)


## Vector containing mean sd for each gamma
print("$test_sd: ")
best$test_sd

## Sd of best gamma
print("$best_sd: ")
best$best_sd

## Best gamma
print("Best gamma: ")
best$best_gamma

proc.time() - ptm
```

## Sources

Hastie Trevor, Tibshirani Robert and Friedman Jerome. 2016. *The Elements of
Statistical Learning: Data mining, Inference and Prediction, Second Edition
(Springer Series in Statistics)*. Springer.
