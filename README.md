## Introduction and Dependencies

We define several functions to do a kernel-based ridge regression, and subsequent cross-validation. This is a technique for general non-linear regression using kernels. A complete treatment of the theory behind this can be found in Hastie, Tibshirani and Friedman (2016), specifically section 5.8. A linear kernel can be used for a generalized kernel-based linear regression. The only dependency is :

* `kernlab` available in
[CRAN](https://cran.r-project.org/web/packages/kernlab/index.html).

Please consult the `kernlab`
[manual](https://cran.r-project.org/web/packages/kernlab/kernlab.pdf) to become
familiar with this package, and available kernels it provides.

## Usage

Suppose we have a dataset with n observations, k descriptors and 1 response
variable. Then we have an n x 1 response vector, `Y` and an n x k descriptor
matrix, `X`. It is recommended that `X` be normalized by column (by descriptor).

For a given kernel, we refer to the n x n kernel matrix as `G`.

### `klm`

`klm(G, y, lambda)`

* `G` - the n x n kernel matrix
* `y` - the n x 1 response vector
* `lambda` - a meta-parameter, scale factor. Values too large will result in
robust models at the expense of predictive performance. Values too small will
result in overfitting (or error in computation if extremely small).

##### Output of `klm`

Outputs a list consisting of

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

Outputs a vector consisting of the predicted values of the test set.


## Sources

Hastie Trevor, Tibshirani Robert and Friedman Jerome. 2016. *The Elements of
Statistical Learning: Data mining, Inference and Prediction, Second Edition
(Springer Series in Statistics)*. Springer.
