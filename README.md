# ABIP For Machine Learning Problems

ABIP(`ADMM-Based Interior Point Method`) is a first order method for convex optimization [1]. In this project, we develop a customized ABIP solver for large scale LP in machine learning. Our main work includes

1.  We present a friendly modeling interface, which gives plenty of flexibility for customizing ABIP to different tasks.

2. We developed customized linear systems solvers for a variety of machine learning applications. Our customized solvers exploit the problem structure very efficiently and significantly improve the performance of ABIP on those applications. 

It is written in MATLAB and perform well on some machine learning problems.
ABIP-ML is developed by Wenhao Gu (https://github.com/WenhaoGu2000) while he was visiting RIIS, SHUFE. Contact: g1751246188@gmail.com

## Examples
This package uses the ABIP algorithm to solve linear programs. We provide examples of general linear program and two machine learning problems, which can be solved efficiently due to the special block structure of coefficient matrices.

Before using it, add `abip` to you MATLAB path. Here we use string variable `ABIP_HOME` to denote the path where you download the package.
```
addpath(ABIP_HOME+'/abip-machine-learning/abip')
```
### General LP
We solve general linear programs in the following standard form:

<img src="https://latex.codecogs.com/svg.latex?\min&space;c^Tx&space;\quad&space;s.t.&space;\quad&space;Ax=b,&space;\quad&space;x\geq&space;0" title="\min c^Tx \quad s.t. \quad Ax=b, \quad x\geq 0" />

We provide a test dataset `ABIP_HOME/datasets/ADLITTLE.mat`. The following demo shows how to solve general LP by ABIP.
```
% load dataset
load('./datasets/ADLITTLE.mat');
data = struct('A', Problem.A, 'b', Problem.b, 'c', Problem.c);
% set params
params = struct();
% solve
[x,y,s,info] = abip(data,params);
```
Here the output `x` is the primal optimal solution and `y,s` is the dual optimal solution. 

### <img src="https://latex.codecogs.com/svg.latex?\large&space;\ell_1" title="\large \ell_1" /> SVM
<img src="https://latex.codecogs.com/svg.latex?\ell_1" title="\ell_1" /> SVM is used in classification problems. Given training data `X` and label `y` of -1 and 1, we aim to find the supporting hyperplane, i.e., to solve the following problem [2]

<img src="https://latex.codecogs.com/svg.latex?\min_{\beta_0,&space;\beta}&space;\quad&space;\sum_{i=1}^{n}[1-y_i(\beta_0&plus;x_i\beta)]_&plus;&space;&plus;&space;\lambda&space;\Vert&space;\beta&space;\Vert_1" title="\min_{\beta_0, \beta} \quad \sum_{i=1}^{n}[1-y_i(\beta_0+x_i\beta)]_+ + \lambda \Vert \beta \Vert_1" />

We provide a test dataset `ABIP_HOME/datasets/australian.mat`. You can use the following MATLAB code.

```
% load dataset
load('./datasets/australian.mat');
data = struct('X', Problem.X, 'y', Problem.y);
data.scalar = 0.1; % (Optional) Tuning parameter \lambda, with default value 1.0
% set params
params = struct('Problem', 'L1_SVM'); % specify the abstract class it used
% solve
[x,y,s,info] = abip(data,params);
```


### Dantzig Selector
The Dantzig Selector is an estimator for linear regression. Given the data `X` and response `y` , we aim to fit the model `y=X*r+z` where `z` is i.d.d. noise vector.  The problem is formulated as follows [3]

<img src="https://latex.codecogs.com/svg.latex?\min_{\beta}&space;\quad&space;\Vert&space;\beta&space;\Vert_1&space;\quad&space;s.t.&space;\quad&space;\Vert&space;X^T(y-X\beta)&space;\Vert_{\infty}&space;\leq&space;\lambda" title="\min_{\beta} \quad \Vert \beta \Vert_1 \quad s.t. \quad \Vert X^T(y-X\beta) \Vert_{\infty} \leq \lambda" />

We provide a test dataset `ABIP_HOME/datasets/bodyfat.mat`

```
% load dataset
load('./datasets/bodyfat.mat');
data = struct('X', Problem.X, 'y', Problem.y);
data.scalar = 0.1; % (Optional) Tuning parameter \lambda, with default value 1.0
% set params
params = struct('Problem', 'Dantzig_selector'); % specify the abstract class it used
% solve
[x,y,s,info] = abip(data,params);
```

## Customized Linear System Solver
According to ABIP, in each inner iteration, we need to solve linear systems with the same coefficient matrices (ABIP is matrix-free), so peculiar block-structure of the coefficient matrix can be used. Moreover, for some problems (e.g. Dantzig selector), there's no need to work out the matrix `A` explicitly. So we provide abstract class APIs for efficiently solving problems with different structures. The abstract classes for general LP, <img src="https://latex.codecogs.com/svg.latex?\ell_1" title="\ell_1" /> SVM, and Dantzig selector are in the folder `ABIP_HOME/abip/+model`. Moreover, block informations of the coefficient matrix, tuning parameters and residuals are stored in structs  `work`, `settings`, `barrier` and  `residual`  

`params_reset(obj, params)` This function is used to set defalut values and tuning parameters.

`A_times(obj, x)` Compute `A*x`. Note that matrix `A` does not have an explicit form.

`AT_times(obj, x)` Compute `A'*x`.

`normalize_data(obj, work, setting)` Data normalization. Some special normalization is applied when problems are different.

`Q_times(obj, x)` Compute `Q'*x`, where `Q` equals to `[O, A, -b; -A^T, O, c; b^T, -c^T, 0]`

`solve_lin_sys(obj, work, rhs, setting, init, warm_start, k)` Solve linear system `M*x=rhs`. The matrix `M` equals to `[\rho*I, A; A^T, -I]`. `warm_start` can be used if we use conjugate gradient method(CG).

`no_normalize_settings(obj, work, setting)` This function is used if normalization is disabled and additional informations should be added. (Usually this function is empty)

`individual_sparsity_ratio(obj)` Compute the (approximated) sparisity ratio of matrix `A`. It can also be set to constant if the real sparsity is difficult to calculate or unimportant.

`individual_data_work_settings(obj,work,setting,part)` Add informations to `data` and `work`, such as the LDL factorization of `I+A*A'`.

## Bibliography
[[1](https://arxiv.org/abs/1805.12344)] Tianyi Lin, Shiqian Ma, Yinyu Ye & Shuzhong Zhang(2020) *An ADMM-based interior point method for large-scale linear programming*, Optimization Methods and Software.  
[[2](https://papers.nips.cc/paper/2003/hash/49d4b2faeb4b7b9e745775793141e2b2-Abstract.html)] J. Zhu, S. Rosset, R. Tibshirani, and T.J. Hastie. *1-norm support vector machines*. In *Advances in neural information processing systems (NIPS)*, pages 49-56, 2004.  
[[3](https://arxiv.org/abs/math/0506081)] E. Candes and T. Tao. *The dantzig selector: Statistical estimation when p is much larger than n*, Ann. Statist. **35**(12, 2007)2313-2351.  
