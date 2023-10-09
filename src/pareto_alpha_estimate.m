function alpha = pareto_alpha_estimate(X)
%% Description 
% Estimates the scaling parameter alpha in the normal Pareto distribution
% using maximum likelihood estimation (MLE)
% MLE formula is \hat{\alpha} = \frac{N}{\sum_{i=1}^N \ln \frac{x_i}{x_{\min }}}
% The code will skip NaNs in X

%% Inputs
% X - the vector of datapoints. Should be one dimensional.

%% Ouputs 
% alpha - the estimator for the scaling exponent alpha

%% Setup 
% NA

%% Calculation
N = numel(X);
xmin = min(X, [], 'all');
alpha = N / (sum(log(X ./ xmin)));

end