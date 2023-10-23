function R = weighted_partialcorrs(X, Y, w)
%% Description 
% Calculates the exponentially weighted correlation matrix R of a matrix X,
% given Y
%% Inputs
% X: a dt-by-N matrix of data points (e.g. returns)
% Y: a dt-by-1 factor vector 
% w: a 1-by-dt vector of weights that sum to one
%% Ouputs 
% R: the N-by-N correlation matrix. 
%% Setup 
dt = size(X, 1);

%check that Y is dt-by-1

if ~isequal(size(Y), [dt, 1])
    error("Y is not dt-by-1")
end

%% Calculation
%calculate residuals
residuals = single_factor_residuals(X, Y);
%calculate correlations 
R = weightedcorrs(residuals, w);
end