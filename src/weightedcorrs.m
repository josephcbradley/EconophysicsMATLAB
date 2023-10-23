function R = weightedcorrs(X, w)
%% Description 
% Calculates the exponentially weighted correlation matrix R of a matrix X
%% Inputs
% X: a dt-by-N matrix of data points (e.g. returns)
% w: a 1-by-dt vector of weights that sum to one
%% Ouputs 
% R: the N-by-N correlation matrix. 
%% Setup 
dt = size(X, 1);
%check that weights has the correct length 
if ~isequal(size(w), [1, dt])
    error("The length of the weights vector does not match the length of X.")
end
%check that X contains no missing data 
if ~isempty(find(isnan(X), 1))
    error("There are NaN values in X!")
end
%check that weights contains no missing data
if ~isempty(find(isnan(w), 1))
    error("There are NaN values in weights!")
end

%% Calculation
%create weighted mean of each column 
column_means = w * X;
%subtract means
X_less_mean = X - column_means;
%covariance calculation 
%lhs is X' with each row weighted
R = (X_less_mean' .* w) * X;
%normalize
D = diag(diag(R));
Dinvsqrt = inv(sqrt(D));
R = Dinvsqrt * R * Dinvsqrt;
%symmetry 
R = 0.5 * (R + R');
%check pos def
[~, p] = chol(R);
if p ~= 0
    eigs = eig(R);
    min_eig = min(eigs(~isnan(eigs)));
    warning(compose("R is not positive semi-definite - smallest eigenvalue is %.6f", min_eig));
end
end