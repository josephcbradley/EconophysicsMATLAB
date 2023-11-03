function C = weighted_pearson_corrs(X, w)
%% Description 
% Calculates the exponentially weighted correlation matrix R of a matrix X
%% Inputs
% X: a dt-by-N matrix of data points (e.g. returns)
% w: a 1-by-dt vector of weights that sum to one
%% Ouputs 
% R: the N-by-N correlation matrix. 
%% Setup 


arguments
    X (:, :) double
    w (:, :) double
end


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
C = (X_less_mean' .* w) * X;
%normalize
D = diag(C) + eps; %add tiny quantity eps to ensure things are not singular
C = C ./ sqrt(D * D');
%make symmetric 
C = (C + C') / 2;
%guarantee pos_semi_def 
C = nearest_pos_semi_definite(C);
end