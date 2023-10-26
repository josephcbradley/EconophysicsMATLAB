function tau = weighted_kendall_corrs(X, w, i2, i1)
%% Description 
% Calculates the exponentially weighted correlation matrix R of a matrix X
%% Inputs
% X: a dt-by-N matrix of data points (e.g. returns)
% w: a 1-by-dt vector of weights that sum to one
%% Ouputs 
% C: the N-by-N correlation matrix. 
%% Setup 
%check that weights has the correct length 
% if ~isequal(size(w), [1, dt])
%     error("The length of the weights vector does not match the length of X.")
% end
% %check that X contains no missing data 
% if ~isempty(find(isnan(X), 1))
%     error("There are NaN values in X!")
% end
% %check that weights contains no missing data
% if ~isempty(find(isnan(w), 1))
%     error("There are NaN values in weights!")
% end
arguments
        X (:, :) double
        w (:, :) double
        i2
        i1
end
[~, N] = size(X);

%% Calculation
% Indexes for all dt * (dt - 1) / 2 combinations without repetition: % {1, 2}, {1, 3}, ..., {1, dt}, {2, 3}, {2, 4}, ..., {2, dt}, ..., {dt - 1, dt} 
%[i2, i1] = find(tril(ones(dt, 'uint8'), -1)); 
% Signs of differences between variables at time i2 and at time i1 
tau = sign(X(i2, :) - X(i1, :)); 
% Number of concordant/discordant pairs (weighted) 
tau = tau' * (tau .* repmat(w, 1, N)); 
% Must be exactly symmetric 
tau = 0.5 * (tau + tau'); 
% Number of pairs minus number of ties (weighted) 
temp = diag(tau); 
% Matrix of Weighted Correlation Coefficients 
tau = tau ./ sqrt(temp * temp');
end
