function pvals = bootstrap_correlation_signifiance(X, real_C, nreps, C_handle)
%% Description 
%calculate the p values for the correlations in X. 
%uses a boostrap technique but repeatedly shuffling the rows in X
%returns things as a matrix - quite useful

%% Inputs
% X - the T-by-N data matrix (e.g. of returns)
% real_C - the N-by-N matrix of calculate correlation coefficients 
% nreps - the number of bootstraps to calculate 
% C_handle - a function handle that accepts and array of the same dimensions
% as X and calculates 

%% Ouputs 
% pvals - a N-by-N matrix of p-values for each value of C(i, j)
%% Setup 
if ~isempty(find(isnan(X), 1))
    error("NaNs in X!")
end

N = size(X, 2);

if size(real_C, 1) ~= N
    error("Sizes of X and C do not match")
end

%% Calculation
%for each r, want to track how many times abs(r_true) > abs(r_rand)
stats = zeros(size(real_C));
for i = 1:nreps
    %get a random perm 
    shuffled_X = shuffle_matrix_along_columns(X);
    shuffled_C = C_handle(shuffled_X);
    %update_stats
    stats = stats + (abs(real_C) > abs(shuffled_C));
end
%stats are number of times real C is greater than expected 
stats = stats ./ nreps;

%pval is the prob. that shuffled r is >= real R. which is 1-frac. times
%that real R > shuffled r. So 1 - stats
pvals = 1 - stats;
end
   