function [rho, n_significant_corrs] = calculate_rolling_mean_corr(R, w, dt, t_range)
%% Description 
% Calcualtes the rolling mean of the upper triangle of the correlation matrix
% C, calculated from rolling windows over X

%% Inputs
% R - the T-by-N matrix of returns 
% w - the 1-by-dt vector of weights
% dt - the length of the window over which correlations are measured
% t_range - the vector of time points t at which correlations will be
% measured

%% Ouputs 
% rho - the n_windows-by-1 vector of the mean of the correlation 
% coefficients

%% Setup 
%check that w and window_size match 
if ~isequal(size(w, 2), dt)
    error("The weights vector is not equal to the window size")
end

%check that the first value of t_range is at least dt
if t_range(1) < dt 
    error('The first values of t is less than dt.')
end

%check that the last value of t is within the size of R
if t_range(end) > size(R, 1)
    error('The final value for t is too large given R.')
end

%allocate memory
rho = NaN(length(t_range), 1);
n_significant_corrs = NaN(length(t_range), 1);
%% Calculation
for k = 1:length(t_range) %index for windows 
    t = t_range(k);
    local_tspan = t-dt+1:t;
    %get subset of data used for correlation coefficients 
    X = R(local_tspan, :);
    %correlation coefficients
    C = weightedcorrs(X, w);
    %remove anything not significant 
    pvals = bootstrap_correlation_signifiance(X, C, 100, @(x) weightedcorrs(x, w));
    %calculate significant corrs 
    significant = pvals < 0.05;
    n_significant_corrs(k) = sum(utri_to_vec(significant));
    C(pvals >= 0.05) = NaN;
    rho(k) = mean(utri_to_vec(C), 'omitnan'); %mean of upper triangle
end
end