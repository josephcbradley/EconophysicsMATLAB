function [rho, n_significant_corrs] = calculate_rolling_mean_partial_corr(R, Y, dt, t_range, options)
%% Description 
% Calcualtes the rolling mean of the upper triangle of the correlation matrix
% C, calculated from rolling windows over X

%% Inputs
% R - the T-by-N matrix of returns 
% Y - the T-by-1 vector of factor returns
% dt - the length of the window over which correlations are measured
% t_range - the vector of time points t at which correlations will be
% measured

%% Ouputs 
% rho - the n_windows-by-1 vector of the mean of the correlation 
% coefficients
% n_significant_corrs - the number of coefficients in the upper triangle of
% the correlation matrix that are significant at the 5% level

%% Setup 

arguments
        R (:, :) double
        Y (:, 1) double
        dt (1, 1) double
        t_range (1, :) double
        options.CorrHandle (1, 1) function_handle = @(r, y) weighted_partialcorrs(r, y, generate_expweights(dt, dt / 3))
        options.RemoveInsignificant (1, 1) logical = true
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
    C = options.CorrHandle(X, Y(local_tspan));
    %remove anything not significant 
    if options.RemoveInsignificant
        pvals = bootstrap_correlation_signifiance(X, C, 100, @(r) options.CorrHandle(r, Y(local_tspan)));
        %calculate significant corrs 
        significant = pvals < 0.05;
        n_significant_corrs(k) = sum(utri_to_vec(significant));
        C(pvals >= 0.05) = NaN;
    end
    rho(k) = mean(utri_to_vec(C), 'omitnan'); %mean of upper triangle
end
end