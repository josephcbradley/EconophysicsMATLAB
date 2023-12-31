function [rho, n_significant_corrs] = calculate_rolling_mean_corr(R, dt, t_range, options)
%% Description 
% Calcualtes the rolling mean of the upper triangle of the correlation matrix
% C, calculated from rolling windows over X

%% Inputs
% R - the T-by-N matrix of returns 
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
        dt (1, 1) double
        t_range (1, :) double
        options.CorrHandle (1, 1) function_handle = @(r) weighted_pearson_corrs(r, generate_expweights(dt, dt/3))
        options.RemoveInsignificant (1, :) logical = false
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
    C = options.CorrHandle(X);
    C = force_symmetric(C); %bc MATLAB sometimes doesn't!
    pvals = bootstrap_correlation_signifiance(X, C, 100, options.CorrHandle);
    %calculate significant corrs 
    significant = pvals < 0.05;
    n_significant_corrs(k) = sum(utri_to_vec(significant));
    if options.RemoveInsignificant %remove anything not significant 
        C(pvals >= 0.05) = 0;
        C = nearest_pos_semi_definite(C);
        if ~isreal(C)
            warning("C is not real!")
        end
    end
    rho(k) = mean(utri_to_vec(C), 'omitnan'); %mean of upper triangle
end
end