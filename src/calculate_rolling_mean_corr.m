function rho = calculate_rolling_mean_corr(R, w, n_windows, dt_step, window_size)
%% Description 
% Calcualtes the rolling mean of the upper triangle of the correlation matrix
% C, calculated from rolling windows over X

%% Inputs
% R - the T-by-N matrix of returns 
% w - the 1-by-window_size
% n_windows - the scalar representing the number of windows over which to 
% calculate the correlations
% dt_step - the time between windows 
% window_size - the size of each window

%% Ouputs 
% mean_corrs - the n_windows-by-1 vector of the mean of the correlation 
% coefficients

%% Setup 
%check that w and window_size match 
if ~isequal(size(w, 2), window_size)
    error("The weights vector is not equal to the window size")
end

%allocate memory
rho = NaN(n_windows, 1);

%% Calculation
for tau = 1:n_windows %index for windows 
    %Calculate starting time
    t = (dt_step * tau) - (dt_step - 1);
    %get subset of data used for correlation coefficients 
    X = R(t:t+window_size-1, :);
    %correlation coefficients
    C = weightedcorrs(X, w);
    %remove anything not significant 
    pvals = bootstrap_correlation_signifiance(X, C, 100, @(x) weightedcorrs(x, w));
    C(pvals >= 0.05) = NaN;
    rho(tau) = mean(utri_to_vec(C), 'omitnan'); %mean of upper triangle
end
end