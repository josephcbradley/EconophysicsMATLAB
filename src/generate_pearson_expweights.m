function w = generate_pearson_expweights(dt, theta)
%% Description 
% Generate exponential weights for statistics, as described in Pozzi et al 
% 2012 (https://doi.org/10.1140/epjb/e2012-20697-x)
%% Inputs
% dt: the length of the time period 
% theta: the charateristic time of the correlations. A good default is dt/3
%% Ouputs 
%w: the 1-by-dt vector of weights
%% Setup 
%NA
%% Calculation
% Calculate constant w0
w0 = (1 - exp( -1 / theta)) / (1 - exp( -dt / theta)); 
% Calculate exponential weights w, for Pearson
w = arrayfun(@(t) w0 * exp((t - dt) / theta), 1:dt);
% For safety, ensure the sum of weights is 1
w = w / sum(w);
end