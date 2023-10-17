function alpha = rank_tail_estimator(X)
%% Description 
%Calculates the rank tail estimator of the vector X. 
%% Inputs
% X - a vector of data
%% Ouputs 
% alpha - the scaling parameter of the tails 
%% Setup 
%check x is all positive 
if ~isempty(find(X < 0, 1))
    error('Negative values present in X!')
end
T = length(X);

%% Calculation
    s = sort(X);
    Ts = arrayfun(@(t) t/T, flip(1:T));
    p = polyfit(log(s), log(Ts), 1);
    alpha = -p(1);
end