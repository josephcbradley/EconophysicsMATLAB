function flag = anynan(X)
%% Description 
%Checks to see whether that are any NaN data in the array X
%% Inputs
% X: an array
%% Ouputs 
% flag: true if there are any NaNs in X, false otherwise.
%% Setup 
%None

%% Calculation
flag = ~isempty(find(isnan(X), 1));
end