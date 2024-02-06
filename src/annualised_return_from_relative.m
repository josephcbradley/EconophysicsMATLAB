function output = annualised_return_from_relative(X)
%% Description 
%Annualised return from relative returns
%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
N_not_nan = sum(~isnan(X));
cumret = prod(1 + X, 'omitnan');
output = 252 * (cumret^(1/N_not_nan) - 1);
end