function output = annualised_return_from_logs(X)
%% Description 
%Annualised return from log returns
%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
output = annualised_return_from_relative(exp(X) - 1);
end