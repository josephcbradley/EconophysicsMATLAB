function output = annualised_risk_from_logs(X)
%% Description 
%Annualised risk from log returns
%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
output = annualised_risk_from_relative(exp(X) - 1);
end