function output = annualised_risk_from_relative(X)
%% Description 
%Annualised risk from relative returns
%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
output = sqrt(252) * std(X, 'omitnan');
end