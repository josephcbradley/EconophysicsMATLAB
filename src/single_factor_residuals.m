function [residuals, alphas, betas] = single_factor_residuals(X, Y)
%% Description 
%Calculate the residuals of regressing each colum of X on Y
%Optionally, also returns the alphas and betas of this regression, i.e.
%X(:, i) = alpha(i) + beta(i)*Y(:) + epsilon(:, i)

%% Inputs
% X - a dt-by-N matrix of inputs 
% Y - a dt-by-1 factor

%% Ouputs 
% residuals - the dt-by-N matrix of residuals
% alphas - the N-by-1 vector of constant terms 
% betas - the N-by-1 vector of slope term terms 

%% Setup 
[dt, N] = size(X);

%compatible sizes
if ~isequal(size(Y), [dt, 1])
    error("Y is not equal to dt-by-1")
end

%check for NaNs
if anynan(X)
    warning("NaNs in X!")
elseif anynan(Y)
    warning("NaNs in Y!")
end


%% Calculation
%allocate
residuals = NaN(size(X));
alphas = NaN(N, 1);
betas = NaN(N, 1);

for i = 1:N
    available_data = ~isnan(X(:, i)) & ~isnan(Y);
    p = polyfit(Y(available_data), X(available_data, i), 1);
    residuals(:, i) = X(:, i) - polyval(p, Y);
    alphas(i) = p(1);
    betas(i) = p(2);
end
end