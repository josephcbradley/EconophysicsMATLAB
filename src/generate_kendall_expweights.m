function [w, i2, i1] = generate_kendall_expweights(dt, theta)
%% Description 

%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
% Calculate constant w0 
w0 = (exp(1/theta) + 1) * ((exp(1/theta) - 1)^2) / exp(2/theta) / (1 - exp(-dt/theta)) / (1 - exp((1 - dt)/theta)); 
% Indexes for all dt * (dt - 1) / 2 combinations without repetition: % {1, 2}, {1, 3}, ..., {1, dt}, {2, 3}, {2, 4}, ..., {2, dt}, ..., {dt - 1, dt} 
[i2, i1] = find(tril(ones(dt, 'uint8'), -1)); % Calculate exponential weights w, for Kendall 
w(:, 1) = w0 * exp((i1 + i2 - 2 * dt) / theta); % Ensure sum of weights is 1 
w=w/sum(w);
end