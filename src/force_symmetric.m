function As = force_symmetric(A)
%% Description 
% Forces A to be symmetric when MATLAB can't do it!
%% Inputs
% A - a matrix
%% Ouputs 
% As - a that has been forced to be symmetric
%% Setup 
arguments
    A (:, :) double {must_be_real}
end

%% Calculation
As = (A + A') / 2;
end