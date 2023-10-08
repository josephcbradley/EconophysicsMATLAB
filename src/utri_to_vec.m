function vec = utri_to_vec(A)
%% Description 
% Extract the upper triangle of a matrix (not including diagonal)

%% Inputs
% A - a matrix

%% Ouputs 
% vec - the length N*(N-1)/2 vectorcorresponding to the upper triangle of
% the matrix
%% Setup 
%NA

%% Calculation
m = triu(true(size(A)), 1);
vec = A(m);
end

