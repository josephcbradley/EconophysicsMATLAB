function X_F = nearest_pos_semi_definite(B)
%% Description 
%Calculates the nearest positive semi-definite matrix to B. 
%Nearess us judged by the Frobenius norm, using the method and code given
%in:
%https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
%% Inputs
% B - a real numeric matrix
arguments
    B (:, :) double {must_be_real}
end

%% Ouputs 
%X_F - nearest positive semi-definite matrix to B

%% Setup 
% None

%% Calculation
[Q, d] = eig(B, 'vector'); %spectral decomposition
X_F = Q * (max(d, eps) .* Q'); %ensure all eigenvalues are at least eps
X_F = force_symmetric(X_F);

end