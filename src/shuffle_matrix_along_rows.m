function shuffled_X = shuffle_matrix_along_rows(X)
%% Description 
%Shuffles the rows of X. Useful for bootstrap calculations.
%% Inputs
%X - an array
%% Ouputs 
%shuffled_X = X with the rows independently shuffled
%% Setup 
%NA

%% Calculation
    shuffled_X = X;
    s = rand(size(X));
    [~, perm] = sort(s, 2);
    for r = 1:size(X, 1)
        shuffled_X(r, :) = shuffled_X(r, perm(r, :));
    end
end