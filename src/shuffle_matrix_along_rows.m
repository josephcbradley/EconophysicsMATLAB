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
    
    for r = 1:size(X, 1)
        s = rand(1, size(X, 2));
        [~, perm] = sort(s);
        shuffled_X(r, :) = shuffled_X(r, perm);
    end
end