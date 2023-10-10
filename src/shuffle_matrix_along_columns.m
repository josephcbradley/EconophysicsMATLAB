function shuffled_X = shuffle_matrix_along_columns(X)
%% Description 
%Shuffles the columns of X. Useful for bootstrap calculations.
%% Inputs
%X - an array
%% Ouputs 
%shuffled_X = X with the columns independently shuffled
%% Setup 
%NA

%% Calculation
    shuffled_X = X;
    s = rand(size(X));
    [~, perm] = sort(s, 1);
    for r = 1:size(X, 2)
        shuffled_X(:, r) = X(perm(:, r), r);
    end
end