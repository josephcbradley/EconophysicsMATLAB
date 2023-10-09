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
    
    for r = 1:size(X, 2)
        s = rand(size(X, 1), 1);
        [~, perm] = sort(s);
        shuffled_X(:, r) = X(perm, r);
    end
end