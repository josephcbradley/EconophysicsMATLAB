function [zeta, SD] = lsq_estimator(Y, S)
arguments
    Y (:, 1) double {must_be_real}
    S (:, 1) double {must_be_real, must_have_same_length(Y, S)}
end

    N = length(Y);
    Ss = sum(S);
    Sy = sum(Y);
    Ssy = dot(S, Y);
    Sss = dot(S, S);
    Syy = dot(Y, Y);
    zeta = (N * Ssy - Ss * Sy) / (N * Sss - Ss^2);
    se2 = (1 / (N * (N - 2))) * (N * Syy - Sy^2 - zeta^2 * (N * Sss - Ss^2));
    SD = sqrt((N * se2) / (N * Sss - Ss^2));
end

function flag = must_have_same_length(Y, S)
    flag = isequal(length(Y), length(S));
end