function range_buffer = zeta_estimator_range(X, tau_range, q_range)

arguments
    X(:, 1) double {must_be_real}
    tau_range (:, 1) {must_be_integer}
    q_range (:, 1) {must_be_real}
end

    L = length(q_range);
    range_buffer = NaN(L, 2);
    for i = 1:L
        q = q_range(i);
        [range_buffer(i, 1), range_buffer(i, 2)] = zeta_estimator(X, tau_range, q);
    end
end