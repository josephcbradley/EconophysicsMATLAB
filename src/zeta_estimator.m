function [zeta, SD] = zeta_estimator(X, tau_range, q)
    %ζ(q) = qH(q) where H is the GHE 
    %let Q = E[|X(t+ tau) - X(t)|^q] = K(q)tau^ζ(q)
    %take logarithms and rearrange
    %ln(Q) = ζ(q)ln(tau) + ln(K(q))
    %to estimate ζ(q), consider the linear regression of
    %Y = β₀ + β₁S where 
    %Y = ln(Q), S = ln(tau), β₀ = ln(K(q)), and β₁ = ζ(q)
    %return β₁ and its std error

    %form Y 
    N = length(tau_range);
    Y = NaN(N, 1);
    S = NaN(N, 1);

    %calculate regression data
    for i = 1:N
        tau = tau_range(i);
        Y(i) = log(qth_abs_moment(X, tau, q));
        S(i) = log(tau);
    end

    %simple regression formulae
    [zeta, SD] = lsq_estimator(Y, S);

end