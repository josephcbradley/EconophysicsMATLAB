function Q = qth_abs_moment(X, tau, q) 

arguments
    X (:, 1) double {must_be_real}
    tau (1, 1) double {must_be_integer, tau_small_enough(tau, X)}
    q (1, 1) double {must_be_real}
end
    %Calculate Q = E[|X(t+ tau) - X(t)|^q]
    T = length(X);
    %Q = E[|X(t+ tau) - X(t)|^q]
    Q = 0.0;
    %with gap τ
    for t = 1:T-tau
         Q = Q + abs(X(t + tau) - X(t))^q;
    end
    Q = Q / (T - tau);
end


%local validation 
function flag = tau_small_enough(tau, X)
    flag = tau < length(X);
end