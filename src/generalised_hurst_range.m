function [Hs, SEs] = generalised_hurst_range(X, tau_range, q_range)
    % L = length(q_range);
    % N = length(tau_range)
    % Y = Vector{Float64}(undef, N)
    % S = similar(Y)

    buffer = zeta_estimator_range(X, tau_range, q_range);
    buffer(:, 1) = buffer(:, 1) ./ q_range; %divide by q 
    % If the process does not actually scale in the assumed way, the linear regression method will throw values 
    % for H outside (0, 1). To make this clear, we will set any such values of the buffer to NaN 

    %set Hs 
    buffer(isnan(buffer(:, 1)), 1) = NaN;
    %set SDs
    buffer(isnan(buffer(:, 1)), 2) = NaN;
    Hs = buffer(:, 1);
    SEs = buffer(:, 2);
end