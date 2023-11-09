function [Hs, SEs] = generalised_hurst_range(X, tau_range, q_range)
    % L = length(q_range);
    % N = length(tau_range)
    % Y = Vector{Float64}(undef, N)
    % S = similar(Y)

    [zetas, SEs] = zeta_estimator_range(X, tau_range, q_range);
    Hs = zetas ./ q_range; %divide by q 
    % If the process does not actually scale in the assumed way, the linear regression method will throw values 
    % for H outside (0, 1). To make this clear, we will set any such values of the buffer to NaN 

    %set Hs 
    Hs(isnan(zetas)) = NaN;
    %set SDs
    SEs(isnan(zetas)) = NaN;
end