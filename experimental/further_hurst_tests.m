%% Hurst tests
%These are long tests for validating the Hurst exponent methods but take
%too long to be part of a normal test suite
addpath('src/')
%% Student's t 
L = 1000000;
n_trials = 10000;
n_range = 3:5;
tstudent_estimated_Bs = NaN(n_trials, length(n_range));
tau_range = 30:250;
q_range = 0.1:0.1:2;
for t = 1:n_trials
    for n_idx = 1:length(n_range)
        n = n_range(n_idx);
        students_walk = cumsum(trnd(n, L, 1));
        zetas = zeta_estimator_range(students_walk, tau_range, q_range);
        mdl = polyfitn(q_range, zetas, 2);
        tstudent_estimated_Bs(t, n_idx) = mdl.Coefficients(1);
    end
end
mus = mean(tstudent_estimated_Bs, 1);
sigmas = std(tstudent_estimated_Bs, [], 1);
disp(compose("Mean B for t-student's with n = 3 is %.2f", mus(n_range == 3)))
disp(compose("Mean B for t-student's with n = 4 is %.2f", mus(n_range == 4)))
disp(compose("Mean B for t-student's with n = 5 is %.2f", mus(n_range == 5)))
% 
% figure
% errorbar(n_range, mus, sigmas)
% xlabel("$q$", 'interpreter', 'latex')       
% ylabel("$B$", 'interpreter', 'latex')
% title("Multiscaling parameter for $t$-Student's distribution with varying DOF", 'interpreter', 'latex')

%% MFRW 
lambda2 = 0.03;
autocor_length = 1000;
sigma = 1;
epsilon_distribution =  makedist('Normal','mu', 0 ,'sigma', sigma);
omega_distribution = makedist('Normal','mu', -lambda2 * log(autocor_length), 'sigma', lambda2 * log(autocor_length));
MFRW_estimated_Bs = NaN(n_trials, 1);
MFRW_estimated_Hs = NaN(n_trials, 1);

for t = 1:n_trials
    epsilons = random(epsilon_distribution, L, 1);
    omegas = random(omega_distribution, L, 1);
    MFRW_increments = epsilons .* omegas;
    MFRW = cumsum(MFRW_increments);
    [mfrw_zetas, mfrw_z_SEs] = zeta_estimator_range(MFRW, tau_range, q_range);
    mfrw_mdl = polyfitn(q_range, mfrw_zetas, 2);
    MFRW_estimated_Bs(t) = mfrw_mdl.Coefficients(1);
    MFRW_estimated_Hs(t) = mfrw_zetas(q_range == 1);
end

mus = mean(MFRW_estimated_Bs);
sigmas = std(MFRW_estimated_Bs);

disp(compose("Mean B for MFRW is %.2f", mus))