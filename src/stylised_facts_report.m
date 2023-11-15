function stylised_facts_report(R, P, dates)
%% Check inputs
arguments
    R (:, :) double {must_be_real}
    P(:, :) double {must_be_real}
    dates (:, 1) datetime {lengths_must_match(R, dates)}
end

%% Check for Econometrics toolbox 
econometrics_available = license('test', 'econometrics_toolbox');
if ~econometrics_available
    warning("The Econometrics toolbox is not available - some tests will be skipped.")
end
%% Set up dir for plots
current_date_formatted = datetime('now','TimeZone','local','Format','yyyyMMddHHmmss');
report_name = input("Please enter the title of your report or dataset:\n", "s");
report_name = clean_report_name(report_name);
folder_name = report_name + "_stylised_facts_" + string(current_date_formatted);

if isfolder(folder_name)
    error("A folder for this report already seems to exist. Stopping.")
else 
    mkdir(folder_name)
    old_wd = pwd;
    cd(folder_name)
end

% setup up plot counter 

current_plot_number = 1;

%% Stationarity & Market Mode
[T, N] = size(R);
disp(compose("R has %i time points and %i stocks.", T, N))
if econometrics_available
    %check for stationarity 
    adf_results = zeros(N, 1);
    for i = 1:N 
        if sum(~isnan(R(:, i))) > 10 %check that we have enough data points
            %for ADF test
            adf_results(i) = adftest(R(:, i));
        end
        if adf_results(i) == 0
            warning(compose("ADF test for series %i is non-stationary.", i))
        end
    end
else
    warning("Skipping ADF tests as the Econometrics toolbox is not available.")
end

market_mode = mean(R, 2);
market_mode_plot = figure('visible','off');
plot(dates, market_mode)
xlabel("Date")
ylabel("Market mode")
title("Market mode")
exportgraphics(market_mode_plot, compose("%05i", current_plot_number) + "_market_mode.png", 'resolution', 300)
current_plot_number = current_plot_number + 1;

%% Heavy tails
kurtoses = NaN(N, 1);
for i = 1:N
    kurtoses(i) = kurtosis(R(:, i));
end

%plot 
kurtoses_plot = figure('visible','off');
histogram(kurtoses, 'Normalization', 'probability')
xlabel("$\kappa$", 'Interpreter', 'latex')
ylabel("$P(\kappa)$", 'Interpreter', 'latex')
title("Distribution of kurtoses")
exportgraphics(kurtoses_plot, compose("%05i", current_plot_number) + "_kurtoses_histogram.png", 'resolution', 300)
current_plot_number = current_plot_number + 1;

positive_tails = NaN(N, 1);
negative_tails = NaN(N, 1);

parfor i = 1:N
    positive_indices = R(:, i) > 0;
    positive_tails(i) = plfit(R(positive_indices, i), 'sample', 100);
    negative_indicies = R(:, i) < 0;
    negative_tails(i) = plfit(abs(R(negative_indicies, i)), 'sample', 100);
end

%check what largest and smallest exponent is 
max_exponent = max(max(positive_tails), max(negative_tails));
if max_exponent > 7
    warning("There are tails greater than 7 in the data - this is unusual for stock returns and should be checked.")
end

min_exponent = min(min(positive_tails), min(negative_tails));
if min_exponent < 1
    warning("There are tails less than 1 in the data - this is unusual for stock returns and should be checked.")
end

%plot distribution 
positive_tails_plot = figure('visible','off');
histogram(positive_tails, 'BinWidth', 0.1, 'Normalization', 'probability')
xlim([1 7])
xlabel("$\alpha$", 'Interpreter', 'latex')
ylabel("$P(\alpha)$", 'Interpreter', 'latex')
title("Distribution of positive tail exponents")
exportgraphics(positive_tails_plot, compose("%05i", current_plot_number) + "_positive_tails_histogram.png", 'resolution', 300)
current_plot_number = current_plot_number + 1;

negative_tails_plot = figure('visible','off');
histogram(negative_tails, 'BinWidth', 0.1, 'Normalization', 'probability')
xlim([1 7])
xlabel("$\alpha$", 'Interpreter', 'latex')
ylabel("$P(\alpha)$", 'Interpreter', 'latex')
title("Distribution of negative tail exponents")
exportgraphics(negative_tails_plot, compose("%05i", current_plot_number) + "_negative_tails_histogram.png", 'resolution', 300)
current_plot_number = current_plot_number + 1;

% aggregational Gaussianity
dt_range = 1:10;
mean_positive_tail = NaN(length(dt_range), 1);
std_positive_tail = NaN(length(dt_range), 1);
for dt_idx = 1:length(dt_range)
    dt = dt_range(dt_idx);
    Rdt = movsum(R, [dt 1], 1, "Endpoints", "discard");
    tails = NaN(N, 1);
    parfor i = 1:N
        pos_indices = Rdt(:, i) > 0;
        tails(i) = plfit(Rdt(pos_indices, i), 'sample', 100);
    end
    mean_positive_tail(dt_idx) = mean(tails);
    std_positive_tail(dt_idx) = std(tails);
end

agg_gaussianity_plot = figure('visible','off');
plot(dt_range, mean_positive_tail)
xlabel("$\Delta t$", 'Interpreter', 'latex')
ylabel("$\alpha^\star(\Delta t)$", 'Interpreter', 'latex')
title("Aggregational Gaussianity: mean $\alpha(\Delta t)$ of positive tail", 'Interpreter', 'latex')
exportgraphics(agg_gaussianity_plot, compose("%05i", current_plot_number) + "_aggregational_gaussianity.png", 'resolution', 300)
current_plot_number = current_plot_number + 1;

%Conditional heavy tails
if econometrics_available
    garch_1_1 = garch(1, 1);
    %garch_params = NaN(N, 3);
    conditional_kurtoses = NaN(N, 1);
    for i = 1:N
        est_mdl = estimate(garch_1_1, R(:, i), 'Display', 'off');
        % garch_params(i, 1) = est_mdl.Constant;
        % garch_params(i, 2) = est_mdl.GARCH{1};
        % garch_params(i, 3) = est_mdl.ARCH{1};
        v = infer(est_mdl, R(:, i));
        res = (R(:, i) - est_mdl.Offset)./sqrt(v);
        conditional_kurtoses(i) = kurtosis(res);
    end
    conditional_kurtoses_plot = figure('visible','off');
    histogram(conditional_kurtoses, 'Normalization', 'probability')
    xlabel("$\kappa$", 'Interpreter', 'latex')
    ylabel("$P(\kappa)$", 'Interpreter', 'latex')
    title("Distribution of conditional kurtoses")
    exportgraphics(conditional_kurtoses_plot, compose("%05i", current_plot_number) + "_conditional_kurtoses_histogram.png", 'resolution', 300)
    current_plot_number = current_plot_number + 1;
else
    warning("Skipping GARCH tests as the Econometrics toolbox is not available.")
end

%% Abscence of autocorrelations 
tau_max = 40;
Corr = zeros(N, tau_max + 1);
Corr_abs = zeros(N, tau_max + 1);
Corr_sq = zeros(N, tau_max + 1);

for i=1:N
    Corr(i, :) = autocorr(R(:, i), 'NumLags', tau_max);
    Corr_abs(i, :) = autocorr(abs(R(:, i)), 'NumLags', tau_max);
    Corr_sq(i, :) = autocorr(abs(R(:, i)) .^ 2, 'NumLags', tau_max);
end

%N.B. rows index stocks, columns index lags
av_corr = mean(Corr, 1);
sigma_corr = std(Corr, [], 1);
av_corr_abs = mean(Corr_abs, 1);
sigma_corr_abs = std(Corr_abs, [], 1);
av_corr_sq = mean(Corr_sq, 1);
sigma_corr_sq = std(Corr_sq, [], 1);

autocorrs_figure = figure('Visible','off');
%plot(0:tau_max, av_corr,"-",LineWidth=1.0,MarkerSize=12)
%hold on
%errorbar(0:tau_max, av_corr,sigma_corr,".-",LineWidth=1.0,MarkerSize=9)
hold on 
errorbar(0:tau_max, av_corr, sigma_corr, ".-");
%legend(["\Deltat = "+num2str(dt)+" log-returns"])
%title("Average autocorrelation function of daily log-returns")
% hold on
% errorbar(0:tau_max, av_corr_abs,sigma_corr_abs,".-",LineWidth=1.0,MarkerSize=9);
%errorbar(0:tau_max, av_corr_sq,sigma_corr_sq,".-",LineWidth=1.0,MarkerSize=9);
errorbar(0:tau_max, av_corr_sq, sigma_corr_sq, ".-");
errorbar(0:tau_max, av_corr_abs, sigma_corr_abs, ".-");
hold off
xlabel("$\tau$", 'interpreter', 'latex')
ylabel("$\rho^\star(\tau)$", 'interpreter', 'latex')
title("Average autocorrelations of different functions of the returns")
legend(["$\rho(r(i, t),r(i, t+ \tau))$","$\rho(|r(i, t)|^2,|r(t+\tau)|^2)$","$\rho(|r(i, t)|,|r(t+\tau)|)$"], 'interpreter', 'latex')
exportgraphics(autocorrs_figure, compose("%05i", current_plot_number) + "_decay_of_autocorrs.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

% %% Multiscaling
% q_range = (0.1:0.1:2)';
% tau_range = (30:250)';
% mean_price = mean(P, 2);
% [zeta, SEs] = zeta_estimator_range(mean_price, tau_range, q_range);
% %fit 2nd order polynomial to qHs (which are a function of qs
% mdl = polyfitn(q_range, zeta, 2);
% B = mdl.Coefficients(1);
% B_SE = mdl.ParameterStd(1);
% A = mdl.Coefficients(2);
% A_SE = mdl.ParameterStd(2);
% 
% market_mode_scaling_figure = figure('Visible', 'off');
% errorbar(q_range, zeta, SEs)
% xlabel("$q$", 'Interpreter', 'latex')
% ylabel("$qH(q)$", 'Interpreter', 'latex')
% title("Generalised Hurst exponent for the market mode" + newline + ...
%     compose("$B = %.5f \\pm %.5f, A = %.5f \\pm %.5f$", B, B_SE, A, A_SE), ...
%     'Interpreter','latex')
% exportgraphics(market_mode_scaling_figure, compose("%05i", current_plot_number) + "_mktmode_mutliscaling.png", "Resolution", 300)
% current_plot_number = current_plot_number + 1;
% 
% %calculate Hurst exponent for all stocks 
% zetas = NaN(N, length(q_range)); 
% SE_zetas = NaN(N, length(q_range));
% Bs = NaN(N, 1);
% SE_Bs = NaN(N, 1);
% for i = 1:N 
%     [zetas(i, :), SE_zetas(i, :)] = generalised_hurst_range(P(:, i), tau_range, q_range);
%     mdl = polyfitn(q_range, zetas(i, :), 2);
%     Bs(i) = mdl.Coefficients(1);
%     SE_Bs(i) = mdl.ParameterStd(1);
% end
% 
% H1_distribution_figure = figure('Visible', 'off'); 
% histogram(zetas(:, q_range == 1), 'Normalization','probability')
% xlabel("$H(1)$", 'Interpreter', 'latex')
% ylabel("$P(H(1))$", 'Interpreter', 'latex')
% title("Histribution of Hurst exponents across stocks")
% exportgraphics(H1_distribution_figure, compose("%05i", current_plot_number) + "_H1_distribution.png", "Resolution", 300)
% current_plot_number = current_plot_number + 1;
% 
% B_distribution_figure = figure("Visible", "off");
% histogram(Bs, 'Normalization', 'probability')
% xlabel("$B$", 'Interpreter','latex')
% ylabel("$P(B)$", 'Interpreter','latex')
% title("Distribution of multiscaling parameter $B$ across stocks", 'Interpreter','latex')
% exportgraphics(B_distribution_figure, compose("%05i", current_plot_number) + "_B_distribution.png", "Resolution", 300)
% current_plot_number = current_plot_number + 1;

%% Momentum 
% J = 252;
% K = 63;
% L_max = 200;
% L_results = NaN(L_max, 1);
% for L = 1:L_max
% 
%     %ranks for J day returns (backward)
%     R_J = [NaN(J - 1, N); movmean(R, [J - 1, 0], 1, 'omitmissing', 'Endpoints','discard')];
%     ranks_J = NaN(T, N);
%     for t = J:T
%         [~, p] = sort(R_J(t, :), 'ascend');
%         for i = 1:N
%             ranks_J(t, i) = find(p == i);
%         end
%     end
% 
%     %ranks for K day returns (forward)
%     R_K = [movmean(R, [0, K - 1], 1, 'omitmissing', 'Endpoints','discard'); NaN(K - 1, N)];
%     ranks_K = NaN(T, N);
%     for t = 1:T-K+1
%         [~, p] = sort(R_K(t, :), 'ascend');
%         for i = 1:N
%             ranks_K(t, i) = find(p == i);
%         end
%     end
% 
% 
%     C = corr(R_J(J:T-L-K+1, :), R_K(J+L:T-K+1, :));
%     L_results(L) = mean(diag(C));
% end
% momentum_figure = figure('Visible','off'); 
% plot(1:L_max, L_results)
% xlabel("$L$", 'Interpreter','latex')
% ylabel("Mean correlation coefficient between $J$ and $K$", 'Interpreter','latex')
% title("Momentum rank correlations")
% exportgraphics(momentum_figure, compose("%05i", current_plot_number) + "_momentum_rank_corrs.png", "Resolution", 300)
% current_plot_number = current_plot_number + 1;

%% Interdependence 
%eigenvalue distribution 
C = corr(R);
[V, d] = eig(C, 'vector');
eigenvalue_distribution_histogram = figure('Visible', 'off');
histogram(d, 'BinWidth', 0.1, 'Normalization', 'probability')
xlabel("$\lambda$", 'Interpreter', 'latex')
ylabel("$P(\lambda)$", 'Interpreter', 'latex')
title("Eigenvalue distribution of correlation matrix")
xlim([0 5])
axes('Position',[.4 .4 .5 .5])
box on
histogram(d, 'Normalization','pdf', 'BinWidth', 0.1)
xls = xlim;
xlim([5 xls(2)])
ylim([0 0.05])
xlabel("$\lambda$", 'Interpreter', 'latex')
ylabel("$P(\lambda)$", 'Interpreter', 'latex')
exportgraphics(eigenvalue_distribution_histogram, compose("%05i", current_plot_number) + "_lambda_distribution.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

%effect of market mode - see correlation coefficients
P = bootstrap_correlation_signifiance(R, C, 100, @(x) corr(x));
C_vec = utri_to_vec(C);
P_vec = utri_to_vec(P);

resids = single_factor_residuals(R, market_mode);
resids_C = corr(resids);
resids_P = bootstrap_correlation_signifiance(resids, resids_C, 100, @(x) corr(x));
resids_C_vec = utri_to_vec(resids_C);
resids_P_vec = utri_to_vec(resids_P);

bin_width = 0.01;
detrending_corrs_figure = figure('Visible', 'off');
tiledlayout("flow")
nexttile
hold on
histogram(C_vec(P_vec < 0.05), 'FaceColor', 'blue', 'BinWidth', bin_width, ...
    'Normalization', 'count')
histogram(C_vec(P_vec >= 0.05), 'FaceColor', 'yellow', 'BinWidth', bin_width, ...
    'Normalization', 'count')
hold off
xlim([-0.5 1])
xlabel("$C_{i, j}$", 'Interpreter', 'latex')
ylabel("Number of $C_{i, j}$", 'Interpreter', 'latex')
legend("Significant", "Insignificant")
title("Distribution of correlation coefficients - log returns")
nexttile
hold on
histogram(resids_C_vec(resids_P_vec < 0.05), 'FaceColor', 'blue', 'BinWidth', bin_width, ...
    'Normalization', 'count')
histogram(resids_C_vec(resids_P_vec >= 0.05), 'FaceColor', 'yellow', 'BinWidth', bin_width, ...
    'Normalization', 'count')
hold off
xlim([-0.5 1])
xlabel("$C_{i, j}$", 'Interpreter', 'latex')
ylabel("Number of $C_{i, j}$", 'Interpreter', 'latex')
legend("Significant", "Insignificant")
title("Distribution of correlation coefficients - residuals from market mode")
exportgraphics(detrending_corrs_figure, compose("%05i", current_plot_number) + "_detrended_corrs_distribution.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

%epps affect 
%1. mean value of corr
dt_range = 2:252;
dt_mean_C = NaN(length(dt_range), 1);
dt_std_C = NaN(length(dt_range), 1);
dt_largest_eigenvalue = NaN(length(dt_range), 1);
for dt_idx = 1:length(dt_range)
    dt = dt_range(dt_idx);
    Rdt = movsum(R, [0 dt - 1], 1, 'omitnan', 'Endpoints', 'discard');
    C = corr(Rdt);
    dt_mean_C(dt_idx) = mean(utri_to_vec(C));
    dt_std_C(dt_idx) = std(utri_to_vec(C));
    [~, d] = eig(C, 'vector');
    dt_largest_eigenvalue(dt_idx) = max(d);
end

mean_corr_figure = figure('Visible','off');
errorbar(dt_range, dt_mean_C, dt_std_C)
xlabel("$\Delta t$", 'Interpreter','latex')
ylabel("Mean of $C_{i, j}$", 'Interpreter','latex')
title("Size of mean correlation coefficient as $\Delta t$ grows", 'Interpreter','latex')
exportgraphics(mean_corr_figure, compose("%05i", current_plot_number) + "_mean_corr.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

largest_eigenvalue_figure = figure('Visible','off');
plot(dt_range, dt_largest_eigenvalue/N)
xlabel("$\Delta t$", 'Interpreter','latex')
ylabel("$\lambda_{max}$", 'Interpreter','latex')
title("Size of largest eigenvalue as $\Delta t$ increases", 'Interpreter','latex')
exportgraphics(largest_eigenvalue_figure, compose("%05i", current_plot_number) + "_largest_eigenvalue.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

%% Persistent topology
%calculate an MST and a PMFG for every T 
dt = 1000;
w = generate_pearson_expweights(dt, dt / 3);
t_range = dt:30:T;
adj_matrix_size = N * (N - 1) * 0.5;
mst_edges = zeros(adj_matrix_size, length(t_range), 'uint64');
pmfg_edges = zeros(adj_matrix_size, length(t_range), 'uint64');
corr_coeffs_over_time = zeros(adj_matrix_size, length(t_range));
for t_idx = 1:length(t_range)
    t = t_range(t_idx);
    local_tspan = t-dt+1:t;
    X = R(local_tspan, :);
    C = weighted_pearson_corrs(X, w); 
    corr_coeffs_over_time(:, t_idx) = utri_to_vec(C);
    %C = force_symmetric(C); %shouldn't be needed for our weighted corrs
    D = sqrt_dissimilarity(C);
    M = mst(sparse(D));
    mst_edges(:, t_idx) = full(utri_to_vec(M));
    P = PMFG_T2s(C + 1);
    pmfg_edges(:, t_idx) = full(utri_to_vec(P));
end

meta_C = corr(corr_coeffs_over_time);
cls = [0 1];
metacorrs_figure = figure('Visible','off');
heatmap(meta_C)
colormap("turbo")
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 30) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("Heatmap of metacorrelations")
exportgraphics(metacorrs_figure, compose("%05i", current_plot_number) + "_metacorrs_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

mst_corr = zeros(length(t_range)); %square!
pmfg_corr = zeros(length(t_range)); %square!
for i = 2:length(t_range)
    for j = 1:i-1
        mst_corr(i, j) = (1/(N-1)) * raw_ES_from_upper(mst_edges(:, i), mst_edges(:, j));
        pmfg_corr(i, j) = (1 / (3*N - 6)) * raw_ES_from_upper(pmfg_edges(:, i), pmfg_edges(:, j));
    end
end

mst_corr = mst_corr + mst_corr';
pmfg_corr = pmfg_corr + pmfg_corr';
for i = 1:length(t_range)
    mst_corr(i, i) = 1;
    pmfg_corr(i, i) = 1;
end

pmfg_ES_figure = figure('Visible','off');
heatmap(pmfg_corr)
colormap(turbo)
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 30) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("ES for PMFG edges")
exportgraphics(pmfg_ES_figure, compose("%05i", current_plot_number) + "_pmfg_ES_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

mst_ES_figure = figure('Visible','off');
heatmap(mst_corr)
colormap(turbo)
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 30) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("ES for MST edges")
exportgraphics(mst_ES_figure, compose("%05i", current_plot_number) + "_mst_ES_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

mst_residuals_edges = zeros(adj_matrix_size, length(t_range), 'uint64');
pmfg_residuals_edges = zeros(adj_matrix_size, length(t_range), 'uint64');
residuals_corr_coeffs_over_time = zeros(adj_matrix_size, length(t_range));
for t_idx = 1:length(t_range)
    t = t_range(t_idx);
    local_tspan = t-dt+1:t;
    X = resids(local_tspan, :);
    C = weighted_pearson_corrs(X, w); 
    residuals_corr_coeffs_over_time(:, t_idx) = utri_to_vec(C);
    %C = force_symmetric(C); %shouldn't be needed for our weighted corrs
    D = sqrt_dissimilarity(C);
    M = mst(sparse(D));
    mst_residuals_edges(:, t_idx) = full(utri_to_vec(M));
    P = PMFG_T2s(C + 1);
    pmfg_residuals_edges(:, t_idx) = full(utri_to_vec(P));
end

resids_meta_C = corr(residuals_corr_coeffs_over_time);
resids_metacorrs_figure = figure('Visible','off');
heatmap(resids_meta_C)
colormap("turbo")
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 15) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("Heatmap of metacorrelations for residuals")
exportgraphics(resids_metacorrs_figure, compose("%05i", current_plot_number) + "_resids_metacorrs_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

mst_residuals_corr = zeros(length(t_range)); %square!
pmfg_residuals_corr = zeros(length(t_range)); %square!
for i = 2:length(t_range)
    for j = 1:i-1
        mst_residuals_corr(i, j) = (1/(N-1)) * raw_ES_from_upper(mst_residuals_edges(:, i), mst_residuals_edges(:, j));
        pmfg_residuals_corr(i, j) = (1 / (3*N - 6)) * raw_ES_from_upper(pmfg_residuals_edges(:, i), pmfg_residuals_edges(:, j));
    end
end

mst_residuals_corr = mst_residuals_corr + mst_residuals_corr';
pmfg_residuals_corr = pmfg_residuals_corr + pmfg_residuals_corr';
for i = 1:length(t_range)
    mst_residuals_corr(i, i) = 1;
    pmfg_residuals_corr(i, i) = 1;
end

pmfg_residuals_ES_figure = figure('Visible','off');
heatmap(pmfg_residuals_corr)
colormap(turbo)
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 15) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("ES for PMFG edges from residuals")
exportgraphics(pmfg_residuals_ES_figure, compose("%05i", current_plot_number) + "_pmfg_residuals_ES_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;

mst_residuals_ES_figure = figure('Visible','off');
heatmap(mst_residuals_corr)
colormap(turbo)
caxis(cls)
ax = gca;
ax.XData = dates(t_range);
ax.YData = dates(t_range);
for i = 0:length(t_range)-1
    if mod(i, 30) == 0
        continue
    else
        ax.XDisplayLabels{i+1} = {''};
        ax.YDisplayLabels{i+1} = {''};
    end
end
title("ES for MST edges from residuals")
exportgraphics(mst_residuals_ES_figure, compose("%05i", current_plot_number) + "_mst_residuals_ES_heatmap.png", "Resolution", 300)
current_plot_number = current_plot_number + 1;


%% Cleanup 
cd(old_wd)
end

function lengths_must_match(R, x)
    if ~isequal(size(R, 1), length(x))
        eid = 'Size:notEqual';
        msg = 'Size of returns must be equal to lengths of other vectors.';
        error(eid,msg)
    end
end