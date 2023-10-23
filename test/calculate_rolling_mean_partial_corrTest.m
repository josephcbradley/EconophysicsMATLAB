function tests = calculate_rolling_mean_partial_corrTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    data = testCase.TestData.data;
    assert(min(data) >= -1);
    assert(max(data) <= 1);
    assert(isempty(find(~isreal(data), 1)));
end

function testFunctionTwo(testCase)
    assert(mean(testCase.TestData.AB_rho) < 1e-3)
end

function setup(testCase)  
    %create data
    T = 100000;
    N = 10;
    X = randn(T, N);
    Y = randn(T, 1);
    dt = 100;
    theta = dt / 3;
    gap_dt = 30;
    t_range = dt:gap_dt:T;
    w = generate_expweights(dt, theta);
    handle = @(r, y) weighted_partialcorrs(r, y, w);
    testCase.TestData.data = calculate_rolling_mean_partial_corr(X, Y, dt, t_range, ...
        CorrHandle=handle, RemoveInsignificant = false);
    %data that are partial independent
    testCase.TestData.A = 0.1 + (0.3 * Y) + randn(T, 1);
    testCase.TestData.B = 0.5 + (1.2 * Y) + randn(T, 1);
    testCase.TestData.AB_rho = calculate_rolling_mean_partial_corr([testCase.TestData.A testCase.TestData.B], Y, dt, t_range, ...
        CorrHandle = handle, RemoveInsignificant = false);
end