function tests = calculate_rolling_mean_partial_corrTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    data = testCase.TestData.data;
    assert(min(data) >= -1);
    assert(max(data) <= 1);
    assert(isempty(find(~isreal(data), 1)));
end

function setup(testCase)  
    %create data
    T = 1000;
    N = 10;
    X = randn(T, N);
    Y = randn(T, 1);
    dt = 100;
    theta = dt / 3;
    gap_dt = 30;
    t_range = dt:gap_dt:T;
    w = generate_expweights(dt, theta);
    handle = @(r, y) weighted_partialcorrs(r, y, w);
    testCase.TestData.data = calculate_rolling_mean_partial_corr(X, Y, dt, t_range, CorrHandle=handle);
end