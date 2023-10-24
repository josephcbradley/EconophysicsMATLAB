function tests = calculate_rolling_mean_corrTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check output dimensions
    assert(isequal(size(testCase.TestData.rho), [length(testCase.TestData.t_range) 1]))
end

function testFunctionTwo(testCase)
    %check that ouput is bounded
    assert(max(testCase.TestData.rho) <= 1)
    assert(min(testCase.TestData.rho) >= -1)
    assert(isempty(find(~isreal(testCase.TestData.rho), 1)));
end

function setup(testCase)  
    %create data
    testCase.TestData.T = 1000;
    testCase.TestData.N = 10;
    testCase.TestData.dt_step = 20;
    testCase.TestData.dt = 50;
    testCase.TestData.R = randn(testCase.TestData.T, testCase.TestData.N);
    testCase.TestData.theta = testCase.TestData.dt / 3;
    testCase.TestData.w = generate_expweights(testCase.TestData.dt, ...
        testCase.TestData.theta);
    testCase.TestData.t_range = testCase.TestData.dt:testCase.TestData.dt_step:testCase.TestData.T;
    testCase.TestData.hndl = @(x) weighted_pearson_corrs(x, testCase.TestData.w);
    testCase.TestData.rho = calculate_rolling_mean_corr(testCase.TestData.R, ...
        testCase.TestData.dt, testCase.TestData.t_range, CorrHandle = testCase.TestData.hndl);
end