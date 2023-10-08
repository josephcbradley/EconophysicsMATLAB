function tests = calculate_rolling_mean_corrTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check output dimensions
    assert(isequal(size(testCase.TestData.rho), [testCase.TestData.n_windows 1]))
end

function testFunctionTwo(testCase)
    %check that ouput is bounded
    assert(max(testCase.TestData.rho) <= 1)
    assert(min(testCase.TestData.rho) >= -1)
end

function setup(testCase)  
    %create data
    testCase.TestData.T = 1000;
    testCase.TestData.N = 10;
    testCase.TestData.n_windows = 10;
    testCase.TestData.dt_step = 20;
    testCase.TestData.window_size = 50;
    testCase.TestData.R = randn(testCase.TestData.T, testCase.TestData.N);
    testCase.TestData.theta = testCase.TestData.window_size / 3;
    testCase.TestData.w = generate_expweights(testCase.TestData.window_size, ...
        testCase.TestData.theta);
    testCase.TestData.rho = calculate_rolling_mean_corr(testCase.TestData.R, ...
        testCase.TestData.w, testCase.TestData.n_windows, testCase.TestData.dt_step, ...
        testCase.TestData.window_size);
end