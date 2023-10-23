function tests = single_factor_residualsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check sizes
    assert(isequal(size(testCase.TestData.X), [testCase.TestData.T, testCase.TestData.N]))
    
end

function testFunctionTwo(testCase)
    %check means are zero
    assert(isequal(abs(mean(testCase.TestData.X_residuals, 1)) < 1e-3, [true true]))
end

function testFunctionThree(testCase)
    %check means are zero
    C = corr(testCase.TestData.X_residuals);
    rho_estimated = C(1, 2);
    rho_actual = testCase.TestData.sigma(1, 2);
    assert(abs(rho_estimated - rho_actual) < 1e-2)
end

function testFunctionFour(testCase)
    raw_C = corr([testCase.TestData.A testCase.TestData.B]);
    raw_rho = raw_C(1, 2);
    residual_C = corr(testCase.TestData.AB_residuals);
    residual_rho = residual_C(1, 2);
    assert(raw_rho > 0)
    assert(residual_rho < 1e-2)
end

function setup(testCase)  
    %create data
    testCase.TestData.T = 100000; 
    testCase.TestData.N = 2;
    testCase.TestData.mu = [0 0];
    testCase.TestData.sigma = [1 0.5; 0.5 1];
    testCase.TestData.X = mvnrnd(testCase.TestData.mu, testCase.TestData.sigma, testCase.TestData.T);
    testCase.TestData.Y = randn(testCase.TestData.T, 1);
    testCase.TestData.X_residuals = single_factor_residuals(testCase.TestData.X, testCase.TestData.Y);
    testCase.TestData.A = 0.1 + (0.3 * testCase.TestData.Y) + randn(testCase.TestData.T, 1);
    testCase.TestData.B = 0.5 + (1.2 * testCase.TestData.Y) + randn(testCase.TestData.T, 1);
    testCase.TestData.AB_residuals = single_factor_residuals([testCase.TestData.A testCase.TestData.B], testCase.TestData.Y);
end