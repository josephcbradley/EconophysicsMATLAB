function tests = generalised_hurst_exponentTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    X = testCase.TestData.gaussian_walk;
    tol = 0.1;
    assert(abs(generalised_hurst_exponent(X, 1:19, 1) - 0.5) <= tol) 
end

function setup(testCase)  
    testCase.TestData.gaussian_walk = cumsum(randn(10000, 1));
end