function tests = zeta_estimatorTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    tol = 0.1;
    assert(abs(testCase.TestData.z_g - 0.5) < tol)
end

function setup(testCase)  
    %create data
    testCase.TestData.gaussian_walk = cumsum(randn(10000, 1));
    [testCase.TestData.z_g, testCase.TestData.z_s] = zeta_estimator(testCase.TestData.gaussian_walk, 1:19, 1);
end