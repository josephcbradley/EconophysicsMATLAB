function tests = lsq_estimatorTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    beta = lsq_estimator(testCase.TestData.y, testCase.TestData.x);
    tol = 0.001;
    assert(abs(testCase.TestData.beta_matrix - beta) < tol)
end

function setup(testCase)  
    %create data
    testCase.TestData.x = (-1:0.1:1)';
    testCase.TestData.y = 2*testCase.TestData.x + 3 + randn(length(testCase.TestData.x), 1);
    testCase.TestData.beta_matrix = (testCase.TestData.x' * testCase.TestData.x) \ (testCase.TestData.x' * testCase.TestData.y); 
end