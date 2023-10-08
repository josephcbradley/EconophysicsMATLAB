function tests = weightedcorrsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check that C is symmetric 
    assert(issymmetric(testCase.TestData.C))
end

function testFunctionTwo(testCase)
    %check that C is positive definite 
    [~, flag] = chol(testCase.TestData.C);
    assert(flag == 0)
end

function testFunctionThree(testCase)
    %check that C is NxN
    N = testCase.TestData.N;
    assert(isequal(size(testCase.TestData.C), [N N]))
end

function setup(testCase)  
    %create data
    testCase.TestData.dt = 1000;
    testCase.TestData.theta = testCase.TestData.dt / 3;
    testCase.TestData.N = 100;
    testCase.TestData.X = randn(testCase.TestData.dt, testCase.TestData.N);
    testCase.TestData.w = generate_expweights(testCase.TestData.dt, testCase.TestData.theta);
    testCase.TestData.C = weightedcorrs(testCase.TestData.X, testCase.TestData.w);
end