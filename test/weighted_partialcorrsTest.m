function tests = weighted_partialcorrsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    assert(isequal(size(testCase.TestData.partial_C), [testCase.TestData.N, testCase.TestData.N]))
end

function testFunctionTwo(testCase)
    rho = testCase.TestData.partial_C(1, 2);
    assert(rho < 0.1)
end

function testFunctionThree(testCase)
    rho = testCase.TestData.raw_C(1, 2);
    assert(rho > 0.1)
end



function setup(testCase)  
    %create data
    testCase.TestData.dt = 100000;
    testCase.TestData.theta = testCase.TestData.dt / 3;
    testCase.TestData.N = 2;
    %testCase.TestData.X = randn(testCase.TestData.dt, testCase.TestData.N);
    testCase.TestData.w = generate_expweights(testCase.TestData.dt, testCase.TestData.theta);
    %testCase.TestData.C = weightedcorrs(testCase.TestData.X, testCase.TestData.w);
    testCase.TestData.Y = randn(testCase.TestData.dt, 1);
    testCase.TestData.A = 0.1 + (0.3 * testCase.TestData.Y) + randn(testCase.TestData.dt, 1);
    testCase.TestData.B = 0.2 + (1.2 * testCase.TestData.Y) + randn(testCase.TestData.dt, 1);
    testCase.TestData.partial_C = weighted_partialcorrs([testCase.TestData.A testCase.TestData.B], testCase.TestData.Y, testCase.TestData.w);
    testCase.TestData.raw_C = weightedcorrs([testCase.TestData.A testCase.TestData.B], testCase.TestData.w);
end