function tests = generate_expweightsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check that weights are normalised
    assert(ismembertol(sum(testCase.TestData.w), 1, 1e-6))
end

function testFunctionTwo(testCase)
    %check that weights are increasing
    greater_than_prev = arrayfun(@(i) testCase.TestData.w(i) > testCase.TestData.w(i - 1), 2:testCase.TestData.dt);
    assert(sum(greater_than_prev) == testCase.TestData.dt - 1)
end

function testFunctionThree(testCase)
    %check dimensions 
    assert(isequal(size(testCase.TestData.w), [1 testCase.TestData.dt]))
end

function setup(testCase)  
    %create data
    testCase.TestData.dt = 1000;
    testCase.TestData.theta = testCase.TestData.dt / 3;
    testCase.TestData.w = generate_expweights(testCase.TestData.dt, testCase.TestData.theta);
end
