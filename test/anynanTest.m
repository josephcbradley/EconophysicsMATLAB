function tests = anynanTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    assert(~anynan(testCase.TestData.X))
end

function testFunctionTwo(testCase)
    testCase.TestData.X(3) = NaN;
    assert(anynan(testCase.TestData.X))
end

function setup(testCase)  
    %create data
    testCase.TestData.X = [1 2 3];
end