function tests = qth_abs_momentTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    assert(qth_abs_moment(testCase.TestData.A, 1, 1) == 1)
end

function testFunctionTwo(testCase)
    assert(qth_abs_moment(testCase.TestData.A, 2, 1) == 2)
end

% TODO - add function checking that it throws w/ tau >= length(X)

function setup(testCase)  
    %create data
    testCase.TestData.A = (1:100)';
end