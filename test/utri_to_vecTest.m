function tests = utri_to_vecTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    assert(isequal(utri_to_vec(testCase.TestData.X), [1; 6; 7]))
end

function setup(testCase)  
    %create data
    testCase.TestData.X = magic(3);
end