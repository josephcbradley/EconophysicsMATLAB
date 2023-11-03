function tests = force_symmetricTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check output dimensions
    assert(issymmetric(force_symmetric(testCase.TestData.X)))
end

function setup(testCase)  
    testCase.TestData.X = magic(5);
end