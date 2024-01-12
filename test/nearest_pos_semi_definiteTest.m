function tests = nearest_pos_semi_definiteTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %check that X is actually different
    assert(~isequal(testCase.TestData.B, testCase.TestData.X))
end

function testFunctionTwo(testCase)
    %check that X is actually pos semi def
    l = min(eig(testCase.TestData.X));
    assert((l >= 0) || (abs(l) < 1e-10))
end

function testFunctionThree(testCase)
    %check that X is real
    assert(isreal(testCase.TestData.X))
end

function setup(testCase)  
    %create data
    testCase.TestData.B = magic(5);
    testCase.TestData.B = force_symmetric(testCase.TestData.B);
    testCase.TestData.X = nearest_pos_semi_definite(testCase.TestData.B);
end