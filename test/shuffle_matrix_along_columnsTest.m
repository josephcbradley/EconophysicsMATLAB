function tests = shuffle_matrix_along_columnsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %means should be equal
    assert(isequal(mean(testCase.TestData.X, 1), mean(testCase.TestData.shuffled_X, 1)))
end

function testFunctionTwo(testCase)
    assert(isequal(size(testCase.TestData.X), size(testCase.TestData.shuffled_X)))
end

function setup(testCase)  
    %create data
    testCase.TestData.X = magic(3);
    testCase.TestData.shuffled_X = shuffle_matrix_along_columns(testCase.TestData.X);
end