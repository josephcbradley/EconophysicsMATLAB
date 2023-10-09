function tests = shuffle_matrix_along_rowsTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    %means should be equal
    assert(isequal(mean(testCase.TestData.X, 2), mean(testCase.TestData.shuffled_X, 2)))
end

function testFunctionTwo(testCase)
    assert(isequal(size(testCase.TestData.X), size(testCase.TestData.shuffled_X)))
end

function setup(testCase)  
    %create data
    testCase.TestData.X = magic(3);
    testCase.TestData.shuffled_X = shuffle_matrix_along_rows(testCase.TestData.X);
end