function tests = tile_byTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    tiling = tile_by(testCase.TestData.X, testCase.TestData.n_tiles);
    correct_X = [1 1 2 2; 1 1 2 2; 1 1 2 2];
    assert(isequal(correct_X, tiling))
end

function testFunctionTwo(testCase)
    tiling = tile_by(testCase.TestData.Y, testCase.TestData.n_tiles);
    assert(isnan(tiling(3, 2)) && isnan(tiling(3, 3)))
end



function setup(testCase)  
    %create data
    testCase.TestData.X = [1 2 3 4; 5 6 7 8; 9 10 11 12];
    testCase.TestData.Y = [1 2 3 4; 5 6 7 8; 9 NaN NaN 12];
    testCase.TestData.n_tiles = 2;
end