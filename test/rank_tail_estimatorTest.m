function tests = rank_tail_estimatorTest
    tests = functiontests(localfunctions);
end

function testFunctionOne(testCase)
    data = testCase.TestData.data;
    alpha = testCase.TestData.alpha; 
    assert(ismembertol(rank_tail_estimator(data), alpha, 1e-1)); %slow to converge!
end

function setup(testCase)  
    %create data
    %For MATLAB generalised Pareto, xmin = sigma/k and alpha = 1/k, and 
    %theta = sigma/k
    %thus if alpha = 2.5 and xmin = 1, k = 0.2, sigma = 0.2, and theta = 
    % 1
    %gppdf(x,k,sigma,theta) 
    xmin = 1;
    alpha = 2.5;
    k = 1/alpha;
    sigma = xmin * k;
    theta = sigma / k;
    testCase.TestData.alpha = alpha;
    N = 1000;
    testCase.TestData.data = gprnd(k,sigma,theta, N, 1);
end