function output = factor_mimicking_portfolio(factor, returns, n_tiles)
%% Description 

%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
tiles = tile_by(factor, n_tiles);
output = NaN(size(returns));
output(tiles == n_tiles) = 1;
output(tiles == 1) = -1;
output = mean(output .* returns, 2, 'omitnan');

end