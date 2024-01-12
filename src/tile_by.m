function tiling = tile_by(data, n_tiles)
%% Description 
%Given an TxN matrix of features, tile cross sectionally
%% Inputs
%data - TxN data matrix
%n_tiles - number of tiles to be sorted into
%% Ouputs 
%tiling - finished tiling
%% Setup 
if n_tiles < 2
    error("n_tiles < 2 not supported")
else
    p = arrayfun(@(i) i/(n_tiles), 1:n_tiles-1);
end

%% Calculation
Q = quantile(data, p, 2);
tiling = NaN(size(data));
[T, N] = size(data);
for i = 1:N
    for t = 1:T
        if isnan(data(t, i))
            continue
        else
            tiling(t, i) = sum(data(t, i) > Q(t, :)) + 1;
        end

    end
end



end


