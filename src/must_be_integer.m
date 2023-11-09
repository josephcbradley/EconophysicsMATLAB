function flag = must_be_integer(x)
    flag = floor(x) == x;
    %no guarantee this works for Inf, NaN, etc...
end