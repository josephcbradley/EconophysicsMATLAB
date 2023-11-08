function val = raw_ES_from_upper(Ta_upper, Tb_upper)
    %calculate shared edges in Ta and Tb
    val = sum((Ta_upper ~= 0) & (Tb_upper ~= 0));
end