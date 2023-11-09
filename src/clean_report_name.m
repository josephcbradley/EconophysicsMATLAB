function new_name = clean_report_name(old_name)
    new_name = replace(old_name, {'-', ' '}, '_');
    new_name = lower(new_name);
end

