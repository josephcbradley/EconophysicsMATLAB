function must_be_numeric(A)
%% Description 
%Checks that A is numeric. Throws if not

%% Inputs
% A - any numeric array
arguments 
    A (:, :)
end
%% Ouputs 
% None
%% Setup 
% None

%% Calculation
if ~isnumeric(A)
    eid = 'Class:NotNumeric';
    msg = 'A is not numeric.';
    error(eid, msg)
end

end