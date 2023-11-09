function must_be_real(A)
%% Description 
%Checks that A is real. Throws if not

%% Inputs
% A - any numeric array
arguments 
    A (:, :) {must_be_numeric}
end
%% Ouputs 
% None
%% Setup 
% None

%% Calculation
if ~isreal(A) | isnan(A) | isinf(A)
    eid = 'Class:NotReal';
    msg = 'A is not real.';
    error(eid, msg)
end
end