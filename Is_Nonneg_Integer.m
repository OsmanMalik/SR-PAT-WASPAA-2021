function retBool = Is_Nonneg_Integer(x)
%Is_Positive_Integer    Check if a value is a positive integer
%
% INPUTS:
%   x           - The value to be checked.
%
% OUTPUTS:
%   retBool     - Return value is "true" if x is a nonnegative integer, 
%                 false otherwise.

if isreal(x) && mod(x,1) == 0 && x >= 0
    retBool = true;
else
    retBool = false;
end

end
