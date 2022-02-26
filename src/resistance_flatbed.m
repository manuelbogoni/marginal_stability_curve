%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRICTION RESISTANCE FOR FLAT BED (Keulegan 1938, Einstein 1950)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Cf, dCD, dCT ] = resistance_flatbed(dg,D)

if nargin == 1
    RelRough = dg;
else
    RelRough = dg./D;
end

% friction coeffient
Cf = ( 6 + 2.5 * log(1./(2.5*RelRough)) ).^(-2);

% derivatives
dCD = -5 * ( 6 + 2.5 * log(1./(2.5*RelRough)) ).^(-3);
dCT = 0;

% end of function
return