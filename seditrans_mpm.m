%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEDIMENT TRANSPORT INTENSITY (Meyer-Peter & Muller 1948)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phi, dphiD, dphiT ] = seditrans_mpm(theta, flag)

switch flag
    
% original formula
    case 1
        thetacr = 0.047;
        phi0 = 8;
        A = 1.5;
      
% modified formula
    case 2
        thetacr = 0.0495;
        phi0 = 3.97;
        A = 1.5;

% error
    otherwise
        error('Wrong flag for Meyer-Peter Muller formula!');
end
        

% sediment transport intensity
phi = phi0 * (max(0, theta-thetacr))^A;

% derivatives
dphiD = 0;
dphiT = A * max(0, phi/(theta-thetacr));

% end of function
return