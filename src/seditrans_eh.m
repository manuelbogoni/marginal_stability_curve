%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEDIMENT TRANSPORT INTENSITY by Engelund and Hansen (1967)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phi, dphiD, dphiT ] = seditrans_eh(theta, Cf, dCD, dCT)

% constants
A = 0.05;
B = 2.50;

% sediment transport intensity
phi  = A/Cf * theta^B;

% derivatives
if nargin >= 3
    dphiD = -phi/Cf * dCD;
    dphiT = -phi/Cf * dCT + B/theta*phi;
else
    dphiD = 0;
    dphiT = 0;
end
      
% end of the function
return