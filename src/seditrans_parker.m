%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEDIMENT TRANSPORT INTENSITY (Parker 1982, 1990)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phi, dphiD, dphiT ] = seditrans_parker(theta)

thetar = 0.0386;
csi = theta / thetar;
A = 0.00218;

if csi < 1
	B = 14.2;
	G0 = A * (theta / thetar)^B;
    dG0 = B/theta * G0;
elseif csi >= 1 && csi <= 1.65
	B = 14.2;
	C = 9.28;
	G0 = A * exp(B * (theta/thetar-1) - C * (theta/thetar-1)^2);
    dG0 = G0/thetar * (B + 2*C - 2*C*theta/thetar);
elseif csi > 1.65
	B = 5474;
	C = 0.853;
	E = 4.5;
	G0 = A * B * (1 - C*thetar/theta)^E;
    dG0 = A*B*C*E*thetar/theta^2 * (1 - C*thetar/theta)^(E-1);
end

% sediment transport intensity
a = 1.5;
phi = G0 * theta^a;

% derivatives
dphiD = 0;
dphiT = dG0 * theta^a + a * phi/theta;
      
% end of function
return