%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRICTION RESISTANCE FOR DUNE-COVERED BED (Engelund-Hansen 1967)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Cf, dCD, dCT, rpic ] = resistance_dunebed_eh(ds, theta, rpic0)

% experimental function
if theta <= 2.4385
	thetaA = 0.06 + 0.4 * theta.^2;
	dthetaA = 0.8 * theta;
else
    thetaA = theta;
    dthetaA = 1;
end

% friction coeffient Cf
Cf = (6 + 2.5 * log((thetaA./theta)./(2.5*ds))).^(-2) .* (thetaA./theta).^(-1);

% derivative of the friction coefficient with respect to the scaled depth
dCD = -5 * (thetaA./theta).^(-1) .* ...
    (6 + 2.5 * log((thetaA./theta)./(2.5*ds))).^(-3);

% derivative of the friction coefficient with respect to the Shields number
% dCT = dCf1 * Cf2 + Cf1 * dCf2 where dCf1 = d(Cf1) and dCf2 = d(Cf2)
dCT = Cf ./ theta .* (1-theta./thetaA .* dthetaA) .* ...
    (1+5*(6 + 2.5 * log((thetaA./theta)./(2.5*ds))).^(-1));

% Talmon parameter for transverse bed slope
if nargin == 3
    rpic = rpic0./sqrt(thetaA./theta);
else
    rpic = [];
end

% check wheter friction is lower than the flat scenario
for j = 1:length(Cf)
    if Cf(j) < ((6 + 2.5 * log(1./(2.5*ds(j))) )^(-2))
        [ Cf(j), dCD(j), dCT(j) ] = resistance_flatbed(ds(j));
        if nargin == 3
            rpic(j) = rpic0(j);
        else
            rpic = [];
        end
    end
end

% end of function
return
