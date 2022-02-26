%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARGINAL STABILITY CURVE LAMBDA-BETA FOR FINITE AMPLITUDE ALTERNATE BARS
% Colombini et al. JFM 1987
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the code
clear; close all; clc

%% INPUT DATA

% set constants
rho = 1000;                 % specific weight of water (kg/m3)
rhos = 2650;                % specific weight of sediment (kg/m3)
g = 9.806;                  % gravity acceleration (m/s2)
nu = 1.01*10^(-6);          % kinematic viscosity (m2/s)
p = 0.3;                    % porosity of the bed (-)

% set bed configuration
flag_bed = input('Flat bed (1) or dune-covered bed (2) ? ');

% set energy slope
S = input('Energy slope (-)  S = ');

% set flow depth
D0 = input('Depth (m)        D0 = ');

% set characteristic sediment size
dg = input('Grain size (mm)  dg = ');
dg = dg/1000;

% set Talmon parameter
rpic = 0.3;

%% DIMENSIONLESS PARAMETERS

tic
fprintf('----------------------------------------------------\n');

% specific gravity
Delta = rhos/rho-1;

% Shields number
theta0 = D0*S / (Delta * dg);
fprintf('Shields number             tau* = %9.5f \n', theta0);

% dimensionless grain size
ds = dg / D0;
fprintf('scaled grain size          ds   = %9.5f \n', ds);

% Reynolds particle number
Rp = sqrt(Delta*g) * dg^1.5 / nu;
fprintf('Reynolds particle number   Rp   = %9.5f \n', Rp);

%% FLOW CHARACTERISTICS

% compute the friction coefficient, the sediment transport intensity and
% their derivatives
if flag_bed == 1
    [ Cf0, dCD, dCT ] = resistance_flatbed(ds);
    [ phi0, dphiD, dphiT ] = seditrans_mpm(theta0, 1);
    %[ phi0, dphiD, dphiT ] = seditrans_parker(theta0);
elseif flag_bed == 2
    [ Cf0, dCD, dCT, rpic ] = resistance_dunebed_eh(ds, theta0, rpic);
    [ phi0, dphiD, dphiT ] = seditrans_eh(theta0, Cf0, dCD, dCT);
end
CD = dCD / Cf0;
CT = theta0 / Cf0 * dCT;
phiD = dphiD / phi0;
phiT = theta0 / phi0 * dphiT;

% compute Colombini's terms
s1 = 2/(1-CT);
s2 = CD/(1-CT);
f1 = 2*phiT/(1-CT);
f2 = phiD + CD*phiT / (1-CT);
R = rpic / sqrt(theta0);

fprintf('Friction coefficient       Cf   = %9.5f \n', Cf0);
fprintf('Sediment load intensity    phi  = %9.5f \n', phi0);

% Froude number
F0 = sqrt(Delta * ds * theta0 / Cf0);
fprintf('Froude number              Fr   = %9.5f \n', F0);

% uniform flow velocity
U0 = sqrt(g*D0*S/Cf0);
fprintf('Uniform flow velocity      U    = %9.5f m/s\n', U0);

% specific discharge
q = U0 * D0;
fprintf('Specific discharge         q    = %9.5f m2/s\n', q);

% specific sediment flux
qs = phi0 * Rp * nu;
qskg = qs * rhos;
fprintf('Specific sediment flux     qs   = %9.5f kg/sm\n', qskg);

% term Q0
Q0 = sqrt(Delta*g) * dg^1.5 / ((1-p)*D0*U0);

%% STABILITY MARGINAL CURVE

% initialize storage
beta_v = [];
lambda_v = [];

% set the values of beta
beta = (1:0.020:30)';
nbeta = length(beta);

% set the values of lambda
lambda = [0.01:0.0001:1 1:0.005:5]';
nlambda = length(lambda);

% initialize the counter
jzero = 0;

% loop over the values of beta
for jbeta = 1:nbeta
        
% compute the matrix coefficients
    a11 = 1i * lambda + beta(jbeta)*Cf0*s1;
    a13 = beta(jbeta) * Cf0 * (s2-1);
    a14 = 1i * lambda;
    a22 = 1i* lambda + beta(jbeta)*Cf0;
    a24 = pi/2;
    a31 = 1i * lambda;
    a32 = -pi/2;
    a33 = 1i * lambda;
    a41 = 1i * lambda * Q0 * phi0 * f1;
    a42 = -pi/2 * Q0 * phi0;
    b1 = Q0 * phi0 * (1i*lambda*f2 - pi^2 / 4 * R/beta(jbeta));
    b2 = F0^2 * Q0 * phi0 * pi^2 / 4 * R / beta(jbeta);
        
% compute mu = Omega - i*omega
    DEN = a11.*(a22.*a33*F0^2 - a24.*a32) - a31.*a22.*(a13.*F0^2 + a14);
    NUM = a11.*a24.*(a33.*a42-a32.*b1) + a22.*b2.*(a13.*a31-a11.*a33) ...
        - a13.*a24.*(a31.*a42-a32.*a41) - a14.*a22.*(a31.*b1-a33.*a41);
    mu = NUM./DEN;
        
% since mu is a downbound function (i.e. the concavity is downward), there
% may be zero, one, or two zero-crossing points

% find the first zero crossing point, i.e. the point where the
% function changes from negative to positive values
    idx1 = find(real(mu)>0, 1,'first');
    if ~isempty(idx1) && idx1 > 1
            
% weigh the values to find the exact zero-crossing
        lambdaZero = ...
            (lambda(idx1)*abs(lambda(idx1-1)) + ...
            lambda(idx1-1)*abs(lambda(idx1))) / ...
             (abs(lambda(idx1-1)) + abs(lambda(idx1)));

% update the counter and store the current values
        jzero = jzero + 1;
        aux = [lambda_v; lambdaZero];
        lambda_v = aux;
        aux = [beta_v; beta(jbeta)];
        beta_v = aux;
    end
        
% find the last zero crossing point, i.e. the point where the
% function changes from positive to negative values
    idx2 = find(real(mu)>0, 1,'last');
    if ~isempty(idx2) && idx2 < nlambda

% weigh the values to find the exact zero-crossing
        lambdaZero = ...
            (lambda(idx2)*abs(lambda(idx2+1)) + ...
            lambda(idx2+1)*abs(lambda(idx2))) / ...
             (abs(lambda(idx2+1)) + abs(lambda(idx2)));

% update the counter and store the current values
        jzero = jzero + 1;
        aux = [lambda_v; lambdaZero];
        lambda_v = aux;
        aux = [beta_v; beta(jbeta)];
        beta_v = aux;
    end        
        
% end of the loop over betas
end

% print on screen
fprintf('----------------------------------------------------\n');
if jzero
    fprintf('Number of found zero-crossing points : %u \n', jzero);
else
    fprintf('No zero-crossing found! \n')
end

% sort the marginal stability curve
[lambda_v, idx3] = sort(lambda_v, 'ascend');
beta_v = beta_v(idx3);

% find the critical value
[beta_c, idx4] = min(beta_v);
lambda_c = lambda_v(idx4);

toc;

%% PLOT

% new figure
figure1 = figure('units', 'normalized', 'position', [0.1 0.3 0.8 0.6]);

% plot the marginal stability curve
plot(lambda_v, beta_v, 'k', 'linewidth', 2)
hold on

% get limits
xl = xlim;
yl = [min(beta) max(beta)];

% plot the critical value
plot([xl(1); lambda_c], [beta_c; beta_c], 'r--')
plot([lambda_c; lambda_c], [yl(1); beta_c], 'r--')
plot(lambda_c, beta_c, 'marker', 'o', 'markersize', 6, ...
    'markeredgecolor', 'k', 'markerfacecolor', 'r');
text(lambda_c, beta_c*1.06, ...
    strcat('(',num2str(lambda_c,'%3.2f'),', ', num2str(beta_c,'%3.2f'),')'), ...
    'horizontalalignment', 'left', 'color', 'r');

% set axes
axis square
ylim(yl)
set(gca, 'fontname', 'times new roman', 'fontsize', 15);
xlabel('wavenumber $\lambda = 2\pi / L_*$', 'interpreter', 'latex', ...
    'fontname', 'times new roman', 'fontsize', 20);
ylabel('aspect ratio $\beta = B_0^* / D_0^*$', 'interpreter', 'latex', ...
    'fontname', 'times new roman', 'fontsize', 20);

% title
if flag_bed == 1
    bed_type = 'flat bed';
elseif flag_bed == 2
    bed_type = 'dune-covered bed';
end
str1 = 'Marginal stability curve'; 
str2 = strcat(' (', bed_type, ...
 	', $\theta$ = ', num2str(theta0,'%4.3f'), ...
 	', $d_s$ = ', num2str(ds,'%3.2e'), ')');
% title({str1, str2}, 'interpreter', 'latex', ...
%     'fontname', 'times new roman', 'fontsize', 25);
title(strcat(str1, str2), 'interpreter', 'latex', ...
    'fontname', 'times new roman', 'fontsize', 25);