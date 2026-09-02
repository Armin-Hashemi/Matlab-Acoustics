
% If there's an air gap behind the porous layer before the rigid wall,
% the rigid-wall impedance seen behind the layer becomes -j*cot(k0*d_air)
% and Zs must be found via transfer-matrix / input-impedance recursion
% instead of the simple formula above (left as d_air = 0 case here).

clear
clc
f = linspace(50, 8000, 500);   % [Hz]
omega = 2*pi*f;
% Air properties
rho0   = 1.213;        % air density [kg/m^3]
c0     = 342.2;        % speed of sound [m/s]
% Porous material parameters -----------------
sigma_f = 20000;       % flow resistivity [Pa.s/m^2] (Rayls/m)
D       = 0.05;         % thickness of porous layer [m]
d_air   = 0;             % air gap behind porous layer [m] (0 = layer against rigid wall)

% Delany-Bazley dimensionless frequency parameter
% X = rho0*f/sigma                                    (Delany & Bazley 1970)
X = rho0*f/sigma_f;

% Characteristic impedance of the porous material (normalized by rho0*c0)
% R/(rho0*c0) = 1 + 0.0571*X^-0.754                    (D&B 1970, Fig. 1)
% X/(rho0*c0) =     -0.087*X^-0.732                    (D&B 1970, Fig. 2)
% Zc/(rho0*c0) = R/(rho0*c0) + j*X/(rho0*c0)
zc = 1 + 0.0571*X.^(-0.754) - 1i*0.087*X.^(-0.732);

% Propagation constant (complex wavenumber), normalized by k0 = omega/c0
% alpha_p = (c0/omega)*alpha = 0.189*X^-0.595          (D&B 1970, Fig. 3)
% beta_p  = (c0/omega)*beta  = 1 + 0.098*X^-0.700       (D&B 1970, Fig. 4)
% gamma = alpha + j*beta = (omega/c0)*(alpha_p + j*beta_p)
alpha_p = 0.189*X.^(-0.595);   % attenuation constant term
beta_p  = 1 + 0.098*X.^(-0.700);  % phase constant term
k0 = omega./c0;
gamma = k0.*(alpha_p + 1i*beta_p);  % complex propagation constant

% Surface impedance of porous layer of thickness D backed by rigid wall
% Z = Z0*coth(gamma*l)                  (D&B 1970, "Applications" section)
Zs = zc./tanh(gamma*D);


% Normal-incidence absorption coefficient
% alpha_n = 1 - |(Z-rho0*c0)/(Z+rho0*c0)|^2             (D&B 1970, "Applications" section)
alpha = 1 - abs((Zs - 1)./(Zs + 1)).^2;

figure('Color','w');
semilogx(f, alpha, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Absorption coefficient \alpha');
title('Porous Absorber - Normal Incidence Absorption (Delany-Bazley model)');
ylim([0 1]);
xlim([f(1) f(end)]);
% Annotate peak
[alpha_max, idx_max] = max(alpha);
hold on;
plot(f(idx_max), alpha_max, 'ro', 'MarkerFaceColor','r');
text(f(idx_max), alpha_max-0.05, ...
    sprintf('  f_{peak} = %.0f Hz, \\alpha_{max} = %.3f', f(idx_max), alpha_max));
%% ---------------- Print summary -----------------------------------
fprintf('--- Porous Material Parameters ---\n');
fprintf('Flow resistivity sigma = %.0f Pa.s/m^2\n', sigma_f);
fprintf('Thickness D = %.1f mm\n', D*1e3);
fprintf('--- Result ---\n');
fprintf('Peak absorption alpha_max = %.3f at f = %.1f Hz\n', alpha_max, f(idx_max));
