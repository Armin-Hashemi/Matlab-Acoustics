% porous_absorption_JCA.m
% Sound absorption coefficient of a rigid-backed porous layer using the
% full Johnson-Champoux-Allard (JCA) model with porosity, tortuosity,
% and viscous/thermal characteristic lengths.
%
% Source: Cox, T.J. and D'Antonio, P., "Acoustic Absorbers and
% Diffusers: Theory, Design and Application," 3rd ed., CRC Press,
% Chapter 6, 6.5.3 Semi-phenomenological models (Eqs. 6.19-6.25), p. 206.
%
% Author: Armin Hashemii
clear
clc
f = linspace(50, 8000, 500);   % [Hz]
omega = 2*pi*f;
% Air properties
rho0   = 1.213;        % air density [kg/m^3]
c0     = 343.0;         % speed of sound [m/s]
eta    = 1.84e-5;       % dynamic viscosity of air [Pa.s]
kappa  = 2.41e-2;       % thermal conductivity of air [W/(m.K)]
cp     = 1.01e3;        % specific heat capacity of air at const. pressure [J/(kg.K)]
gam    = 1.4;           % ratio of specific heats
P0     = 101320;        % atmospheric pressure [N/m^2]

% Porous material (JCA) parameters -----------------
sigma_f  = 20000;       % flow resistivity [Pa.s/m^2]
phi      = 0.98;        % porosity
alpha_inf = 1.05;       % tortuosity
Lambda   = 100e-6;      % viscous characteristic length [m]
Lambda_p = 200e-6;      % thermal characteristic length [m]
D        = 0.05;        % thickness of porous layer [m]

% Viscous boundary layer thickness
% delta_v = sqrt(2*eta/(rho0*omega))                        (Eq. 6.22)
delta_v = sqrt(2*eta./(rho0*omega));

% Thermal boundary layer thickness
% delta_h = sqrt(2*kappa/(rho0*cp*omega))                   (Eq. 6.23)
delta_h = sqrt(2*kappa./(rho0*cp*omega));

% Prandtl number
% Np = (delta_v/delta_h)^2                                  (Eq. 6.21)
Np = (delta_v./delta_h).^2;

% Effective (dynamic) density of the porous material
% rho_e = (alpha_inf*rho0/phi) * [1 + (sigma*phi)/(j*omega*rho0*alpha_inf)
%          * sqrt(1 + 4j*alpha_inf^2*eta*rho0*omega/(sigma^2*Lambda^2*phi^2))]
%                                                             (Eq. 6.19)
rho_e = (alpha_inf*rho0./phi) .* ...
    (1 + (sigma_f*phi)./(1i*omega*rho0*alpha_inf) .* ...
    sqrt(1 + (4i*alpha_inf^2*eta*rho0.*omega)/(sigma_f^2*Lambda^2*phi^2)));

% Effective (dynamic) bulk modulus of air in the material
% Ke = (gamma*P0/phi) * { gamma - (gamma-1)*[1 + (8*eta)/(j*Lambda'^2*Np*omega*rho0)
%       * sqrt(1 + j*rho0*omega*Np*Lambda'^2/(16*eta))]^-1 }^-1
%                                                             (Eq. 6.20)
inner_sqrt = sqrt(1 + (1i*rho0.*omega.*Np*Lambda_p^2)/(16*eta));
bracket = 1 + (8*eta)./(1i*Lambda_p^2.*Np.*omega*rho0) .* inner_sqrt;
K_e = (gam*P0/phi) .* (gam - (gam-1)*bracket.^(-1)).^(-1);

% Characteristic impedance of the material
% zc = sqrt(Ke*rho_e)                                        (Eq. 6.24)
zc = sqrt(K_e.*rho_e);

% Propagation wavenumber
% k = omega*sqrt(rho_e/Ke)                                   (Eq. 6.25)
k = omega.*sqrt(rho_e./K_e);

% Surface impedance of rigidly-backed layer of thickness D
Zs = -1i*zc.*cot(k*D);

% Normal-incidence energy absorption coefficient
alpha = 1 - abs((Zs - rho0*c0)./(Zs + rho0*c0)).^2;

figure('Color','w');
semilogx(f, alpha, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Absorption coefficient \alpha');
title('Porous Absorber - Normal Incidence Absorption (JCA model)');
ylim([0 1]);
xlim([f(1) f(end)]);
% Annotate peak
[alpha_max, idx_max] = max(alpha);
hold on;
plot(f(idx_max), alpha_max, 'ro', 'MarkerFaceColor','r');
text(f(idx_max), alpha_max-0.05, ...
    sprintf('  f_{peak} = %.0f Hz, \\alpha_{max} = %.3f', f(idx_max), alpha_max));
%% ---------------- Print summary -----------------------------------
fprintf('--- JCA Material Parameters ---\n');
fprintf('Flow resistivity sigma = %.0f Pa.s/m^2\n', sigma_f);
fprintf('Porosity phi = %.2f\n', phi);
fprintf('Tortuosity alpha_inf = %.2f\n', alpha_inf);
fprintf('Lambda = %.1f um, Lambda_p = %.1f um\n', Lambda*1e6, Lambda_p*1e6);
fprintf('Thickness D = %.1f mm\n', D*1e3);
fprintf('--- Result ---\n');
fprintf('Peak absorption alpha_max = %.3f at f = %.1f Hz\n', alpha_max, f(idx_max));
