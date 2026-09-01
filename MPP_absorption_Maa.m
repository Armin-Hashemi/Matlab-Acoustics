% Reference papers used for verification:
% 1. Maa, D.-Y. (1975). "Theory and Design of Microperforated Panel
%    Sound-Absorbing Constructions." Scientia Sinica, XVIII, 55-71.
% 2. Maa, D.-Y. (1989). "Microperforated-Panel Wideband Absorbers."
%    Noise Control Engineering Journal, 29(3), 77-84.
% 3. Maa, D.-Y. (1998). "Potential of Microperforated Panel Absorber."
%    J. Acoust. Soc. Am., 104(5), 2861-2866.
% 
% Note: papers 1 and 2 agree with each other on every coefficient below.
% Paper 3 (1998) contains two apparent typos relative to both. The script
% follows the 1975/1989 (majority, mutually consistent) versions.
% 
% Errors fixed:
% -------------
% 1. Perforate constant k (Maa Eq. 6)
%    Wrong:   k = d*sqrt(omega*rho0/(eta/4))
%    Correct: k = d*sqrt(omega*rho0/(4*eta))
%    The eta/4 vs 4*eta swap made k 4x too large (16x too large inside
%    the sqrt). This was a typing error, not a paper discrepancy.
% 
% 2. Resistance coefficient k_r (Maa 1975 Eq.14 / 1989 Eq.10)
%    Wrong:   k_r = sqrt(1+k^2/32) + (sqrt(2)/32)*k*(d/t)
%    Correct: k_r = sqrt(1+k^2/32) + (sqrt(2)/8)*k*(d/t)
%    The 1998 paper prints sqrt(2)/32; both 1975 and 1989 print sqrt(2)/8.
%    Treated as a typo in the 1998 paper.
% 
% 3. Mass reactance coefficient k_m (Maa 1975 Eq.15 / 1989 Eq.11)
%    Wrong:   k_m = 1 + (1+k^2/2)^(-1/2) + 0.85*(d/t)
%    Correct: k_m = 1 + 1/sqrt(9+k^2/2) + 0.85*(d/t)
%    The 1998 paper's denominator is missing the "9" (prints 1+k^2/2
%    instead of 9+k^2/2); both 1975 and 1989 have the "9".

clear
clc
f = linspace(50, 8000, 500);   % [Hz]
omega = 2*pi*f;
% Air properties -----------------
rho0   = 1.213;        % air density [kg/m^3]
c0     = 342.2;        % speed of sound [m/s]
eta    = 1.789e-5;     % dynamic viscosity of air [Pa*s]

% MPP parameters -----------------
t      = 0.5e-3;       % panel thickness [m]
d      = 0.5e-3;       % hole diameter [m]
sigma  = 1E-02;        % perforation ratio
D      = 0.06;         % air cavity depth behind panel [m]
k = d .* sqrt(omega*rho0/(4*eta));   % (Maa Eq. 6)
k_r = sqrt(1 + k.^2/32) + (sqrt(2)/8).*k.*(d/t); % (Maa Eq. 5a)  <-- FIXED: was sqrt(2)/32, should be sqrt(2)/8

% Normalized resistance r 
r = (32*eta*t)./(sigma*rho0*c0*d^2) .* k_r;  % (Maa Eq. 5a)
k_m = 1 + 1./sqrt(9 + k.^2/2) + 0.85*(d/t); % (Maa Eq. 5b)  <-- FIXED: was (1+k^2/2)^-1/2, should be 1/sqrt(9+k^2/2)
omega_m = (omega*t)./(sigma*c0) .* k_m; % (Maa Eq. 5b)

% Cavity reactance term
cavity_term = cot(omega.*D./c0);

% Normal-incidence absorption coefficient -------
alpha = 4*r ./ ( (1+r).^2 + (omega_m - cavity_term).^2 );  % (Maa Eq. 9)

% Plot -----------------
figure('Color','w');
semilogx(f, alpha, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Absorption coefficient \alpha');
title('MPP Absorber - Normal Incidence Absorption (Maa model)');
ylim([0 1]);
xlim([f(1) f(end)]);
% Annotate peak
[alpha_max, idx_max] = max(alpha);
hold on;
plot(f(idx_max), alpha_max, 'ro', 'MarkerFaceColor','r');
text(f(idx_max), alpha_max-0.05, ...
    sprintf('  f_{peak} = %.0f Hz, \\alpha_{max} = %.3f', f(idx_max), alpha_max));
%% ---------------- Print summary -----------------------------------
fprintf('--- MPP Parameters ---\n');
fprintf('Thickness t     = %.3f mm\n', t*1e3);
fprintf('Hole diameter d = %.3f mm\n', d*1e3);
fprintf('Perforation ratio sigma = %.2f %%\n', sigma*100);
fprintf('Cavity depth D  = %.1f mm\n', D*1e3);
fprintf('--- Result ---\n');
fprintf('Peak absorption alpha_max = %.3f at f = %.1f Hz\n', alpha_max, f(idx_max));
