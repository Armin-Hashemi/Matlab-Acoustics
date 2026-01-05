% Sound Transmission Loss with SEA method
% All the information is available in the Sound Insulation by Carl Hopkins.
clear; clc;
cl = 5200;
Lx = 1.5;                         % Panel length (m)
Ly = 1.25;                        % Panel height (m)
S = Lx * Ly;                      % Panel Area (m^2)
h   = 0.006;                      % Panel thickness (m)
E   = 6.37e10;                    % Young’s modulus of the panel material (Pa)
rho = 1.2;                        % Air density (kg/m^3)
nu  = 0.24;                       % Poisson’s ratio of the panel material (-)
c   = 343;                        % Speed of sound in air (m/s)
rhos2 = 2500*h;                   % Surface density of the panel (kg/m^2)
B = E*h^3/(12*(1 - nu^2));        % Bending stiffness of the plate (N·m)
fc = c^2/(2*pi)*sqrt(rhos2/B);    % Coincidence frequency (Hz). Eq(2.201)

% USER-DEFINED PARAMETERS (valuses are ommited during subsitution)
V1 = 1;      % Room volume 1
V3 = 1;      % Room volume 3
T3 = 1;      % Reverberation time
eta1 = 0.01;
eta2 = 0.025;

% Frequency range
f_start = 100;
f_stop  = 6000;
n_max = floor(6*log2(f_stop / f_start));      % number of 1/6-octave steps
f = f_start * 2.^((0:n_max)/6);               % center frequencies (Hz)
TL = zeros(size(f));

for i = 1:length(f)
    % Wavenumber and mu
    k = 2*pi*f(i)/c;
    mu = sqrt(fc/f(i));

    % Radiation efficiency sigma for f<fc
    if f(i) < fc
        Cbc = 1;
        Cob = 1;
        sigma = 2*(Lx + Ly) / (2*pi*mu*k*S*sqrt(mu^2 - 1)) * ...
            ( log((mu + 1)/(mu - 1)) + 2*mu/(mu^2 - 1) ) *( Cbc*Cob - (Cbc*Cob - 1)/mu^8 );

        % U term for non-resonant coupling loss factors. Eq. (4.44)
        U = 1/(2*pi) * (Lx/Ly + Ly/Lx) * log(1 + (Lx/Ly)^2) ...
            - (0.5 + Lx/(pi*Ly)) * log(Lx/Ly) - log(2)/pi - 2/pi * integral(@(t) atan(t)./t, Lx/Ly, 1);
        % Non resonant transmission coefficient. Eq. (4.42)
        tau = (2*rho/(rhos2*k*(1 - 1/mu^4)))^2 * ...
            ( log(k*sqrt(S)) + 0.16 - U + 1/(4*mu^6) * ( ...
            (2*mu^2 - 1)*(mu^2 + 1)^2 * log(mu^2 - 1) + ...
            (2*mu^2 + 1)*(mu^2 - 1)^2 * log(mu^2 + 1) - 4*mu^2 - 8*mu^6*log(mu) ) );
    else
        sigma=1/sqrt(1-mu^2);
        U=0;
        tau=0;
    end

    % Modal densities
    n1 = 4*pi*f(i)^2 * V1 / c^3; % Eq. (1.59)
    n2 = sqrt(3)*S / (h*cl);
    n3 = 4*pi*f(i)^2 * V3 / c^3;

    % Non resonant Coupling loss factors. Eq. (4.25)
    eta13 = c*S/(8*pi*f(i)*V1) * tau;
    eta31 = c*S/(8*pi*f(i)*V3) * tau;

    % Resonant Coupling loss factors. Eq. (4.21)
    eta21 = rho*c*sigma/(2*pi*f(i)*rhos2);
    eta12 = eta21 * n2 / n1;
    eta23 = rho*c*sigma/(2*pi*f(i)*rhos2);
    eta32 = n2 * eta23 / n3;
    % Recieved room internal loss factor Eq. (1.107)
    eta3 = 2.2/(f(i)*T3);
    % Sound Transmission Loss
    if f(i) < fc
        arg = (eta3*V3*S*T3/(0.16*V3)/(eta13*V1));
    else
        arg = eta2 .* eta3 .* V3 .* S .* T3 ./ (0.16*V3) ./ (eta12 .* eta23 .* V1);
    end
    TL(i) = 10*log10(arg);
end

% Plot
semilogx(f, TL, 'LineWidth', 1.5)
xlabel('Frequency (Hz)')
ylabel('STL (dB)')
title('Sound Transmission Loss vs Frequency')
grid on