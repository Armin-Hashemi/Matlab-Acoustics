clear
clc

rho_g = 2500;     % Glass density (kg/m^3)
h = 0.006;        % Glass thickness (m)
m = rho_g*h;      % Surface density (kg/m^2)
c = 340;          % Speed of sound in air (m/s)
rho = 1.2;        % Air density (kg/m^3)
fc = 2079;        % Coincidence frequency (Hz)
f = 50:1:fc/2;    % Frequency range (Hz)

STL_normal = zeros(size(f));     % Normal incidence (x = 0)
STL_diffuse = zeros(size(f));    % Diffuse field (angular integration)

for i = 1:length(f)
    w = 2*pi*f(i);

    % 1) Normal incidence (x = 0)
    tau_normal = 1/(1 + (w*m/(2*c*rho))^2);
    STL_normal(i) = 10*log10(1/tau_normal);

    % 2) Diffuse field — Angular integration
    tau_diffuse = @(x) 1./(1 + (w*m*cos(x)./(2*c*rho)).^2);
    sol = integral(@(x) tau_diffuse(x).*sin(2*x), 0, 78*pi/180);
    STL_diffuse(i) = 10*log10(1/sol);
end

% Plot
semilogx(f, STL_normal, 'b', 'LineWidth', 1.5)
hold on
semilogx(f, STL_diffuse, 'r', 'LineWidth', 1.5)

xlabel('Frequency (Hz)')
ylabel('STL (dB)')
title('Sound Transmission Loss')
legend('Normal incidence', 'Diffuse field')
grid on
axis([50 fc/2 10 50])
