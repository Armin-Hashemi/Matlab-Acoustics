% Sound Transmission Loss with Cremer formula
clear
clc

% Panel and material properties
m   = 15;        % Surface density of the panel (kg/m^2)
h   = 0.006;     % Panel thickness (m)
E   = 6.37e10;   % Young’s modulus of the panel material (Pa)
rho = 1.2;       % Air density (kg/m^3)
nu  = 0.24;      % Poisson’s ratio of the panel material (-)
c   = 343;       % Speed of sound in air (m/s)
fc = 2079;       % Coincidence frequency (Hz)
f = 50:1:2*fc;   % Frequency range (Hz)
eta=0.01;        % Damping loss factor of glass (-)
% Bending stiffness
B = E*h^3/(12*(1 - nu^2));   % Bending stiffness of the plate (N·m)

STL_normal = zeros(size(f));     % Normal incidence (x = 0)
STL_diffuse = zeros(size(f));    % Diffuse field (angular integration)

for i = 1:length(f)
    w = 2*pi*f(i);  % Angular frequency (rad/s)

    % 1) Diffuse-field to 78 (angular integration)
    tau_diffuse = @(theta) 1 ./ ( (1 + (m*w - B*w^3.*(sin(theta)/c).^4).^2 ...
       + B^2*eta^2*w^6.*(sin(theta)/c).^8).* (cos(theta)./(2*c*rho)).^2 );

    sol = integral(@(theta) tau_diffuse(theta).*sin(2*theta), 0, 78*pi/180);
    STL_diffuse(i) = 10*log10(1/sol);

    % 2) Diffuse-field to 90 (angular integration)
    tau_diffuse2 = @(theta) 1 ./ ( (1 + (m*w - B*w^3.*(sin(theta)/c).^4).^2 ...
       + B^2*eta^2*w^6.*(sin(theta)/c).^8).* (cos(theta)./(2*c*rho)).^2 );
    sol2 = integral(@(theta) tau_diffuse2(theta).*sin(2*theta), 0, 89*pi/180);
    STL_diffuse2(i) = 10*log10(1/sol2);
end

% Plot STL
semilogx(f, STL_diffuse, 'LineWidth', 1.5)      % Diffuse-field 78
hold on
semilogx(f, STL_diffuse2, 'LineWidth', 1.5)     % Diffuse-field 90
xlabel('Frequency (Hz)')
ylabel('STL (dB)')
title('Sound Transmission Loss (Cremer)')
legend('Diffuse field 78', 'Diffuse field 90')
grid on
axis([50 2*fc 5 50])
