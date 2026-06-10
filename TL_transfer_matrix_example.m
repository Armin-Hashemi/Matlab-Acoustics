clc;
clear;

f = 100;
c = 343;                % Speed of sound (m/s)
rho = 1.225;            % Density of gas (kg/m^3)
S = 0.01;               % Cross-sectional area (m^2)
L = 0.5;                % Length of the duct segment (m)

omega = 2 * pi * f;     % Angular frequency (rad/s)
k = omega / c;          % Wavenumber (m^-1)
Y = S / (rho * c);      % Characteristic admittance (m^3/(N·s))

T_straight = [cos(k * L), 1i * (rho * c / S) * sin(k * L); ...
              1i * (S / (rho * c)) * sin(k * L), cos(k * L)];
T_exp = eye(2);
T_cont = eye(2);

T = T_exp * T_straight * T_cont;

% Extract matrix elements
T11 = T(1,1);
T12 = T(1,2);
T21 = T(2,1);
T22 = T(2,2);

% For equal inlet/outlet media, Y1 = Yn = Y
Y1 = Y;
Yn = Y;

TL = 20 * log10(abs(0.5 * (T11 + (T12 / Y1) + (Yn * T21) + ((Yn / Y1) * T22))));

fprintf('TL = %.3f dB\n', TL);
