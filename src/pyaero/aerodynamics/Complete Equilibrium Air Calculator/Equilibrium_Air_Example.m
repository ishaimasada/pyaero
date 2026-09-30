clear
clc

%% Inputs
h = 30000; % Altitude in meters
M = 8; % Mach number
beta = 10; % Shock angle in degrees
gamma_inf = 1.4; % Cold Air Specific Heat Ratio
TCALC = 1; % Set to 0 for stagnation point heating or 1 for shock-wave relation

%% Determine Mixture Gamma
gamma_mix = equil_air(h, M, beta, gamma_inf, TCALC);