% Ian Doughty
% AE 595CC - Feb. 4, 2026
% This function calculates the density 

function [rho, P_atm, T_atm] = standard_atm(h)
    h_range = 0:50000; % Altitude in meters
    R = 287; % Gas Constant (J/kg*K)
    T = zeros(length(h_range), 1); % Initialize temperature array
    P = zeros(length(h_range), 1); % Initialize pressure array
    for i = 1:length(h_range)
        if i < 11001
            T(i) = 15.04 - 0.00649*i; % Temperature (C)
            P(i) = 101.29 * ((T(i)+273.1)/288.08)^5.256; % Pressure (kPa)
        elseif i >= 11001 && i <= 25001
            T(i) = -56.46; % Temperature (C)
            P(i) = 22.65 * exp(1.73 - 0.000157*i); % Pressure (kPa)
        elseif i > 25001
            T(i) = -131.21 + 0.00299 * i; % Temperature (C)
            P(i) = 2.488 * ((T(i)+273.1)/216.6)^-11.388; % Pressure (kPa)
        end
    end
    
    % Convert degrees celcius to kelvin
    T = T + 273.15;

    % Convert kPa to Pa
    P = P .* 1000;
    
    % Calculate density at altitude
    rho = P./(R.*T);
    
    % Set function outputs
    rho = rho(h+1);
    T_atm = T(h+1);
    P_atm = P(h+1);
end