% Ian Doughty
% AE 595CC - Feb. 4, 2026
% This script plots the molar concentration of each specie across a 200K to
% 15000K temperature range.

%% Inputs
% Standard atmosphere
h = 30000; % Altitude in meters
[rho, P_atm, T_atm] = standard_atm(h); % Calculates the standard atmosphere conditions

% Temperature range
T_vec = linspace(200, 15000, 400);   % Increase points for smoother curves

% Storage
S_store = zeros(length(T_vec), 11);   % mol/kg
X_store = zeros(length(T_vec), 11);   % mole fraction

for k = 1:length(T_vec)
    T = T_vec(k);

    % Thermo equilibrium constants at this T and rho
    [K1, K2, K3, K4, K5, K6, K7] = thermo_data(T, rho);
    
   % Determine model viability temperature ranges based on density
    rho_log = log10(rho); % Log of density
    Tmin = 200; % Minimum required temperature for analysis
    T0 = 1600; % Model 1 viable below this temperature in K
    T1 = 1000 * (3.5 + 10*((rho_log+6)/7)^2.5); % Model 2 viable between T0 and T1
    T2 = 1000 * (5.4 + 3.2*(rho_log+6)/6 + 5.55*(abs(rho_log+6)/6)^3.35); % Model 3 viable between T1 and T2
    Tmax = 15000; % Maximum allowed temperature for analysis
    
    if T >= Tmin && T < T0
        % Model 1
        S = model_1(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar);
    elseif T >= T0 && T < T1
        % Model 2
        S = model_2(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar);
    elseif T >= T1 && T < T2
        % Model 3
        S = model_3(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar);
    elseif T >= T2 && T <= Tmax
        % Model 4
        S = model_4(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar);
    else
        error('Temperature must be between 200K and 15000K')
    end

    % Store sigmas
    S_store(k,:) = S;
end

%% Plot Molar Concentration

% Plot settings
Tscale = 1000;  % Plot x as T/1000

% Manual label positioning
xLabelPos = nan(1,11);
yLabelPos = nan(1,11);

% 1 O2
xLabelPos(1) = 0.8;  yLabelPos(1) = -0.95;
% 2 N2
xLabelPos(2) = 0.8;  yLabelPos(2) = -0.35;
% 3 O
xLabelPos(3) = 13.5; yLabelPos(3) = -1.1;
% 4 NO
xLabelPos(4) = 5.4;  yLabelPos(4) = -2.64;
% 5 N
xLabelPos(5) = 14.1; yLabelPos(5) = -0.25;
% 6 NO+
xLabelPos(6) = 7.9;  yLabelPos(6) = -4;
% 7 e-
xLabelPos(7) = 12.9; yLabelPos(7) = -0.5;
% 8 N+
xLabelPos(8) = 12; yLabelPos(8) = -1.25;
% 9 O+
xLabelPos(9) = 14; yLabelPos(9) = -1.65;
% 10 Ar
xLabelPos(10)= 0.8;  yLabelPos(10)= -1.8;
% 11 Ar+
xLabelPos(11)= 13; yLabelPos(11)= -3.3;

% Create label names
names = {'O_2','N_2','O','NO','N','NO^+','e^-','N^+','O^+','Ar','Ar^+'};

% Color-code species into families
c_air  = [0 0.4470 0.7410]; % O2, N2, Ar
c_neut = [0.4660 0.6740 0.1880]; % O, N, NO
c_ion  = [0.8500 0.3250 0.0980]; % ions + e-

idx_air  = [1 2 10];
idx_neut = [3 4 5];
idx_ion  = [6 7 8 9 11];

colors = zeros(11,3);
colors(idx_air,:)  = repmat(c_air,  numel(idx_air),  1);
colors(idx_neut,:) = repmat(c_neut, numel(idx_neut), 1);
colors(idx_ion,:)  = repmat(c_ion,  numel(idx_ion),  1);

% Compute y-axis
Tplot = T_vec(:)/Tscale;  % 10^3 K
Yplot = nan(size(S_store));

for k = 1:length(Tplot)
    sig = sum(S_store(k,:));
    if sig > 0
        Yplot(k,:) = log10(S_store(k,:)/sig);
    end
end

% Plot molar concentration curves vs. temp
figure; hold on;

for j = 1:11
    if ismember(j, idx_air)
        lw = 2.6;
        ls = '-';
    elseif ismember(j, idx_neut)
        lw = 2.1;
        ls = '-';
    else
        lw = 1.5;
        ls = '-'; 
    end

    plot(Tplot, Yplot(:,j), 'LineWidth', lw, 'LineStyle', ls, 'Color', colors(j,:));
end

xlabel('Temperature (10^3 K)');
ylabel('log_{10}(\sigma_i / \sigma)');
xlim([0 15]);
ylim([-5 0]);
grid off; box on;

% Plot leaders for labels
for j = 1:11
    if isnan(xLabelPos(j)) || isnan(yLabelPos(j))
        continue
    end

    % Find anchor point on the curve at the same x as the label (or closest)
    yj = Yplot(:,j);
    good = isfinite(yj);

    if nnz(good) < 2
        continue
    end

    % Interpolate y at xLabelPos(j) for anchor
    xA = xLabelPos(j);
    yA = interp1(Tplot(good), yj(good), xA, 'linear', NaN);

    % If interpolation failed (outside valid range), use nearest finite point
    if ~isfinite(yA)
        [~, idxNear] = min(abs(Tplot(good) - xA));
        tGood = Tplot(good);
        yGood = yj(good);
        xA = tGood(idxNear);
        yA = yGood(idxNear);
    end

    % Label position
    xL = xLabelPos(j);
    yL = yLabelPos(j);

    % Leader line
    plot([xA xL], [yA yL], '-', 'Color', colors(j,:), 'LineWidth', 0.9);

    % Label text with white background for readability
    text(xL, yL, names{j}, ...
        'Color', colors(j,:), ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

hold off;
title(sprintf('Equilibrium Air Composition (\\rho = %.1e kg/m^3)', rho));