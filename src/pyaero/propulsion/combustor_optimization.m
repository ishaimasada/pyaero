clear; clc; close all;

%%  INPUTS 

% Inlet Station (SLS)
Pt3   = 2002.39e3;    % Pa
Tt3   = 707.3;        % K
Pt4   = 1917.97e3;    % Pa
Tt4   = 1650;         % K
mdot3 = 62.011;       % kg/s
R = 287;   % J/(kg*K)
mdot_fuel = 1.71;     % kg/s
LHV       = 43.1e6;   % J/kg

%  Thermodynamic parameters
dPt_over_Pt = 0.07;
K_cold = 16;
K_hot  = 1.3;
liner_area_frac = 0.6;
K_OTDF = 0.05;

% Acceptable limits
OTDF_min = 0.10;
OTDF_max = 0.25;
theta_i_max_MW = 60;   % MW/(m^3*atm)

% Target combustor length
L_target = 0.4;   % m


%%  CALCULATIONS

TR = Tt4 / Tt3;
wref  = K_cold + K_hot * (TR - 1);
wcold = K_cold;

% Reference area from pressure loss sizing
Aref = sqrt( (R/2) * (mdot3 * sqrt(Tt3) / Pt3)^2 * (wref / dPt_over_Pt) );

% Liner area
Aliner = liner_area_frac * Aref;

% Pressure in atm for combustor loading
Pt3_atm = Pt3 / 101325;

% Choose which area annulus geometry must satisfy
Ageom = Aref;
% Ageom = Aliner;

%%  PARAMETRIC SWEEP RANGES

r_tip_vec = linspace(0.22, 0.28, 140);   % m
Lcomb_vec = linspace(0.35, 0.55, 140);   % m, centered near 0.45 m
[RTIP, LGRID] = meshgrid(r_tip_vec, Lcomb_vec);
RHUB    = nan(size(RTIP));
DL      = nan(size(RTIP));
OTDF    = nan(size(RTIP));
VOL     = nan(size(RTIP));
THETAI  = nan(size(RTIP));   % MW/(m^3*atm)
FEAS    = false(size(RTIP));

%%  LOOP THROUGH GEOMETRY

for i = 1:numel(RTIP)
    r_tip = RTIP(i);
    Lcomb = LGRID(i);

    % Physical feasibility check
    if pi * r_tip^2 <= Ageom
        continue
    end

    % Solve r_hub from annulus area
    r_hub = sqrt(r_tip^2 - Ageom/pi);

    % Annulus height
    Dl = r_tip - r_hub;

    % OTDF from Ferrer-Vidal style correlation
    otdf_val = 1 - exp( 1 / ( -K_OTDF * (Lcomb / Dl) * wcold ) );

    % Notes-based concept volume
    Vol = Aref * Lcomb;

    % Combustor loading
    theta_i = (mdot_fuel * LHV) / (Vol * Pt3_atm) / 1e6;   % MW/(m^3*atm)

    % Store
    RHUB(i)   = r_hub;
    DL(i)     = Dl;
    OTDF(i)   = otdf_val;
    VOL(i)    = Vol;
    THETAI(i) = theta_i;

    % Acceptable region
    if (otdf_val >= OTDF_min) && (otdf_val <= OTDF_max) && (theta_i <= theta_i_max_MW)
        FEAS(i) = true;
    end
end

%% =========================
%  DISPLAY FEASIBLE REGION SUMMARY
%  =========================

n_feas = nnz(FEAS);

fprintf('\n=== PARAMETRIC SWEEP SUMMARY ===\n');
fprintf('Aref used for geometry:                 %.6f m^2\n', Aref);
fprintf('Aliner:                                 %.6f m^2\n', Aliner);
fprintf('Target Lcomb:                           %.6f m\n', L_target);
fprintf('Number of feasible points:              %d\n', n_feas);

if n_feas == 0
    fprintf('No feasible points found in the chosen r_tip and Lcomb ranges.\n');
else
    % Feasible values
    rtip_feas   = RTIP(FEAS);
    Lcomb_feas  = LGRID(FEAS);
    rhub_feas   = RHUB(FEAS);
    Dl_feas     = DL(FEAS);
    OTDF_feas   = OTDF(FEAS);
    thetai_feas = THETAI(FEAS);
    Vol_feas    = VOL(FEAS);

    fprintf('\nFeasible ranges:\n');
    fprintf('r_tip range:                            %.6f to %.6f m\n', min(rtip_feas), max(rtip_feas));
    fprintf('Lcomb range:                            %.6f to %.6f m\n', min(Lcomb_feas), max(Lcomb_feas));
    fprintf('r_hub range:                            %.6f to %.6f m\n', min(rhub_feas), max(rhub_feas));
    fprintf('Dl range:                               %.6f to %.6f m\n', min(Dl_feas), max(Dl_feas));
    fprintf('OTDF range:                             %.4f to %.4f\n', min(OTDF_feas), max(OTDF_feas));
    fprintf('theta_i range:                          %.4f to %.4f MW/(m^3*atm)\n', min(thetai_feas), max(thetai_feas));

    % Pick best design, strongly favor Lcomb near 0.45 m
    OTDF_target_mid = 0.5 * (OTDF_min + OTDF_max);
    L_err = abs(Lcomb_feas - L_target);

    score = L_err + 0.2 * abs(OTDF_feas - OTDF_target_mid) + 0.002 * thetai_feas;
    [~, idx_best] = min(score);

    fprintf('\n=== BEST FEASIBLE DESIGN NEAR Lcomb = %.3f m ===\n', L_target);
    fprintf('r_tip:                                  %.6f m\n', rtip_feas(idx_best));
    fprintf('r_hub:                                  %.6f m\n', rhub_feas(idx_best));
    fprintf('Dl:                                     %.6f m\n', Dl_feas(idx_best));
    fprintf('Lcomb:                                  %.6f m\n', Lcomb_feas(idx_best));
    fprintf('OTDF:                                   %.4f\n', OTDF_feas(idx_best));
    fprintf('theta_i:                                %.4f MW/(m^3*atm)\n', thetai_feas(idx_best));
    fprintf('Volume:                                 %.6f m^3\n', Vol_feas(idx_best));
end

%% =========================
%  PLOTS
%  =========================

% OTDF contour
figure;
contourf(RTIP, LGRID, OTDF, 20, 'LineColor', 'none');
hold on;
contour(RTIP, LGRID, OTDF, [OTDF_min OTDF_min], 'k', 'LineWidth', 1.5);
contour(RTIP, LGRID, OTDF, [OTDF_max OTDF_max], 'k', 'LineWidth', 1.5);
grid on;
xlabel('r_{tip} (m)');
ylabel('L_{comb} (m)');
title('OTDF over r_{tip} and L_{comb}');
colorbar;

% theta_i contour
figure;
contourf(RTIP, LGRID, THETAI, 20, 'LineColor', 'none');
hold on;
contour(RTIP, LGRID, THETAI, [theta_i_max_MW theta_i_max_MW], 'k', 'LineWidth', 1.5);
grid on;
xlabel('r_{tip} (m)');
ylabel('L_{comb} (m)');
title('\theta_i over r_{tip} and L_{comb}, MW/(m^3*atm)');
colorbar;

% Feasible region map
figure;
imagesc(r_tip_vec, Lcomb_vec, FEAS);
set(gca, 'YDir', 'normal');
grid on;
xlabel('r_{tip} (m)');
ylabel('L_{comb} (m)');
title('Feasible Region, 1 = acceptable, 0 = unacceptable');
colorbar;

% Feasible design scatter
figure;
plot(RTIP(FEAS), LGRID(FEAS), '.');
grid on;
xlabel('r_{tip} (m)');
ylabel('L_{comb} (m)');
title('Feasible r_{tip} and L_{comb} combinations');

% Optional, show only designs close to target length
if n_feas > 0
    figure;
    idx_close = FEAS & abs(LGRID - L_target) <= 0.01;   % within +/- 0.01 m
    plot(RTIP(idx_close), LGRID(idx_close), '.');
    grid on;
    xlabel('r_{tip} (m)');
    ylabel('L_{comb} (m)');
    title('Feasible designs within +/- 0.01 m of target L_{comb}');
end