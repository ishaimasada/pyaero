%{
Preliminary (Annular) Combustor Design
Sizing Procedure --> from Dr. Ferrer-Vidal
%}

clear; clc;

%%  INPUTS

% Inlet Station
Tt3 = 647;
Tt4 = 1650;
Pt3 = 647.39e3;
Pt4 = 622.97e3;
mdot3 = 20.9;
mdot_fuel = 0.577;
LHV = 43.1e6;
R = 287;
gamma = 1.4;

%  Thermodynamic parameters
Pt_loss = 0.05;
omega_cold = 16;
K_OTDF = 0.05;
K_hot = 1.3;
liner_area_frac = 0.5;
r_tip = 0.34;
Lcomb = 0.45;
FAR_stoich = 1/15;
phi_PRZ = 1.02;
phi_SEC = 0.60;


%%  CALCULATIONS

FAR_overall = mdot_fuel / mdot3;
phi_overall = FAR_overall / FAR_stoich;
TR = Tt4 / Tt3;
omega_hot = K_hot * (TR - 1);
omega_ref = omega_cold + omega_hot;
Aref = sqrt( (R/2) * (mdot3 * sqrt(Tt3) / Pt3)^2 * (omega_ref / Pt_loss) );
dPt_check = omega_ref * (R/2) * (mdot3 * sqrt(Tt3) / (Aref * Pt3))^2;
rho_t3 = Pt3 / (R * Tt3);
Vref   = mdot3 / (rho_t3 * Aref);
q_ref  = 0.5 * rho_t3 * Vref^2;
Mref   = Vref / sqrt(gamma * R * Tt3);
Aliner = liner_area_frac * Aref;

if pi * r_tip^2 <= Aref; error('r_tip too small. Increase r_tip so that pi*r_tip^2 > Aref = %.6f m^2', Aref); end

r_hub = sqrt(r_tip^2 - Aref / pi);
Dl    = (r_tip - r_hub);

mdot_air_PRZ    = mdot_fuel / (phi_PRZ * FAR_stoich);
mdot_air_to_SEC = mdot_fuel / (phi_SEC * FAR_stoich);
mdot_air_SEC    = mdot_air_to_SEC - mdot_air_PRZ;
mdot_air_DIL    = mdot3 - mdot_air_to_SEC;

frac_PRZ = (mdot_air_PRZ / mdot3) * 100;
frac_SEC = (mdot_air_SEC / mdot3) * 100;
frac_DIL = (mdot_air_DIL / mdot3) * 100;
OTDF = 1 - exp( 1 / (-K_OTDF * (Lcomb / Dl) * omega_cold) );
Vol = Aref * Lcomb;
tau_res    = Lcomb / Vref;
tau_res_ms = tau_res * 1000;
Pt3_atm = Pt3 / 101325;
theta_i = (mdot_fuel * LHV) / (Vol * Pt3_atm);
theta_L = mdot3 / (Vol * (Pt3_atm^1.8) * 10^(0.00145 * (Tt3 - 400)));
eta_comb = -5.46974e-10 * theta_L^5 + 3.97923e-8  * theta_L^4  - 8.73718e-6  * theta_L^3 + 3.00007e-4  * theta_L^2  - 4.568246e-3 * theta_L   + 99.7;
Aheff = Aref / sqrt(omega_cold);


%%  RESULTS CHECKS

all_pass = true;

if OTDF > 0.25
    fprintf('FLAG: OTDF = %.4f exceeds 0.25. Increase Lcomb or adjust liner sizing.\n', OTDF);
    all_pass = false;
end

if (theta_i / 1e6) > 60
    fprintf('FLAG: theta_i = %.2f MW/(m^3*atm) exceeds 60 SLS limit. Increase volume.\n', theta_i / 1e6);
    all_pass = false;
end

if theta_L > 5
    fprintf('FLAG: theta_L = %.4f exceeds 5 kg/(s*atm^1.8*m^3) SLS stability limit.\n', theta_L);
    all_pass = false;
end

if tau_res_ms < 3
    fprintf('FLAG: Residence time = %.3f ms is below 3 ms minimum. Increase Lcomb.\n', tau_res_ms);
    all_pass = false;
end

if eta_comb < 95
    fprintf('FLAG: Combustion efficiency = %.2f %% below 95%%. Reduce loading or increase volume.\n', eta_comb);
    all_pass = false;
end

if phi_PRZ < 0.8 || phi_PRZ > 1.3
    fprintf('FLAG: phi_PRZ = %.4f outside stable ignition range 0.8-1.3.\n', phi_PRZ);
    all_pass = false;
end

if all_pass
    fprintf('Design passes limitation checks.\n');
end


%% RESULTS

fprintf('\n============================================\n');
fprintf(' Combustor Preliminary Sizing Results\n');
fprintf('============================================\n');
fprintf('\n--- Station Conditions ---\n');
fprintf('Tt3:                                    %.2f K\n',    Tt3);
fprintf('Tt4:                                    %.2f K\n',    Tt4);
fprintf('Pt3:                                    %.3f kPa\n',  Pt3 / 1e3);
fprintf('Pt4:                                    %.3f kPa\n',  Pt4 / 1e3);
fprintf('mdot3:                              %.3f kg/s\n', mdot3);
fprintf('mdot_fuel:                       %.4f kg/s\n', mdot_fuel);
fprintf('FAR (overall):                   %.5f\n',      FAR_overall);
fprintf('phi (overall):                    %.4f\n',      phi_overall);

fprintf('\n--- Loss Coefficients ---\n');
fprintf('TR = Tt4/Tt3:                           %.4f\n', TR);
fprintf('omega_cold:                            %.2f\n',  omega_cold);
fprintf('omega_hot:                             %.4f\n',  omega_hot);
fprintf('omega_ref:                              %.4f\n',  omega_ref);
fprintf('dPt/Pt (target):                       %.4f\n',  Pt_loss);
fprintf('dPt/Pt (check from Aref):     %.4f\n',  dPt_check);

fprintf('\n--- Reference Quantities ---\n');
fprintf('Aref:                                   %.6f m^2\n',    Aref);
fprintf('rho_t3:                               %.4f kg/m^3\n', rho_t3);
fprintf('Vref:                                    %.4f m/s\n',    Vref);
fprintf('q_ref:                                  %.2f Pa\n',     q_ref);
fprintf('Mref:                                   %.5f\n',        Mref);

fprintf('\n--- Liner Geometry ---\n');
fprintf('Aliner (%d%% of Aref):    %.6f m^2\n', round(liner_area_frac*100), Aliner);
fprintf('r_tip:                                    %.4f m\n',    r_tip);
fprintf('r_hub:                                  %.4f m\n',    r_hub);
fprintf('Dl = 2*(r_tip - r_hub):      %.4f m\n',    Dl);

fprintf('\n--- Zone Air Distribution ---\n');
fprintf('PRZ  (phi = %.2f):  mdot = %.4f kg/s   (%.1f%% total)\n', phi_PRZ, mdot_air_PRZ, frac_PRZ);
fprintf('SEC  (phi = %.2f):  mdot = %.4f kg/s   (%.1f%% total)\n', phi_SEC, mdot_air_SEC, frac_SEC);
fprintf('DIL  (remainder):   mdot = %.4f kg/s   (%.1f%% total)\n', mdot_air_DIL, frac_DIL);

fprintf('\n--- Length, Volume, Residence Time ---\n');
fprintf('Lcomb (input):                     %.4f m\n',  Lcomb);
fprintf('OTDF:                                     %.4f\n',    OTDF);
fprintf('Volume (Aref * Lcomb):     %.6f m^3\n', Vol);
fprintf('Residence time:                   %.4f ms\n', tau_res_ms);

fprintf('\n--- Liner Hole Area ---\n');
fprintf('Aheff (total effective hole area):      %.6f m^2\n', Aheff);

fprintf('\n--- Combustor Loading and Efficiency ---\n');
fprintf('Pt3:                                               %.4f atm\n',                  Pt3_atm);
fprintf('theta_i (intensity):                    %.4f MW/(m^3*atm)\n',        theta_i / 1e6);
fprintf('theta_L (stability loading):      %.6f kg/(s*atm^1.8*m^3)\n', theta_L);
fprintf('Combustion efficiency:            %.3f %%\n',                  eta_comb);