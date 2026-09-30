% Ian Doughty
% AE 595CC - Feb. 4, 2026
% This function calculates the equilibrium air conditions over a range of
% temperatures from 200K to 15000K. It also determines the gamma value for
% a given temperature and density.
% NOTE: If using stagnation point heating, T must be found using the NASA
% CEA calculator and matching total enthalpies at the stagnation pressure
% calculated and the correct temperature.

function [gamma_mix] = equil_air(h, M, beta, gamma_inf, TCALC)
    R_air = 287; % Specific Gas Constant(J/kg*K)
    Cp_air = 1005; % Cold Air Specific Heat (J/kg*K)
    [rho, P_atm, T_atm] = standard_atm(h); % Calculates the standard atmosphere conditions
    
    %% Temperature Calculation
    
    if TCALC == 0
        H_o = Cp_air * T_atm; % Static Enthalpy (J/kg)
        a = sqrt(gamma_inf*R_air*T_atm); % Speed of Sound (m/s)
        V = M * a; % Flow Velocity (m/s)
        H = H_o + 0.5*V^2; % Total Freestream Enthalpy (J/kg)
        P_0 = P_atm * ((((gamma_inf+1)^2)*(M^2))/((4*gamma_inf*(M^2))- ...
            (2*(gamma_inf-1))))^(gamma_inf/(gamma_inf-1)) * ((1-gamma_inf+(2*gamma_inf*(M^2))) ...
            /(gamma_inf+1)); % Stagnation Pressure (Pa)
        
        % Use NASA CEA code to get stagnation temperature
        T = 2767.18; % Temperature experienced by vehicle (K)
    
    elseif TCALC == 1
        % Oblique Shock Relations
        M_n = M*sin(deg2rad(beta)); % Normal Mach number
        T = T_atm * ((((2*gamma_inf*M_n^2) - (gamma_inf-1)) * ((gamma_inf-1)* ...
            (M_n^2) + 2)) / (((gamma_inf+1)^2)*(M_n^2)));
    end
    
    %% Initial Sigmas
    % Specie concentrations in air
    y_O2 = 0.2096;
    y_N2 = 0.7811;
    y_Ar = 0.0093;
    
    % Molecular weights
    MW_O2 = 31.998;
    MW_N2 = 28.0134;
    MW_Ar = 39.948;
    
    % Mixture weight
    M_mix = y_O2*MW_O2 + y_N2*MW_N2 + y_Ar*MW_Ar; % kg/kmol
    
    % Molecular concentrations (kmol of atoms per kg)
    S_O2 = 2*y_O2/M_mix; 
    S_N2 = 2*y_N2/M_mix;
    S_Ar = y_Ar/M_mix;
    
    %% Model Implementation
    % Function to calculate the equilibrium coeffs for each reaction
    [K1, K2, K3, K4, K5, K6, K7, Cp_vec] = thermo_data(T, rho); 
    
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
    
    %% Gamma Calculation
    gamma_mix = mixture_gamma(S, Cp_vec);
    
    %% Model 1 Function
    function S = model_1(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar)
           
        % Initialize Sigma Array
        S = zeros(1, 11);
    
        % Set Loop Parameter
        max_iters = 500;
        tol = 1e-6;
        j = 0;
        res1 = Inf;
        res2 = Inf;
        
        % Initial Guesses
        S1p = S_O2/2;
        S2p = S_N2/2;
        S_10 = S_Ar;
    
        while (res1 > tol || res2 > tol) && j < max_iters
            j = j+1;
    
            % Compute Minor Species
            S_3 = sqrt(K1*S1p);
            S_5 = sqrt(K3*S2p);
            
            % Safe Ke computation for very low T
            if (K1 <= 0) || (K2 <= 0) || (K3 <= 0)
                Ke = 0;   
            else
                logKe = 0.5*(log(K2) - log(K1) - log(K3));
                logKe = max(min(logKe, 700), -700); % Clip to avoid overflow/underflow in exp
                Ke = exp(logKe);
            end
           
            S_4 = Ke*S_5*S_3;
        
            % Update Major Species
            S_1 = (S_O2 - S_3 - S_4)/2;
            S_2 = (S_N2 - S_4 - S_5)/2;
    
            % Convergence Check
            res1 = abs(S_1 - S1p);
            res2 = abs(S_2 - S2p);
    
            % Update Previous Guess
            S1p = S_1;
            S2p = S_2;
        end
        
        % Compute Trace Species (robust at very low T)
        % If any ionization equilibrium constants underflow to ~0, ions are negligible
        Kmin = realmin;  
        if (K4 < Kmin) || (K5 < Kmin) || (K6 < Kmin) || (K7 < Kmin)
            S_6  = 0;
            S_7  = 0;
            S_8  = 0;
            S_9  = 0;
            S_11 = 0;
            S_10 = S_Ar;
        else
            arg = (S_4/K4) + (S_5/K5) + (S_3/K6) + (S_10/K7);
        
            % If arg is ~0 (or negative from roundoff), ions are negligible
            if arg <= 0
                S_6  = 0;
                S_7  = 0;
                S_8  = 0;
                S_9  = 0;
                S_11 = 0;
                S_10 = S_Ar;
            else
                S_7  = sqrt(arg);
                S_6  = S_4/(K4*S_7);
                S_8  = S_5/(K5*S_7);
                S_9  = S_3/(K6*S_7);
                S_11 = S_10/(K7*S_7);
            end
        end
    
        % Generate outputs
        S = [S_1, S_2, S_3, S_4, S_5, S_6, S_7, S_8, S_9, S_10, S_11];
    end
    
    %% Model 2 Function
    function S = model_2(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar)
        % Initialize Sigma Array
        S = zeros(1, 11);
    
        % Set Loop Parameters
        max_iters = 500;
        tol = 1e-6;
        j = 0;
        res3 = Inf;
    
        % Constant used in NASA model
        Ke = sqrt(K2/(K1*K3));
    
        % Initial guesses
        S_6  = 0;   % NO+
        S_8  = 0;   % N+
        S_9  = 0;   % O+
        S_10 = S_Ar; % Ar
        S_11 = 0;   % Ar+
    
        % Newton initialization. B4 and B3 are constant in the Newton polynomial
        B4 = (8 - 2*K2) / K1;
        B3 = 8 - K2 - 2*K3*Ke;
    
        % Create constant for simplicity
        C1 = K1*K3*Ke;
    
        % Initial guess for S3
        BB = K1*(1.0 + sqrt((K2*S_N2)/(2.0*K1)));
        S3 = (-BB + sqrt(BB*BB + 8.0*S_O2*K1)) / 4.0;
    
        while (res3 > tol) && (j < max_iters)
            j = j + 1;
    
            % Adjust totals by current ion amounts
            SO1 = S_O2 - S_6 - S_9;
            SN1 = S_N2 - S_6 - S_8;
    
            % Quartic coefficients for S3
            B2 = -8.0*SO1 + 2.0*K1 + (SO1 - SN1)*K2 - C1;
            B1 = SO1*(-4.0*K1 + C1);
            B0 = 2.0*SO1*SO1*K1;
    
            % Quartic and derivative
            F  = B0 + S3*(B1 + S3*(B2 + S3*(B3 + S3*B4)));
            FP = B1 + S3*(2.0*B2 + S3*(3.0*B3 + S3*(4.0*B4)));
    
            % Newton step
            dS3 = -F/FP;
            S3  = S3 + dS3;
    
            % Keep S3 physical
            if S3 < 0
                S3 = 0;
            end
    
            % Compute major species
            denom = Ke*S3;
            if denom <= 0 
                denom = eps; % Keeps S3 non-zero
            end
            S_5 = (SO1 - S3 - 2.0*(S3*S3)/K1) / denom;
    
            S_4 = Ke*S_5*S3;
            S_2 = (SN1 - S_4 - S_5)/2.0;
            S_1 = (SO1 - S_4 - S3)/2.0;
    
            % Compute the trace species
            S_3 = S3;
            S_7 = sqrt((S_4/K4) + (S_5/K5) + (S_3/K6) + (S_10/K7));
            S_6 = S_4/(K4*S_7);
            S_8 = S_5/(K5*S_7);
            S_9 = S_3/(K6*S_7);
            S_11 = S_10/(K7*S_7);
            S_10 = S_Ar - S_11;
    
            % Convergence check
            res3 = abs(dS3);
        end
    
        % Generate outputs
        S = [S_1, S_2, S_3, S_4, S_5, S_6, S_7, S_8, S_9, S_10, S_11];
    end
    
    %% Model 3 Function
    function S = model_3(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar)
    
        % Initialize Sigma Array
        S = zeros(1, 11);
    
        % Set loop parameters
        max_iters = 500;
        tol = 1e-6;
        j = 0;
        res3 = Inf;
    
        % Constants
        Ke = sqrt(K2/(K1*K3));
    
        % Initialize (NASA)
        S_4  = 0;      % NO
        S_6  = 0;      % NO+
        S_8  = 0;      % N+
        S_10 = S_Ar;   % Ar
        S_11 = 0;      % Ar+
    
        % Initial guesses for quartic
        B1 = K3/2.0;
        B2 = K1/2.0;
    
        % Initial guess for S5 (atomic N)
        S_5 = (-B1 + sqrt(B1*B1 + 4.0*B1*S_N2)) / 2.0;
    
        % Initial guess for S3P (atomic O)
        S3p = (-B2 + sqrt(B2*B2 + 4.0*B2*S_O2)) / 2.0;
    
        while (res3 > tol) && (j < max_iters)
            j = j + 1;
    
            % Adjust N total for current NO, NO+, N+
            SN1 = S_N2 - S_4 - S_6 - S_8;
    
            % Newton update for S5 
            F  = -B1*SN1 + S_5*(B1 + S_5);
            FP = B1 + 2.0*S_5;
    
            dS5 = -F/FP;
            S_5 = S_5 + dS5;
    
            if S_5 < 0
                S_5 = 0;
            end
    
            % N2 from mass balance
            S_2 = (SN1 - S_5)/2.0;
    
            % Compute remaining species
            S_1 = (S3p*S3p)/K1;        % O2
            S_4 = Ke*S3p*S_5;          % NO
    
            % Compute ions/electrons
            S_7 = sqrt((S3p/K6) + (S_4/K4) + (S_5/K5) + (S_10/K7));
            S_6 = S_4/(K4*S_7);
            S_8 = S_5/(K5*S_7);
            S_9 = S3p/(K6*S_7);
            S_11 = S_10/(K7*S_7);
    
            % Argon
            S_10 = S_Ar - S_11;
    
            % Atomic oxygen from oxygen mass balance
            S_3 = S_O2 - 2.0*S_1 - S_4 - S_6 - S_9;
    
            % Convergence check based on S3p
            dS3 = S3p - S_3;
            res3 = abs(dS3);
    
            % Update guess
            S3p = S_3;
        end
    
        % Generate outputs
        S = [S_1, S_2, S_3, S_4, S_5, S_6, S_7, S_8, S_9, S_10, S_11];
    end
    
    %% Model 4 Function (NASA "CUBIC" = Model IV)
    function S = model_4(K1, K2, K3, K4, K5, K6, K7, S_O2, S_N2, S_Ar)
    
        % Initialize Sigma Array
        S = zeros(1, 11);
    
        % Set loop parameters
        max_iters = 500;
        tol = 1e-6;
        j = 0;
        res7 = Inf;
    
        % Constants
        Ke = sqrt(K2/(K1*K3));
    
        % Initialize loop
        S_1 = 0.0;   % O2
        S_2 = 0.0;   % N2
        S_4 = 0.0;   % NO
        S_6 = 0.0;   % NO+
        S_10 = S_Ar; % Ar
        S_11 = 0.0;  % Ar+
    
        % Initial guess for electron concentration
        S_7 = (sqrt(1.0 + 4.0*S_N2*K5) - 1.0) / (2.0*K5);
    
        % Initial O+ and Ar+ guesses
        S_9  = S_O2 / (1.0 + K6*S_7);
        S_11 = S_Ar / (1.0 + K7*S_7);
    
        % Update electron guess
        S_7 = S_7 + S_9 + S_11;
    
        % Create constants to simplify
        B3 = K5*K6;
        C2 = K5 + K6;
    
        while (res7 > tol) && (j < max_iters)
            j = j + 1;
    
            SO1 = S_O2 - 2.0*S_1 - S_4 - S_6;
            SN1 = S_N2 - 2.0*S_2 - S_4 - S_6;
    
            S6P = S_6 + S_11;
    
            % Cubic coefficients
            B2 = C2 - B3*S6P;
            B1 = 1.0 - SN1*K6 - SO1*K5 - S6P*C2;
            B0 = -SN1 - SO1 - S6P;
    
            F  = B0 + S_7*(B1 + S_7*(B2 + S_7*B3));
            FP = B1 + S_7*(2.0*B2 + S_7*(3.0*B3));
    
            % Newton step for electron
            dS7 = -F/FP;
            S_7 = S_7 + dS7;
    
            % Keep electron sigma non-negative and real
            if S_7 < 0
                S_7 = 0;
            end
    
            % Compute major species
            S_3  = SO1*S_7*K6 / (1.0 + K6*S_7);      % O
            S_5  = SN1*S_7*K5 / (1.0 + K5*S_7);      % N
            S_11 = S_Ar / (1.0 + K7*S_7);            % Ar+
            S_10 = S_Ar - S_11;                      % Ar
            S_9  = S_3 / (K6*S_7);                   % O+
            S_8  = S_5 / (K5*S_7);                   % N+
    
            % Compute minor species
            S_1 = (S_3*S_3)/K1;                      % O2
            S_2 = (S_5*S_5)/K3;                      % N2
            S_4 = Ke*S_3*S_5;                        % NO
            S_6 = S_4 / (K4*S_7);                    % NO+
    
            % Convergence check
            res7 = abs(dS7);
        end
    
        % Output vector in the NASA species order:
        S = [S_1, S_2, S_3, S_4, S_5, S_6, S_7, S_8, S_9, S_10, S_11];
    end
    
    %% Gamma Calculation Function
    function gamma_mix = mixture_gamma(S, Cp_vec)
    % S: 1x11 in kmol/kg (your sigma outputs)
    % Cp_molar: 1x11 in J/(mol*K)
    
        % MW [kg/kmol] in your species order:
        MW = [31.998, 28.0134, 15.999, 30.006, 14.0067, 30.006, ...
              0.00054858, 14.0067, 15.999, 39.948, 39.948];
    
        Ru_kmol = 8314.462618;    % J/(kmol*K)
    
        % Convert Cp to mass-specific: J/(kg*K)
        % Cp_mass = Cp_molar [J/mol-K] * (1000 mol/kmol) / (MW kg/kmol)
        Cp_mass = Cp_vec .* (1000 ./ MW);
    
        % Species gas constants (mass-specific): J/(kg*K)
        R_mass = Ru_kmol ./ MW;
    
        % Cv per species (mass-specific): J/(kg*K)
        Cv_mass = Cp_mass - R_mass;
    
        % Mass fractions from S (kmol/kg)
        mass_i = S .* MW;          % (kmol/kg)*(kg/kmol) = dimensionless mass per kg mixture
        Y = mass_i ./ sum(mass_i); % mass fractions
    
        % Mixture Cp, Cv (mass basis)
        Cp_mix = sum(Y .* Cp_mass);
        Cv_mix = sum(Y .* Cv_mass);
    
        gamma_mix = Cp_mix / Cv_mix;
    end
end