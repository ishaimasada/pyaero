% Ian Doughty
% AE 595CC - Feb. 4, 2026
% This function calculates the equilibrium constants in terms of mole
% numbers for each reaction.

function [K1, K2, K3, K4, K5, K6, K7, Cp_vec] = thermo_data(T, rho)

    % Inputs
    P_o = 101325; % Reference pressure (Pa)
    R = 8314.46; % Universal gas constant (kJ/mol*K)
    dN13 = 1; % Increase in mole number for reactions 1 and 3
    dN4_7 = -1; % Increase in mole number for reactions 4,5,6, and 7
    dN2 = 0; % Increase in mole number for reaction 2

    % Read-in thermo data from excel file
    aO2 = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C3:I7');
    aN2 = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C9:I13');
    aO = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C15:I19');
    aNO = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C21:I25');
    aN = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C27:I31');
    aNO_pos = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C33:I37');
    aE_neg = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C39:I43');
    aN_pos = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C45:I49');
    aO_pos = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C51:I55');
    aAr = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C57:I61');
    aAr_pos = readmatrix('Prabhu and Erickson Thermo Database.xlsx', 'Range','C63:I67');

    % Calculate Enthalpy and Entropy for each specie at each temperature range
    if T >= 200 && T < 800
        H_O2 = aO2(1,1) + aO2(1,2)*T/2 + aO2(1,3)*(T^2)/3 + aO2(1,4)*(T^3)/4 + aO2(1,5)*(T^4)/5 + aO2(1,6)/T;
        H_N2 = aN2(1,1) + aN2(1,2)*T/2 + aN2(1,3)*(T^2)/3 + aN2(1,4)*(T^3)/4 + aN2(1,5)*(T^4)/5 + aN2(1,6)/T;
        H_O  = aO(1,1)  + aO(1,2)*T/2  + aO(1,3)*(T^2)/3  + aO(1,4)*(T^3)/4  + aO(1,5)*(T^4)/5  + aO(1,6)/T;
        H_NO = aNO(1,1) + aNO(1,2)*T/2 + aNO(1,3)*(T^2)/3 + aNO(1,4)*(T^3)/4 + aNO(1,5)*(T^4)/5 + aNO(1,6)/T;
        H_N  = aN(1,1)  + aN(1,2)*T/2  + aN(1,3)*(T^2)/3  + aN(1,4)*(T^3)/4  + aN(1,5)*(T^4)/5  + aN(1,6)/T;
        H_NO_pos = aNO_pos(1,1) + aNO_pos(1,2)*T/2 + aNO_pos(1,3)*(T^2)/3 + aNO_pos(1,4)*(T^3)/4 + aNO_pos(1,5)*(T^4)/5 + aNO_pos(1,6)/T;
        H_E_neg  = aE_neg(1,1)  + aE_neg(1,2)*T/2  + aE_neg(1,3)*(T^2)/3  + aE_neg(1,4)*(T^3)/4  + aE_neg(1,5)*(T^4)/5  + aE_neg(1,6)/T;
        H_N_pos  = aN_pos(1,1)  + aN_pos(1,2)*T/2  + aN_pos(1,3)*(T^2)/3  + aN_pos(1,4)*(T^3)/4  + aN_pos(1,5)*(T^4)/5  + aN_pos(1,6)/T;
        H_O_pos  = aO_pos(1,1)  + aO_pos(1,2)*T/2  + aO_pos(1,3)*(T^2)/3  + aO_pos(1,4)*(T^3)/4  + aO_pos(1,5)*(T^4)/5  + aO_pos(1,6)/T;
        H_Ar = aAr(1,1) + aAr(1,2)*T/2 + aAr(1,3)*(T^2)/3 + aAr(1,4)*(T^3)/4 + aAr(1,5)*(T^4)/5 + aAr(1,6)/T;
        H_Ar_pos = aAr_pos(1,1) + aAr_pos(1,2)*T/2 + aAr_pos(1,3)*(T^2)/3 + aAr_pos(1,4)*(T^3)/4 + aAr_pos(1,5)*(T^4)/5 + aAr_pos(1,6)/T;

        S_O2 = aO2(1,1)*log(T) + aO2(1,2)*T + aO2(1,3)*(T^2)/2 + aO2(1,4)*(T^3)/3 + aO2(1,5)*(T^4)/4 + aO2(1,7);
        S_N2 = aN2(1,1)*log(T) + aN2(1,2)*T + aN2(1,3)*(T^2)/2 + aN2(1,4)*(T^3)/3 + aN2(1,5)*(T^4)/4 + aN2(1,7);
        S_O  = aO(1,1)*log(T)  + aO(1,2)*T  + aO(1,3)*(T^2)/2  + aO(1,4)*(T^3)/3  + aO(1,5)*(T^4)/4  + aO(1,7);
        S_NO = aNO(1,1)*log(T) + aNO(1,2)*T + aNO(1,3)*(T^2)/2 + aNO(1,4)*(T^3)/3 + aNO(1,5)*(T^4)/4 + aNO(1,7);
        S_N  = aN(1,1)*log(T)  + aN(1,2)*T  + aN(1,3)*(T^2)/2  + aN(1,4)*(T^3)/3  + aN(1,5)*(T^4)/4  + aN(1,7);
        S_NO_pos = aNO_pos(1,1)*log(T) + aNO_pos(1,2)*T + aNO_pos(1,3)*(T^2)/2 + aNO_pos(1,4)*(T^3)/3 + aNO_pos(1,5)*(T^4)/4 + aNO_pos(1,7);
        S_E_neg  = aE_neg(1,1)*log(T)  + aE_neg(1,2)*T  + aE_neg(1,3)*(T^2)/2  + aE_neg(1,4)*(T^3)/3  + aE_neg(1,5)*(T^4)/4  + aE_neg(1,7);
        S_N_pos  = aN_pos(1,1)*log(T)  + aN_pos(1,2)*T  + aN_pos(1,3)*(T^2)/2  + aN_pos(1,4)*(T^3)/3  + aN_pos(1,5)*(T^4)/4  + aN_pos(1,7);
        S_O_pos  = aO_pos(1,1)*log(T)  + aO_pos(1,2)*T  + aO_pos(1,3)*(T^2)/2  + aO_pos(1,4)*(T^3)/3  + aO_pos(1,5)*(T^4)/4  + aO_pos(1,7);
        S_Ar = aAr(1,1)*log(T) + aAr(1,2)*T + aAr(1,3)*(T^2)/2 + aAr(1,4)*(T^3)/3 + aAr(1,5)*(T^4)/4 + aAr(1,7);
        S_Ar_pos = aAr_pos(1,1)*log(T) + aAr_pos(1,2)*T + aAr_pos(1,3)*(T^2)/2 + aAr_pos(1,4)*(T^3)/3 + aAr_pos(1,5)*(T^4)/4 + aAr_pos(1,7);

        Cp_O2     = aO2(1,1)     + aO2(1,2)*T     + aO2(1,3)*T^2     + aO2(1,4)*T^3     + aO2(1,5)*T^4;
        Cp_N2     = aN2(1,1)     + aN2(1,2)*T     + aN2(1,3)*T^2     + aN2(1,4)*T^3     + aN2(1,5)*T^4;
        Cp_O      = aO(1,1)      + aO(1,2)*T      + aO(1,3)*T^2      + aO(1,4)*T^3      + aO(1,5)*T^4;
        Cp_NO     = aNO(1,1)     + aNO(1,2)*T     + aNO(1,3)*T^2     + aNO(1,4)*T^3     + aNO(1,5)*T^4;
        Cp_N      = aN(1,1)      + aN(1,2)*T      + aN(1,3)*T^2      + aN(1,4)*T^3      + aN(1,5)*T^4;
        Cp_NO_pos = aNO_pos(1,1) + aNO_pos(1,2)*T + aNO_pos(1,3)*T^2 + aNO_pos(1,4)*T^3 + aNO_pos(1,5)*T^4;
        Cp_E_neg  = aE_neg(1,1)  + aE_neg(1,2)*T  + aE_neg(1,3)*T^2  + aE_neg(1,4)*T^3  + aE_neg(1,5)*T^4;
        Cp_N_pos  = aN_pos(1,1)  + aN_pos(1,2)*T  + aN_pos(1,3)*T^2  + aN_pos(1,4)*T^3  + aN_pos(1,5)*T^4;
        Cp_O_pos  = aO_pos(1,1)  + aO_pos(1,2)*T  + aO_pos(1,3)*T^2  + aO_pos(1,4)*T^3  + aO_pos(1,5)*T^4;
        Cp_Ar     = aAr(1,1)     + aAr(1,2)*T     + aAr(1,3)*T^2     + aAr(1,4)*T^3     + aAr(1,5)*T^4;
        Cp_Ar_pos = aAr_pos(1,1) + aAr_pos(1,2)*T + aAr_pos(1,3)*T^2 + aAr_pos(1,4)*T^3 + aAr_pos(1,5)*T^4;

    elseif T >= 800 && T < 3000
        H_O2 = aO2(2,1) + aO2(2,2)*T/2 + aO2(2,3)*(T^2)/3 + aO2(2,4)*(T^3)/4 + aO2(2,5)*(T^4)/5 + aO2(2,6)/T;
        H_N2 = aN2(2,1) + aN2(2,2)*T/2 + aN2(2,3)*(T^2)/3 + aN2(2,4)*(T^3)/4 + aN2(2,5)*(T^4)/5 + aN2(2,6)/T;
        H_O  = aO(2,1)  + aO(2,2)*T/2  + aO(2,3)*(T^2)/3  + aO(2,4)*(T^3)/4  + aO(2,5)*(T^4)/5  + aO(2,6)/T;
        H_NO = aNO(2,1) + aNO(2,2)*T/2 + aNO(2,3)*(T^2)/3 + aNO(2,4)*(T^3)/4 + aNO(2,5)*(T^4)/5 + aNO(2,6)/T;
        H_N  = aN(2,1)  + aN(2,2)*T/2  + aN(2,3)*(T^2)/3  + aN(2,4)*(T^3)/4  + aN(2,5)*(T^4)/5  + aN(2,6)/T;
        H_NO_pos = aNO_pos(2,1) + aNO_pos(2,2)*T/2 + aNO_pos(2,3)*(T^2)/3 + aNO_pos(2,4)*(T^3)/4 + aNO_pos(2,5)*(T^4)/5 + aNO_pos(2,6)/T;
        H_E_neg  = aE_neg(2,1)  + aE_neg(2,2)*T/2  + aE_neg(2,3)*(T^2)/3  + aE_neg(2,4)*(T^3)/4  + aE_neg(2,5)*(T^4)/5  + aE_neg(2,6)/T;
        H_N_pos  = aN_pos(2,1)  + aN_pos(2,2)*T/2  + aN_pos(2,3)*(T^2)/3  + aN_pos(2,4)*(T^3)/4  + aN_pos(2,5)*(T^4)/5  + aN_pos(2,6)/T;
        H_O_pos  = aO_pos(2,1)  + aO_pos(2,2)*T/2  + aO_pos(2,3)*(T^2)/3  + aO_pos(2,4)*(T^3)/4  + aO_pos(2,5)*(T^4)/5  + aO_pos(2,6)/T;
        H_Ar = aAr(2,1) + aAr(2,2)*T/2 + aAr(2,3)*(T^2)/3 + aAr(2,4)*(T^3)/4 + aAr(2,5)*(T^4)/5 + aAr(2,6)/T;
        H_Ar_pos = aAr_pos(2,1) + aAr_pos(2,2)*T/2 + aAr_pos(2,3)*(T^2)/3 + aAr_pos(2,4)*(T^3)/4 + aAr_pos(2,5)*(T^4)/5 + aAr_pos(2,6)/T;

        S_O2 = aO2(2,1)*log(T) + aO2(2,2)*T + aO2(2,3)*(T^2)/2 + aO2(2,4)*(T^3)/3 + aO2(2,5)*(T^4)/4 + aO2(2,7);
        S_N2 = aN2(2,1)*log(T) + aN2(2,2)*T + aN2(2,3)*(T^2)/2 + aN2(2,4)*(T^3)/3 + aN2(2,5)*(T^4)/4 + aN2(2,7);
        S_O  = aO(2,1)*log(T)  + aO(2,2)*T  + aO(2,3)*(T^2)/2  + aO(2,4)*(T^3)/3  + aO(2,5)*(T^4)/4  + aO(2,7);
        S_NO = aNO(2,1)*log(T) + aNO(2,2)*T + aNO(2,3)*(T^2)/2 + aNO(2,4)*(T^3)/3 + aNO(2,5)*(T^4)/4 + aNO(2,7);
        S_N  = aN(2,1)*log(T)  + aN(2,2)*T  + aN(2,3)*(T^2)/2  + aN(2,4)*(T^3)/3  + aN(2,5)*(T^4)/4  + aN(2,7);
        S_NO_pos = aNO_pos(2,1)*log(T) + aNO_pos(2,2)*T + aNO_pos(2,3)*(T^2)/2 + aNO_pos(2,4)*(T^3)/3 + aNO_pos(2,5)*(T^4)/4 + aNO_pos(2,7);
        S_E_neg  = aE_neg(2,1)*log(T)  + aE_neg(2,2)*T  + aE_neg(2,3)*(T^2)/2  + aE_neg(2,4)*(T^3)/3  + aE_neg(2,5)*(T^4)/4  + aE_neg(2,7);
        S_N_pos  = aN_pos(2,1)*log(T)  + aN_pos(2,2)*T  + aN_pos(2,3)*(T^2)/2  + aN_pos(2,4)*(T^3)/3  + aN_pos(2,5)*(T^4)/4  + aN_pos(2,7);
        S_O_pos  = aO_pos(2,1)*log(T)  + aO_pos(2,2)*T  + aO_pos(2,3)*(T^2)/2  + aO_pos(2,4)*(T^3)/3  + aO_pos(2,5)*(T^4)/4  + aO_pos(2,7);
        S_Ar = aAr(2,1)*log(T) + aAr(2,2)*T + aAr(2,3)*(T^2)/2 + aAr(2,4)*(T^3)/3 + aAr(2,5)*(T^4)/4 + aAr(2,7);
        S_Ar_pos = aAr_pos(2,1)*log(T) + aAr_pos(2,2)*T + aAr_pos(2,3)*(T^2)/2 + aAr_pos(2,4)*(T^3)/3 + aAr_pos(2,5)*(T^4)/4 + aAr_pos(2,7);

        Cp_O2     = aO2(2,1)     + aO2(2,2)*T     + aO2(2,3)*T^2     + aO2(2,4)*T^3     + aO2(2,5)*T^4;
        Cp_N2     = aN2(2,1)     + aN2(2,2)*T     + aN2(2,3)*T^2     + aN2(2,4)*T^3     + aN2(2,5)*T^4;
        Cp_O      = aO(2,1)      + aO(2,2)*T      + aO(2,3)*T^2      + aO(2,4)*T^3      + aO(2,5)*T^4;
        Cp_NO     = aNO(2,1)     + aNO(2,2)*T     + aNO(2,3)*T^2     + aNO(2,4)*T^3     + aNO(2,5)*T^4;
        Cp_N      = aN(2,1)      + aN(2,2)*T      + aN(2,3)*T^2      + aN(2,4)*T^3      + aN(2,5)*T^4;
        Cp_NO_pos = aNO_pos(2,1) + aNO_pos(2,2)*T + aNO_pos(2,3)*T^2 + aNO_pos(2,4)*T^3 + aNO_pos(2,5)*T^4;
        Cp_E_neg  = aE_neg(2,1)  + aE_neg(2,2)*T  + aE_neg(2,3)*T^2  + aE_neg(2,4)*T^3  + aE_neg(2,5)*T^4;
        Cp_N_pos  = aN_pos(2,1)  + aN_pos(2,2)*T  + aN_pos(2,3)*T^2  + aN_pos(2,4)*T^3  + aN_pos(2,5)*T^4;
        Cp_O_pos  = aO_pos(2,1)  + aO_pos(2,2)*T  + aO_pos(2,3)*T^2  + aO_pos(2,4)*T^3  + aO_pos(2,5)*T^4;
        Cp_Ar     = aAr(2,1)     + aAr(2,2)*T     + aAr(2,3)*T^2     + aAr(2,4)*T^3     + aAr(2,5)*T^4;
        Cp_Ar_pos = aAr_pos(2,1) + aAr_pos(2,2)*T + aAr_pos(2,3)*T^2 + aAr_pos(2,4)*T^3 + aAr_pos(2,5)*T^4;

    elseif T >= 3000 && T < 6000
        H_O2 = aO2(3,1) + aO2(3,2)*T/2 + aO2(3,3)*(T^2)/3 + aO2(3,4)*(T^3)/4 + aO2(3,5)*(T^4)/5 + aO2(3,6)/T;
        H_N2 = aN2(3,1) + aN2(3,2)*T/2 + aN2(3,3)*(T^2)/3 + aN2(3,4)*(T^3)/4 + aN2(3,5)*(T^4)/5 + aN2(3,6)/T;
        H_O  = aO(3,1)  + aO(3,2)*T/2  + aO(3,3)*(T^2)/3  + aO(3,4)*(T^3)/4  + aO(3,5)*(T^4)/5  + aO(3,6)/T;
        H_NO = aNO(3,1) + aNO(3,2)*T/2 + aNO(3,3)*(T^2)/3 + aNO(3,4)*(T^3)/4 + aNO(3,5)*(T^4)/5 + aNO(3,6)/T;
        H_N  = aN(3,1)  + aN(3,2)*T/2  + aN(3,3)*(T^2)/3  + aN(3,4)*(T^3)/4  + aN(3,5)*(T^4)/5  + aN(3,6)/T;
        H_NO_pos = aNO_pos(3,1) + aNO_pos(3,2)*T/2 + aNO_pos(3,3)*(T^2)/3 + aNO_pos(3,4)*(T^3)/4 + aNO_pos(3,5)*(T^4)/5 + aNO_pos(3,6)/T;
        H_E_neg  = aE_neg(3,1)  + aE_neg(3,2)*T/2  + aE_neg(3,3)*(T^2)/3  + aE_neg(3,4)*(T^3)/4  + aE_neg(3,5)*(T^4)/5  + aE_neg(3,6)/T;
        H_N_pos  = aN_pos(3,1)  + aN_pos(3,2)*T/2  + aN_pos(3,3)*(T^2)/3  + aN_pos(3,4)*(T^3)/4  + aN_pos(3,5)*(T^4)/5  + aN_pos(3,6)/T;
        H_O_pos  = aO_pos(3,1)  + aO_pos(3,2)*T/2  + aO_pos(3,3)*(T^2)/3  + aO_pos(3,4)*(T^3)/4  + aO_pos(3,5)*(T^4)/5  + aO_pos(3,6)/T;
        H_Ar = aAr(3,1) + aAr(3,2)*T/2 + aAr(3,3)*(T^2)/3 + aAr(3,4)*(T^3)/4 + aAr(3,5)*(T^4)/5 + aAr(3,6)/T;
        H_Ar_pos = aAr_pos(3,1) + aAr_pos(3,2)*T/2 + aAr_pos(3,3)*(T^2)/3 + aAr_pos(3,4)*(T^3)/4 + aAr_pos(3,5)*(T^4)/5 + aAr_pos(3,6)/T;

        S_O2 = aO2(3,1)*log(T) + aO2(3,2)*T + aO2(3,3)*(T^2)/2 + aO2(3,4)*(T^3)/3 + aO2(3,5)*(T^4)/4 + aO2(3,7);
        S_N2 = aN2(3,1)*log(T) + aN2(3,2)*T + aN2(3,3)*(T^2)/2 + aN2(3,4)*(T^3)/3 + aN2(3,5)*(T^4)/4 + aN2(3,7);
        S_O  = aO(3,1)*log(T)  + aO(3,2)*T  + aO(3,3)*(T^2)/2  + aO(3,4)*(T^3)/3  + aO(3,5)*(T^4)/4  + aO(3,7);
        S_NO = aNO(3,1)*log(T) + aNO(3,2)*T + aNO(3,3)*(T^2)/2 + aNO(3,4)*(T^3)/3 + aNO(3,5)*(T^4)/4 + aNO(3,7);
        S_N  = aN(3,1)*log(T)  + aN(3,2)*T  + aN(3,3)*(T^2)/2  + aN(3,4)*(T^3)/3  + aN(3,5)*(T^4)/4  + aN(3,7);
        S_NO_pos = aNO_pos(3,1)*log(T) + aNO_pos(3,2)*T + aNO_pos(3,3)*(T^2)/2 + aNO_pos(3,4)*(T^3)/3 + aNO_pos(3,5)*(T^4)/4 + aNO_pos(3,7);
        S_E_neg  = aE_neg(3,1)*log(T)  + aE_neg(3,2)*T  + aE_neg(3,3)*(T^2)/2  + aE_neg(3,4)*(T^3)/3  + aE_neg(3,5)*(T^4)/4  + aE_neg(3,7);
        S_N_pos  = aN_pos(3,1)*log(T)  + aN_pos(3,2)*T  + aN_pos(3,3)*(T^2)/2  + aN_pos(3,4)*(T^3)/3  + aN_pos(3,5)*(T^4)/4  + aN_pos(3,7);
        S_O_pos  = aO_pos(3,1)*log(T)  + aO_pos(3,2)*T  + aO_pos(3,3)*(T^2)/2  + aO_pos(3,4)*(T^3)/3  + aO_pos(3,5)*(T^4)/4  + aO_pos(3,7);
        S_Ar = aAr(3,1)*log(T) + aAr(3,2)*T + aAr(3,3)*(T^2)/2 + aAr(3,4)*(T^3)/3 + aAr(3,5)*(T^4)/4 + aAr(3,7);
        S_Ar_pos = aAr_pos(3,1)*log(T) + aAr_pos(3,2)*T + aAr_pos(3,3)*(T^2)/2 + aAr_pos(3,4)*(T^3)/3 + aAr_pos(3,5)*(T^4)/4 + aAr_pos(3,7);

        Cp_O2     = aO2(3,1)     + aO2(3,2)*T     + aO2(3,3)*T^2     + aO2(3,4)*T^3     + aO2(3,5)*T^4;
        Cp_N2     = aN2(3,1)     + aN2(3,2)*T     + aN2(3,3)*T^2     + aN2(3,4)*T^3     + aN2(3,5)*T^4;
        Cp_O      = aO(3,1)      + aO(3,2)*T      + aO(3,3)*T^2      + aO(3,4)*T^3      + aO(3,5)*T^4;
        Cp_NO     = aNO(3,1)     + aNO(3,2)*T     + aNO(3,3)*T^2     + aNO(3,4)*T^3     + aNO(3,5)*T^4;
        Cp_N      = aN(3,1)      + aN(3,2)*T      + aN(3,3)*T^2      + aN(3,4)*T^3      + aN(3,5)*T^4;
        Cp_NO_pos = aNO_pos(3,1) + aNO_pos(3,2)*T + aNO_pos(3,3)*T^2 + aNO_pos(3,4)*T^3 + aNO_pos(3,5)*T^4;
        Cp_E_neg  = aE_neg(3,1)  + aE_neg(3,2)*T  + aE_neg(3,3)*T^2  + aE_neg(3,4)*T^3  + aE_neg(3,5)*T^4;
        Cp_N_pos  = aN_pos(3,1)  + aN_pos(3,2)*T  + aN_pos(3,3)*T^2  + aN_pos(3,4)*T^3  + aN_pos(3,5)*T^4;
        Cp_O_pos  = aO_pos(3,1)  + aO_pos(3,2)*T  + aO_pos(3,3)*T^2  + aO_pos(3,4)*T^3  + aO_pos(3,5)*T^4;
        Cp_Ar     = aAr(3,1)     + aAr(3,2)*T     + aAr(3,3)*T^2     + aAr(3,4)*T^3     + aAr(3,5)*T^4;
        Cp_Ar_pos = aAr_pos(3,1) + aAr_pos(3,2)*T + aAr_pos(3,3)*T^2 + aAr_pos(3,4)*T^3 + aAr_pos(3,5)*T^4;

    
    elseif T >= 6000 && T < 10000
        H_O2 = aO2(4,1) + aO2(4,2)*T/2 + aO2(4,3)*(T^2)/3 + aO2(4,4)*(T^3)/4 + aO2(4,5)*(T^4)/5 + aO2(4,6)/T;
        H_N2 = aN2(4,1) + aN2(4,2)*T/2 + aN2(4,3)*(T^2)/3 + aN2(4,4)*(T^3)/4 + aN2(4,5)*(T^4)/5 + aN2(4,6)/T;
        H_O  = aO(4,1)  + aO(4,2)*T/2  + aO(4,3)*(T^2)/3  + aO(4,4)*(T^3)/4  + aO(4,5)*(T^4)/5  + aO(4,6)/T;
        H_NO = aNO(4,1) + aNO(4,2)*T/2 + aNO(4,3)*(T^2)/3 + aNO(4,4)*(T^3)/4 + aNO(4,5)*(T^4)/5 + aNO(4,6)/T;
        H_N  = aN(4,1)  + aN(4,2)*T/2  + aN(4,3)*(T^2)/3  + aN(4,4)*(T^3)/4  + aN(4,5)*(T^4)/5  + aN(4,6)/T;
        H_NO_pos = aNO_pos(4,1) + aNO_pos(4,2)*T/2 + aNO_pos(4,3)*(T^2)/3 + aNO_pos(4,4)*(T^3)/4 + aNO_pos(4,5)*(T^4)/5 + aNO_pos(4,6)/T;
        H_E_neg  = aE_neg(4,1)  + aE_neg(4,2)*T/2  + aE_neg(4,3)*(T^2)/3  + aE_neg(4,4)*(T^3)/4  + aE_neg(4,5)*(T^4)/5  + aE_neg(4,6)/T;
        H_N_pos  = aN_pos(4,1)  + aN_pos(4,2)*T/2  + aN_pos(4,3)*(T^2)/3  + aN_pos(4,4)*(T^3)/4  + aN_pos(4,5)*(T^4)/5  + aN_pos(4,6)/T;
        H_O_pos  = aO_pos(4,1)  + aO_pos(4,2)*T/2  + aO_pos(4,3)*(T^2)/3  + aO_pos(4,4)*(T^3)/4  + aO_pos(4,5)*(T^4)/5  + aO_pos(4,6)/T;
        H_Ar = aAr(4,1) + aAr(4,2)*T/2 + aAr(4,3)*(T^2)/3 + aAr(4,4)*(T^3)/4 + aAr(4,5)*(T^4)/5 + aAr(4,6)/T;
        H_Ar_pos = aAr_pos(4,1) + aAr_pos(4,2)*T/2 + aAr_pos(4,3)*(T^2)/3 + aAr_pos(4,4)*(T^3)/4 + aAr_pos(4,5)*(T^4)/5 + aAr_pos(4,6)/T;

        S_O2 = aO2(4,1)*log(T) + aO2(4,2)*T + aO2(4,3)*(T^2)/2 + aO2(4,4)*(T^3)/3 + aO2(4,5)*(T^4)/4 + aO2(4,7);
        S_N2 = aN2(4,1)*log(T) + aN2(4,2)*T + aN2(4,3)*(T^2)/2 + aN2(4,4)*(T^3)/3 + aN2(4,5)*(T^4)/4 + aN2(4,7);
        S_O  = aO(4,1)*log(T)  + aO(4,2)*T  + aO(4,3)*(T^2)/2  + aO(4,4)*(T^3)/3  + aO(4,5)*(T^4)/4  + aO(4,7);
        S_NO = aNO(4,1)*log(T) + aNO(4,2)*T + aNO(4,3)*(T^2)/2 + aNO(4,4)*(T^3)/3 + aNO(4,5)*(T^4)/4 + aNO(4,7);
        S_N  = aN(4,1)*log(T)  + aN(4,2)*T  + aN(4,3)*(T^2)/2  + aN(4,4)*(T^3)/3  + aN(4,5)*(T^4)/4  + aN(4,7);
        S_NO_pos = aNO_pos(4,1)*log(T) + aNO_pos(4,2)*T + aNO_pos(4,3)*(T^2)/2 + aNO_pos(4,4)*(T^3)/3 + aNO_pos(4,5)*(T^4)/4 + aNO_pos(4,7);
        S_E_neg  = aE_neg(4,1)*log(T)  + aE_neg(4,2)*T  + aE_neg(4,3)*(T^2)/2  + aE_neg(4,4)*(T^3)/3  + aE_neg(4,5)*(T^4)/4  + aE_neg(4,7);
        S_N_pos  = aN_pos(4,1)*log(T)  + aN_pos(4,2)*T  + aN_pos(4,3)*(T^2)/2  + aN_pos(4,4)*(T^3)/3  + aN_pos(4,5)*(T^4)/4  + aN_pos(4,7);
        S_O_pos  = aO_pos(4,1)*log(T)  + aO_pos(4,2)*T  + aO_pos(4,3)*(T^2)/2  + aO_pos(4,4)*(T^3)/3  + aO_pos(4,5)*(T^4)/4  + aO_pos(4,7);
        S_Ar = aAr(4,1)*log(T) + aAr(4,2)*T + aAr(4,3)*(T^2)/2 + aAr(4,4)*(T^3)/3 + aAr(4,5)*(T^4)/4 + aAr(4,7);
        S_Ar_pos = aAr_pos(4,1)*log(T) + aAr_pos(4,2)*T + aAr_pos(4,3)*(T^2)/2 + aAr_pos(4,4)*(T^3)/3 + aAr_pos(4,5)*(T^4)/4 + aAr_pos(4,7);

        Cp_O2     = aO2(4,1)     + aO2(4,2)*T     + aO2(4,3)*T^2     + aO2(4,4)*T^3     + aO2(4,5)*T^4;
        Cp_N2     = aN2(4,1)     + aN2(4,2)*T     + aN2(4,3)*T^2     + aN2(4,4)*T^3     + aN2(4,5)*T^4;
        Cp_O      = aO(4,1)      + aO(4,2)*T      + aO(4,3)*T^2      + aO(4,4)*T^3      + aO(4,5)*T^4;
        Cp_NO     = aNO(4,1)     + aNO(4,2)*T     + aNO(4,3)*T^2     + aNO(4,4)*T^3     + aNO(4,5)*T^4;
        Cp_N      = aN(4,1)      + aN(4,2)*T      + aN(4,3)*T^2      + aN(4,4)*T^3      + aN(4,5)*T^4;
        Cp_NO_pos = aNO_pos(4,1) + aNO_pos(4,2)*T + aNO_pos(4,3)*T^2 + aNO_pos(4,4)*T^3 + aNO_pos(4,5)*T^4;
        Cp_E_neg  = aE_neg(4,1)  + aE_neg(4,2)*T  + aE_neg(4,3)*T^2  + aE_neg(4,4)*T^3  + aE_neg(4,5)*T^4;
        Cp_N_pos  = aN_pos(4,1)  + aN_pos(4,2)*T  + aN_pos(4,3)*T^2  + aN_pos(4,4)*T^3  + aN_pos(4,5)*T^4;
        Cp_O_pos  = aO_pos(4,1)  + aO_pos(4,2)*T  + aO_pos(4,3)*T^2  + aO_pos(4,4)*T^3  + aO_pos(4,5)*T^4;
        Cp_Ar     = aAr(4,1)     + aAr(4,2)*T     + aAr(4,3)*T^2     + aAr(4,4)*T^3     + aAr(4,5)*T^4;
        Cp_Ar_pos = aAr_pos(4,1) + aAr_pos(4,2)*T + aAr_pos(4,3)*T^2 + aAr_pos(4,4)*T^3 + aAr_pos(4,5)*T^4;

    
    elseif T >= 10000 && T <= 15000
        H_O2 = aO2(5,1) + aO2(5,2)*T/2 + aO2(5,3)*(T^2)/3 + aO2(5,4)*(T^3)/4 + aO2(5,5)*(T^4)/5 + aO2(5,6)/T;
        H_N2 = aN2(5,1) + aN2(5,2)*T/2 + aN2(5,3)*(T^2)/3 + aN2(5,4)*(T^3)/4 + aN2(5,5)*(T^4)/5 + aN2(5,6)/T;
        H_O  = aO(5,1)  + aO(5,2)*T/2  + aO(5,3)*(T^2)/3  + aO(5,4)*(T^3)/4  + aO(5,5)*(T^4)/5  + aO(5,6)/T;
        H_NO = aNO(5,1) + aNO(5,2)*T/2 + aNO(5,3)*(T^2)/3 + aNO(5,4)*(T^3)/4 + aNO(5,5)*(T^4)/5 + aNO(5,6)/T;
        H_N  = aN(5,1)  + aN(5,2)*T/2  + aN(5,3)*(T^2)/3  + aN(5,4)*(T^3)/4  + aN(5,5)*(T^4)/5  + aN(5,6)/T;
        H_NO_pos = aNO_pos(5,1) + aNO_pos(5,2)*T/2 + aNO_pos(5,3)*(T^2)/3 + aNO_pos(5,4)*(T^3)/4 + aNO_pos(5,5)*(T^4)/5 + aNO_pos(5,6)/T;
        H_E_neg  = aE_neg(5,1)  + aE_neg(5,2)*T/2  + aE_neg(5,3)*(T^2)/3  + aE_neg(5,4)*(T^3)/4  + aE_neg(5,5)*(T^4)/5  + aE_neg(5,6)/T;
        H_N_pos  = aN_pos(5,1)  + aN_pos(5,2)*T/2  + aN_pos(5,3)*(T^2)/3  + aN_pos(5,4)*(T^3)/4  + aN_pos(5,5)*(T^4)/5  + aN_pos(5,6)/T;
        H_O_pos  = aO_pos(5,1)  + aO_pos(5,2)*T/2  + aO_pos(5,3)*(T^2)/3  + aO_pos(5,4)*(T^3)/4  + aO_pos(5,5)*(T^4)/5  + aO_pos(5,6)/T;
        H_Ar = aAr(5,1) + aAr(5,2)*T/2 + aAr(5,3)*(T^2)/3 + aAr(5,4)*(T^3)/4 + aAr(5,5)*(T^4)/5 + aAr(5,6)/T;
        H_Ar_pos = aAr_pos(5,1) + aAr_pos(5,2)*T/2 + aAr_pos(5,3)*(T^2)/3 + aAr_pos(5,4)*(T^3)/4 + aAr_pos(5,5)*(T^4)/5 + aAr_pos(5,6)/T;

        S_O2 = aO2(5,1)*log(T) + aO2(5,2)*T + aO2(5,3)*(T^2)/2 + aO2(5,4)*(T^3)/3 + aO2(5,5)*(T^4)/4 + aO2(5,7);
        S_N2 = aN2(5,1)*log(T) + aN2(5,2)*T + aN2(5,3)*(T^2)/2 + aN2(5,4)*(T^3)/3 + aN2(5,5)*(T^4)/4 + aN2(5,7);
        S_O  = aO(5,1)*log(T)  + aO(5,2)*T  + aO(5,3)*(T^2)/2  + aO(5,4)*(T^3)/3  + aO(5,5)*(T^4)/4  + aO(5,7);
        S_NO = aNO(5,1)*log(T) + aNO(5,2)*T + aNO(5,3)*(T^2)/2 + aNO(5,4)*(T^3)/3 + aNO(5,5)*(T^4)/4 + aNO(5,7);
        S_N  = aN(5,1)*log(T)  + aN(5,2)*T  + aN(5,3)*(T^2)/2  + aN(5,4)*(T^3)/3  + aN(5,5)*(T^4)/4  + aN(5,7);
        S_NO_pos = aNO_pos(5,1)*log(T) + aNO_pos(5,2)*T + aNO_pos(5,3)*(T^2)/2 + aNO_pos(5,4)*(T^3)/3 + aNO_pos(5,5)*(T^4)/4 + aNO_pos(5,7);
        S_E_neg  = aE_neg(5,1)*log(T)  + aE_neg(5,2)*T  + aE_neg(5,3)*(T^2)/2  + aE_neg(5,4)*(T^3)/3  + aE_neg(5,5)*(T^4)/4  + aE_neg(5,7);
        S_N_pos  = aN_pos(5,1)*log(T)  + aN_pos(5,2)*T  + aN_pos(5,3)*(T^2)/2  + aN_pos(5,4)*(T^3)/3  + aN_pos(5,5)*(T^4)/4  + aN_pos(5,7);
        S_O_pos  = aO_pos(5,1)*log(T)  + aO_pos(5,2)*T  + aO_pos(5,3)*(T^2)/2  + aO_pos(5,4)*(T^3)/3  + aO_pos(5,5)*(T^4)/4  + aO_pos(5,7);
        S_Ar = aAr(5,1)*log(T) + aAr(5,2)*T + aAr(5,3)*(T^2)/2 + aAr(5,4)*(T^3)/3 + aAr(5,5)*(T^4)/4 + aAr(5,7);
        S_Ar_pos = aAr_pos(5,1)*log(T) + aAr_pos(5,2)*T + aAr_pos(5,3)*(T^2)/2 + aAr_pos(5,4)*(T^3)/3 + aAr_pos(5,5)*(T^4)/4 + aAr_pos(5,7);

        Cp_O2     = aO2(5,1)     + aO2(5,2)*T     + aO2(5,3)*T^2     + aO2(5,4)*T^3     + aO2(5,5)*T^4;
        Cp_N2     = aN2(5,1)     + aN2(5,2)*T     + aN2(5,3)*T^2     + aN2(5,4)*T^3     + aN2(5,5)*T^4;
        Cp_O      = aO(5,1)      + aO(5,2)*T      + aO(5,3)*T^2      + aO(5,4)*T^3      + aO(5,5)*T^4;
        Cp_NO     = aNO(5,1)     + aNO(5,2)*T     + aNO(5,3)*T^2     + aNO(5,4)*T^3     + aNO(5,5)*T^4;
        Cp_N      = aN(5,1)      + aN(5,2)*T      + aN(5,3)*T^2      + aN(5,4)*T^3      + aN(5,5)*T^4;
        Cp_NO_pos = aNO_pos(5,1) + aNO_pos(5,2)*T + aNO_pos(5,3)*T^2 + aNO_pos(5,4)*T^3 + aNO_pos(5,5)*T^4;
        Cp_E_neg  = aE_neg(5,1)  + aE_neg(5,2)*T  + aE_neg(5,3)*T^2  + aE_neg(5,4)*T^3  + aE_neg(5,5)*T^4;
        Cp_N_pos  = aN_pos(5,1)  + aN_pos(5,2)*T  + aN_pos(5,3)*T^2  + aN_pos(5,4)*T^3  + aN_pos(5,5)*T^4;
        Cp_O_pos  = aO_pos(5,1)  + aO_pos(5,2)*T  + aO_pos(5,3)*T^2  + aO_pos(5,4)*T^3  + aO_pos(5,5)*T^4;
        Cp_Ar     = aAr(5,1)     + aAr(5,2)*T     + aAr(5,3)*T^2     + aAr(5,4)*T^3     + aAr(5,5)*T^4;
        Cp_Ar_pos = aAr_pos(5,1) + aAr_pos(5,2)*T + aAr_pos(5,3)*T^2 + aAr_pos(5,4)*T^3 + aAr_pos(5,5)*T^4;

    else
        error('Temperature must be between 200K and 15000K')
    end
    
    % Calculate Gibb's Free Energy
    G_O2 = H_O2 - S_O2;
    G_N2 = H_N2 - S_N2;
    G_O  = H_O  - S_O;
    G_NO = H_NO - S_NO;
    G_N  = H_N  - S_N;
    G_NO_pos = H_NO_pos - S_NO_pos;
    G_E_neg  = H_E_neg  - S_E_neg;
    G_N_pos  = H_N_pos  - S_N_pos;
    G_O_pos  = H_O_pos  - S_O_pos;
    G_Ar = H_Ar - S_Ar;
    G_Ar_pos = H_Ar_pos - S_Ar_pos;

    % Calculate the change in Gibb's free energy for each reaction
    dG1 = 2*G_O - G_O2;
    dG2 = 2*G_NO - G_N2 - G_O2;
    dG3 = 2*G_N - G_N2;
    dG4 = G_NO - G_NO_pos - G_E_neg;
    dG5 = G_N - G_N_pos - G_E_neg;
    dG6 = G_O - G_O_pos - G_E_neg;
    dG7 = G_Ar - G_Ar_pos - G_E_neg;

    % Calculate the equilibrium constants in terms of partial pressures
    Kp1 = exp(-dG1);
    Kp2 = exp(-dG2);
    Kp3 = exp(-dG3);
    Kp4 = exp(-dG4);
    Kp5 = exp(-dG5);
    Kp6 = exp(-dG6);
    Kp7 = exp(-dG7);

    % Calculate the equilibrium constants in terms of mole numbers
    K1 = Kp1 * (P_o/(rho*R*T))^dN13;
    K2 = Kp2 * (P_o/(rho*R*T))^dN2;
    K3 = Kp3 * (P_o/(rho*R*T))^dN13;
    K4 = Kp4 * (P_o/(rho*R*T))^dN4_7;
    K5 = Kp5 * (P_o/(rho*R*T))^dN4_7;
    K6 = Kp6 * (P_o/(rho*R*T))^dN4_7;
    K7 = Kp7 * (P_o/(rho*R*T))^dN4_7;

    % Extract Cp Values
    Cp_vec = (R/1000).* [Cp_O2, Cp_N2, Cp_O, Cp_NO, Cp_N, Cp_NO_pos, Cp_E_neg, Cp_N_pos, Cp_O_pos, Cp_Ar, Cp_Ar_pos];
end
    