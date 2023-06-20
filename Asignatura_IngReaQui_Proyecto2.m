

% A: Metanol (CH3OH)
% B: Oxígeno (O2)
% C: Formaldehído (CH2O)
% D: Agua (H2O)
% E: Monóxido de Carbono (CO)
% F: Nitrógeno (N2)
% G: Fluido de Control

% VALORES CONSTANTES / INICIALES:

% Pesos moleculares [kg / mol]:
     PM_A = 32.04 * 1000;
     PM_B = 32.00 * 1000;
     PM_C = 30.03 * 1000;
     PM_D = 18.00 * 1000;
     PM_E = 28.01 * 1000;
     PM_F = 28.01 * 1000;

% Parámetros de Ingreso:
     % Fluido de Reacción
     TR_in = 573; % [K]
     PR_in = 2.5; % [atm]
     ConstR = 0.08205746 * 1000; % [atm m3 / mol K]
     FR_in = 100; % [mol / s]
     
     y_Ain = 0.2; 
     y_Cin = 0.001;
     y_Din = 0.002;
     y_Ein = 0.008;
     y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21;
     y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;
     
     FR_Ain = FR_in * y_Ain;
     FR_Bin = FR_in * y_Bin;
     FR_Cin = FR_in * y_Cin;
     FR_Din = FR_in * y_Din;
     FR_Ein = FR_in * y_Ein;
     FR_Fin = FR_in * y_Fin;

     % Fluido de Control
     Q_G = 
     Ro_G =
     U = 200;
     T_Gin = 273; % [K]
     T_Gout = 350; % [K]

     % Reactor:
     di = 2.0574 / 100; % [m]
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 10; % [m]

% RESOLUCIÓN DE SISTEMA DE ECUACIONES:
Inter = linspace(0, Lmax, 5000);
CI = []; % Condiciones Iniciales
Ensayo = ode45(@(x, T, FR_A, FR_B)Derivadas(), Inter, CI);

% PLANTEAMIENTO DE SOLUCIÓN DE ECUACIONES:
function Ec = Derivadas()

% Calores Específicos [KJ /mol K]: 
%    Coeficientes obtenidos de "Principios Elementales de los Procesos Químicos - Felder" - Tabla B.2 
     Cp_A = @(T) (+42.93 * 10^-3) + (+8.3010 * 10^-5)*(T-273.15) + (-1.8700 * 10^-8) * (T-273.15)^2 + (-8.030 * 10^-12) * (T-273.15)^3;
     Cp_B = @(T) (+29.10 * 10^-3) + (+1.1580 * 10^-5)*(T-273.15) + (-0.6076 * 10^-8) * (T-273.15)^2 + (+1.311 * 10^-12) * (T-273.15)^3;
     Cp_C = @(T) (+34.28 * 10^-3) + (+4.2680 * 10^-5)*(T-273.15) + (+0.0000 * 10^-8) * (T-273.15)^2 + (-8.694 * 10^-12) * (T-273.15)^3;
     Cp_D = @(T) (+33.46 * 10^-3) + (+0.6880 * 10^-5)*(T-273.15) + (+0.7604 * 10^-8) * (T-273.15)^2 + (-3.593 * 10^-12) * (T-273.15)^3;
     Cp_E = @(T) (+28.95 * 10^-3) + (+0.4110 * 10^-5)*(T-273.15) + (+0.3548 * 10^-8) * (T-273.15)^2 + (-2.220 * 10^-12) * (T-273.15)^3;
     Cp_F = @(T) (+29.00 * 10^-3) + (+0.2199 * 10^-5)*(T-273.15) + (+0.5723 * 10^-8) * (T-273.15)^2 + (-2.871 * 10^-12) * (T-273.15)^3;
     Cp_G = @(T) 1; % ( * 10^-3) + ( * 10^-5)*(T-273.15) + ( * 10^-8) * (T-273.15)^2 + ( * 10^-12) * (T-273.15)^3;

% Viscosidades [Pa * s]:
    Mu_A = @(T) 0.00003771 * T + (-0.002267);
    Mu_B = @(T) 0.00004980 * T + (+0.005360);
    Mu_C = @(T) 0.00003771 * T + (-0.002267);
    Mu_D = @(T)	0.00003710 * T + (-0.001140);
    Mu_E = @(T) 0.00004160 * T + (+0.005100);
    Mu_F = @(T) 0.00004160 * T + (+0.005100);

% Cinéticas de Reacción: 
     % Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539;
     Eact_k1 = 62661.3851;
     k1 = A_k1 * exp(-Eact_k1 / (ConstR * T));
     A_K1m = 62661.3851;
     Eact_K1m = 133441.5273;
     k1 = A_K1m * exp(-Eact_K1m / (ConstR * T));
     A_k2 = 0.017971;
     Eact_k2 = 32820.4093;
     k2 = A_k2 * exp(-Eact_k2 / (ConstR * T));
     A_K2w = 1.1538;
     Eact_K2w = -16817.5167;
     K2w = A_K2w * exp(-Eact_K2w / (ConstR * T));

     r1 = CR_A * CR_B ^ 0.5 * k1;
     r2 = CR_C * CR_B ^ 0.5 * k2;

     n1 = 1;
     n2 = 1;

% Relaciones Estequiométricas:
     % Variables: FR_A, Fr_B:
     E1 = @(FR_A) FR_Ain - FR_A;
     E2 = @(FR_A, FR_B) 2 * FR_Bin - E1(FR_A) - 2 * FR_B;
     FR_C = @(FR_A) E1(FR_A);
     FR_D = @(FR_A, FR_B) E1(FR_A) + E2(FR_A, E1(FR_B));
     FR_E = @(FR_A, FR_B) E2(FR_A, E1(FR_B));
     FR_F = FR_Fin;

% Conversión de Flujos a Concentraciones:
     % C = (F/Q) = (P/RT)
     CR_A = FR_A / (Q_in * (FR_A / FR_Ain) * (T / TR_in) * ((P * y_A) / (PR_in * y_Ain)));
     CR_B = FR_B / (Q_in * (FR_B / FR_Bin) * (T / TR_in) * ((P * y_B) / (PR_in * y_Bin)));
     CR_C = FR_C / (Q_in * (FR_C / FR_Cin) * (T / TR_in) * ((P * y_C) / (PR_in * y_Cin)));
     CR_D = FR_D / (Q_in * (FR_D / FR_Din) * (T / TR_in) * ((P * y_D) / (PR_in * y_Din)));
     CR_E = FR_E / (Q_in * (FR_E / FR_Ein) * (T / TR_in) * ((P * y_E) / (PR_in * y_Ein)));
     CR_F = FR_F / (Q_in * (FR_F / FR_Fin) * (T / TR_in) * ((P * y_F) / (PR_in * y_Fin)));

% Entalpía de Reacción:

% Paquete de Ecuaciones: El Diferencial se encuentra despejado para ser solucionado:
     Ec = [0; 0; 0; 0; 0; 0];
     % Ecuación 1: Balance de Materia de Metanol (A):
     Ec(1) = 0;
     % Ecuación 2: Balance de Materia de Oxígeno (B):
     Ec(2) = 0;
     % Ecuación 3: Balance de Energía de Fluido de Reacción: Diferencial dT / dx
     Ec(3) = (h_masa * a_p * Ro_p * (T - Ts) + U * pi() * di / A_trans) / (FR_A * Cp_A(T) + FR_B * Cp_B(T) + FR_C(FR_A) * Cp_c(T) ...
              + FR_D(FR_A, FR_B) * Cp_D(T) + FR_E(FR_A, FR_B)* Cp_D(T) + FR_F * Cp_F(T));
     % Ecuación 4: Balance de Energía Interfacial:
     Ec(4) = h_calor - a_p * (Ts - T) - n1 * Hreac1 * r1 - n2 * Hreac2 * r2;
     % Ecuación 5: Balance de Energía de Fluido de Control: Diferencial dT / dx
     Ec(5) = -U * pi() * do * (T_Gin - T_Gout) / (Q_G * Ro_G * Cp_G);
     % Ecuación 6: Balance de Presión:
     Ec(6) = 0;
     % Conjunto de Ecuaciones:
     Ec = [Ec(1); Ec(2); Ec(3); Ec(4); Ec(5); Ec(6);];
     % Balance de Energía de Fluido de Control:
end
    














