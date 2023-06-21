
% A: Metanol (CH3OH)
% B: Oxígeno (O2)
% C: Formaldehído (CH2O)
% D: Agua (H2O)
% E: Monóxido de Carbono (CO)
% F: Nitrógeno (N2)
% G: Fluido de Control

% VALORES CONSTANTES / INICIALES:

% Pesos moleculares [kg / mol]:
     PM_A = 32.04 / 1000; PM_B = 32.00 / 1000; PM_C = 30.03 / 1000; PM_D = 18.00 / 1000; PM_E = 28.01 / 1000; PM_F = 28.01 / 1000;
     
% Parámetros de Ingreso:
     % Fluido de Reacción
     TR_in = 573; % [K]
     PR_in = 2.5; % [atm]
     ConstR = 0.08205746 * 1000; % [atm m3 / mol K]
     
     y_Ain = 0.2; y_Cin = 0.001; y_Din = 0.002; y_Ein = 0.008; y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21; 
          y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;

     fR_in = 0.6; % [kg / s]
     FR_in = fR_in / (y_Ain * PM_A + y_Bin * PM_B + y_Cin * PM_C + y_Din * PM_D + y_Ein * PM_E + y_Fin * PM_F); % [mol / s]
     FR_Ain = FR_in * y_Ain; FR_Bin = FR_in * y_Bin; FR_Cin = FR_in * y_Cin; FR_Din = FR_in * y_Din; FR_Ein = FR_in * y_Ein; 
          FR_Fin = FR_in * y_Fin; % [mol / s]  
     fr_Ain = FR_Ain * PM_A; fr_Bin = FR_Bin * PM_B; fr_Cin = FR_Cin * PM_C; fr_Din = FR_Din * PM_D; fr_Ein = FR_Ein * PM_E;
          fr_Fin = FR_Fin * PM_F; % [kg / s]

     % Fluido de Control
     Q_G = 
     Ro_G =
     U = 200;
     T_Gin = 273; % [K]
     T_Gout = 350; % [K]

     % Características de Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * di^2 * 0.25;
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 10; % [m]

     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     Vacio = 0.390 * 1.740 / (di / dp + 1.140)^2;
     NR_s = 4.5; % [kg / m2 s]
     Vel_s = NR_s;

% RESOLUCIÓN DE SISTEMA DE ECUACIONES CON SOLVER ODE15S: https://www.mathworks.com/help/matlab/ref/ode15s.html
     tspan = linspace(0, Lmax, 5000); % Intervalo de evaluación
     CI = [0, FR_A, FR_E, T_rin, T_cin, T_sin]; % Condiciones Iniciales
     M = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]; % Matriz recomendada
     Options = odeset('Mass', M, 'RelTol', 1e-4, 'AbsTol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6]); % Configuración propuesta
     Solucion = ode15s(@(x, FR_A, FR_E, Tr, Tc, Ts)Derivadas(), tspan, CI, options);
          % Para ecuaciones Diferenciales: dZ / dx = ...
          % Para expresiones Algebráicas: 0 = ...

% PLANTEAMIENTO DE SOLUCIÓN DE ECUACIONES:
function Ec = Derivadas()

% Calores Específicos y deltas CP [KJ /mol K]: 
%    Coeficientes obtenidos de "Principios Elementales de los Procesos Químicos - Felder" - Tabla B.2  
     C1_A = (+42.93 * 10^-3); C2_A = (+8.3010 * 10^-5); C3_A = (-1.8700 * 10^-8); C4_A = (-8.030 * 10^-12);
     C1_B = (+29.10 * 10^-3); C2_B = (+1.1580 * 10^-5); C3_B = (-0.6076 * 10^-8); C4_B = (+1.311 * 10^-12);
     C1_C = (+34.28 * 10^-3); C2_C = (+4.2680 * 10^-5); C3_C = (+0.0000 * 10^-8); C4_C = (-8.694 * 10^-12);
     C1_D = (+33.46 * 10^-3); C2_D = (+0.6880 * 10^-5); C3_D = (+0.7604 * 10^-8); C4_D = (-3.593 * 10^-12);
     C1_E = (+28.95 * 10^-3); C2_E = (+0.4110 * 10^-5); C3_E = (+0.3548 * 10^-8); C4_E = (-2.220 * 10^-12);
     C1_F = (+29.00 * 10^-3); C2_F = (+0.2199 * 10^-5); C3_F = (+0.5723 * 10^-8); C4_F = (-2.871 * 10^-12);
     C1_G = (0); C2_G = (0; C3_G = (0); C4_G = (0);
     
     Cp_A = @(T) C1_A + C2_A * (T-273.15) + C3_A * (T-273.15)^2 + C4_A * (T-273.15)^3;
     Cp_B = @(T) C1_B + C2_B * (T-273.15) + C3_B * (T-273.15)^2 + C4_B * (T-273.15)^3;
     Cp_C = @(T) C1_C + C2_C * (T-273.15) + C3_C * (T-273.15)^2 + C4_C * (T-273.15)^3;
     Cp_D = @(T) C1_D + C2_D * (T-273.15) + C3_D * (T-273.15)^2 + C4_D * (T-273.15)^3;
     Cp_E = @(T) C1_E + C2_E * (T-273.15) + C3_E * (T-273.15)^2 + C4_E * (T-273.15)^3;
     Cp_F = @(T) C1_F + C2_F * (T-273.15) + C3_F * (T-273.15)^2 + C4_F * (T-273.15)^3;
     Cp_G = @(T) 1; % ( * 10^-3) + ( * 10^-5)*(T-273.15) + ( * 10^-8) * (T-273.15)^2 + ( * 10^-12) * (T-273.15)^3;

     DCP_1 = @(T) (C1_D + C1_C - 0.5 * C1_B - C1_A) + (C2_D + C2_C - 0.5 * C2_B - C2_A) * (T-273.15) + ...
                  (C3_D + C3_C - 0.5 * C3_B - C3_A) * (T-273.15)^2 + (C4_D + C4_C - 0.5 * C4_B - C4_A) * (T-273.15)^3;
     DCP_2 = @(T) (C1_E + C1_D - 0.5 * C1_B - C1_C) + (C2_E + C2_D - 0.5 * C2_B - C2_C) * (T-273.15) + ...
                  (C3_E + C3_D - 0.5 * C3_B - C3_C) * (T-273.15)^2 + (C4_E + C4_D - 0.5 * C4_B - C4_C) * (T-273.15)^3;

% Cinéticas de Reacción: 
     % Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539;
     Eact_k1 = 62661.3851;
     A_K1A = 62661.3851;
     Eact_K1A = 133441.5273;
     A_k2 = 0.017971;
     Eact_k2 = 32820.4093;
     A_K2D = 1.1538;
     Eact_K2D = -16817.5167;

     k1 = @(T) A_k1 * exp(-Eact_k1 / (ConstR * T));
     K1A = @(T) A_K1A * exp(-Eact_K1A / (ConstR * T));
     k2 = @(T) A_k2 * exp(-Eact_k2 / (ConstR * T));
     K2D = @(T) A_K2D * exp(-Eact_K2D / (ConstR * T));

     n1 = 1;
     n2 = 1;

% Relaciones Estequiométricas y Flujos Molares:
     % Variables Independientes: FR_A, FR_E:
     E1 = @(FR_A) FR_Ain - FR_A; % Avance de Reacción 1 en función de Metanol (A)
     E2 = @(FR_E) FR_Ein - FR_E; % Avance de Reacción 2 en función de Monóxido de Carbono (E)
     FR_B = @(FR_A, FR_E) FR_Bin - 0.5 * E1(FR_A) - 0.5 * E2(FR_E); 
     FR_C = @(FR_A, FR_E) FR_Cin + E1(FR_A) - E2(FR_E); % Concentración de Formaldehído (C)
     FR_D = @(FR_A, FR_E) FR_Din + E1(FR_A) + E2(FR_E); % Concentración de Agua (D)
     FR_F = FR_Fin; % Concentración de Nitrógeno (F)
     FR = FR_A + FR_B + FR_C + FR_D + FR_E + FR_F;

% Flujos (f), Fluxes (N) y Concentraciones (x) Másicos:
     fR_A = FR_A * PM_A; fR_B = FR_B * PM_B; fR_C = FR_C * PM_C; fR_D = FR_D * PM_D; fR_E = FR_E * PM_E; fR_F = FR_F * PM_F;
     fR = fR_A + fR_B + fR_C + fR_D + fR_E + fR_D;

     NR_A = fR_A / A_transi; NR_B = fR_B / A_transi; NR_C = fR_C / A_transi; NR_D = fR_D / A_transi; NR_E = fR_E / A_transi; 
     NR_F = fR_F / A_transi; NR = NR_A + NR_B + NR_C + NR_D + NR_E + NR_F; 

     x_A = fr_A / fr; x_B = fr_B / fr; x_C = fr_C / fr; x_D = fr_D / fr; x_E = fr_E / fr; x_F = fr_F / fr;

     Q = Q_in * (FR / FR_in) * (T / T_in) * (PR_in / PR);

% Concentraciones (y) molares de gas, caudales parciales, presiones parciales y densidades:
     y_A = FR_A / FR; y_B = FR_B / FR; y_C = FR_C / FR; y_D = FR_D / FR; y_E = FR_E / FR; y_F = FR_F / FR;
     
     Q_A = Q * y_A; Q_B = Q * y_B; Q_C = Q * y_C; Q_D = Q * y_D; Q_E = Q * y_E; Q_F = Q * y_F;

     PR_A = PR * y_A; PR_B = PR * y_B; PR_C = PR * y_C; PR_D = PR * y_D; PR_E = PR * y_E; PR_F = PR * y_F;
     
     Ro_A = fR_A / Q_A; Ro_B = fR_B / Q_B; Ro_C = fR_C / Q_C; Ro_D = fR_D / Q_D; Ro_E = fR_E / Q_E; Ro_F = fR_F / Q_F;
     Ro_R = 1 / (x_A / Ro_A + x_B / Ro_B + x_C / Ro_C + x_D / Ro_D + x_E / Ro_E + x_F / Ro_F);

     Vel_R = Q / A_transi;

% Conversión de Flujos a Concentraciones [mol / m3]:
     CR_A = FR_A / Q; CR_B = FR_B / Q; CR_C = FR_C / Q; CR_D = FR_D / Q; CR_E = FR_E / Q; CR_F = FR_F / Q;

% Viscosidades [Pa * s]:
     Mu_A = @(T) 0.00003771 * T + (-0.002267); Mu_B = @(T) 0.00004980 * T + (+0.005360); Mu_C = @(T) 0.00003771 * T + (-0.002267);
     Mu_D = @(T) 0.00003710 * T + (-0.001140); Mu_E = @(T) 0.00004160 * T + (+0.005100); Mu_F = @(T) 0.00004160 * T + (+0.005100);
     Mu_R = @(T) x_A * Mu_A + x_B * Mu_B + x_C * Mu_C + x_D * Mu_D + x_E * Mu_E + x_F * Mu_F;

% Entalpía de Reacción:

% Caída de presión:
     f = (1 - Vacio) / (Vacio^3) * (1.75 + (150 * (1 - Vacio)) / (dp * NR / Mu_R));

% Paquete de Ecuaciones: El Diferencial se encuentra despejado para ser solucionado:
     Ec = [0; 0; 0; 0; 0; 0];
     % Ecuación 1: Balance de Materia de Metanol (A) - Se encuentra dFA / dx igualado a la expresión:
     Ec(1) = (k1(T_r) * (PR_in * y_Ain (Q_Ain / Q_A) * (F_A / F_Ain) * (T_r/T_rin))) / (1 + K1_A * (PR_in * y_Ain * (Q_Ain / Q_A) * (F_A / F_Ain) * (T_r/T_rin)));
     
     % Ecuación 2: Balance de Materia de Monóxido de Carbono (E) - Se encuentra dFE / dx igualado a la expresión:
     Ec(2) = (k2(T_r) * (PR_in * y_Cin * (Q_Cin / Q_CA) * (F_C / F_Cin) * (T_r/T_rin))) / (1 + K1_A * (PR_in * y_Din * (Q_Din / Q_D) * (F_D / F_Din) * (T_r/T_rin)));
     
     % Ecuación 3: Balance de Energía de Fluido de Reacción - Se encuentra dTr / dx igualado a la expresión:
     Ec(3) = (h_masa * a_p * Ro_p * (T - Ts) + U * pi() * di / A_trans) / (FR_A * Cp_A(T) + FR_B * Cp_B(T) + FR_C(FR_A) * Cp_c(T) ...
              + FR_D(FR_A, FR_B) * Cp_D(T) + FR_E(FR_A, FR_B)* Cp_D(T) + FR_F * Cp_F(T));

     % Ecuación 4: Balance de Energía Interfacial - Se encuentra 0 igualado a la expresión para manejar dTs / dx:
     Ec(4) = 0;

     % Ecuación 5: Balance de Energía de Fluido de Control: - Se encuentra dTc / dx igualado a la expresión:
     Ec(5) = (-U * pi() * do * (T_Gin - T_Gout)) / (Q_G * Ro_G * Cp_G);

     % Ecuación 6: Balance de Presión  - Se encuentra dPR / dz igualado a la expresión
     Ec(6) = -f * Ro_R * Vel_s; dp;

     % Conjunto de Ecuaciones:
     Ec = [Ec(1); Ec(2); Ec(3); Ec(4); Ec(5); Ec(6)];

end



