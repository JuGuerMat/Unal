% A: Metanol (CH3OH)
% B: Oxígeno (O2)
% C: Formaldehído (CH2O)
% D: Agua (H2O)
% E: Monóxido de Carbono (CO)
% F: Nitrógeno (N2)
% G: Fluido de Control

% Pesos moleculares [kg / mol]:
     PM_A = 32.04 / 1000; PM_B = 32.00 / 1000; PM_C = 30.03 / 1000; PM_D = 18.00 / 1000; PM_E = 28.01 / 1000; PM_F = 28.01 / 1000;
     
% Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539;
     Eact_k1 = 62661.3851;
     A_K1A = 62661.3851;
     Eact_K1A = 133441.5273;
     A_k2 = 0.017971;
     Eact_k2 = 32820.4093;
     A_K2D = 1.1538;
     Eact_K2D = -16817.5167;

% Parámetros de Ingreso:
     % Fluido de Reacción
     T_rin = 573; % [K]
     PR_in = 2.5 * 101.3 * 1000; % [Pa]
     ConstR = 8.31446261815324; % [J / mol K]
     dP_in = 0.1;
     
     y_Ain = 0.2; y_Cin = 0.001; y_Din = 0.002; y_Ein = 0.008; y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21; 
          y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;

     fR_in = 0.6; % [kg / s]
     FR_in = fR_in / (y_Ain * PM_A + y_Bin * PM_B + y_Cin * PM_C + y_Din * PM_D + y_Ein * PM_E + y_Fin * PM_F); % [mol / s]
     FR_Ain = FR_in * y_Ain; FR_Bin = FR_in * y_Bin; FR_Cin = FR_in * y_Cin; FR_Din = FR_in * y_Din; FR_Ein = FR_in * y_Ein; 
          FR_Fin = FR_in * y_Fin; % [mol / s]  
     
     Q_in = FR_in * T_rin / (PR_in);
          Q_Ain = Q_in * y_Ain;
          Q_Bin = Q_in * y_Bin;
          Q_Cin = Q_in * y_Cin;
          Q_Din = Q_in * y_Din;
          Q_Ein = Q_in * y_Ein;
          Q_Fin = Q_in * y_Fin;

     fr_Ain = FR_Ain * PM_A; fr_Bin = FR_Bin * PM_B; fr_Cin = FR_Cin * PM_C; fr_Din = FR_Din * PM_D; fr_Ein = FR_Ein * PM_E;
          fr_Fin = FR_Fin * PM_F; % [kg / s]

     % Fluido de Control: HP Hytherm 600 (https://www.hplubricants.in/products/specialties/thermic-fluids/hytherm-500-and-600-thermic-fluid-oil)
     f_G = 0.5; % [kg / s];
     U = 200;
     T_cin = 300; % [K]
     T_cout = 350; % [K]
     CondTerm_G = 0.090 / 1.1629999998093; % [W / m K]

     % Características de Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * di^2 * 0.25;
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 10; % [m]

     % Características del Sistema del Lecho:
     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     ap = ((dp/2)^2 / (3 * (dp/2)^3))*1/Ro_p;
     Vacio = 0.390 * 1.740 / (di / dp + 1.140)^2;
     NR_s = 4.5; % [kg / m2 s]
     Vel_s = NR_s;
     T_sin = 560; % [K]
     h_masa = 50;
     Ro_Lecho = 1500;

% RESOLUCIÓN DE SISTEMA DE ECUACIONES CON SOLVER ODE15S: https://www.mathworks.com/help/matlab/ref/ode15s.html
     Lmax = 5;
     tspan = linspace(0, Lmax, 20); % Intervalo de evaluación
     CI = [FR_Ain; FR_Ein; T_rin; T_cin; T_sin; dP_in]; % Condiciones Iniciales
     %CI = [1, 1, 1, 1, 1, 1]
     M = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0]; % Matriz recomendada
     options = odeset('Mass', M, 'RelTol', 1e-4, 'AbsTol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6]); % Configuración propuesta
     Sol = ode15s(@Derivadas, tspan, CI, options);
          % Para ecuaciones Diferenciales: dZ / dx = ...
          % Para expresiones Algebráicas: 0 = ...

% PLANTEAMIENTO DE SOLUCIÓN DE ECUACIONES:
function Ec = Derivadas(x, y)
% VALORES CONSTANTES / INICIALES:

% Pesos moleculares [kg / mol]:
     PM_A = 32.04 / 1000; PM_B = 32.00 / 1000; PM_C = 30.03 / 1000; PM_D = 18.00 / 1000; PM_E = 28.01 / 1000; PM_F = 28.01 / 1000;
     
% Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539;
     Eact_k1 = 62661.3851;
     A_K1A = 62661.3851;
     Eact_K1A = 133441.5273;
     A_k2 = 0.017971;
     Eact_k2 = 32820.4093;
     A_K2D = 1.1538;
     Eact_K2D = -16817.5167;

% Parámetros de Ingreso:
     % Fluido de Reacción
     T_rin = 573; % [K]
     PR_in = 2.5; % [atm]
     ConstR = 0.08205746 * 1000; % [atm m3 / mol K]
     dP_in = 0;
     
     y_Ain = 0.2; y_Cin = 0.001; y_Din = 0.002; y_Ein = 0.008; y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21; 
          y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;

     fR_in = 0.6; % [kg / s]
     FR_in = fR_in / (y_Ain * PM_A + y_Bin * PM_B + y_Cin * PM_C + y_Din * PM_D + y_Ein * PM_E + y_Fin * PM_F); % [mol / s]
     FR_Ain = FR_in * y_Ain; FR_Bin = FR_in * y_Bin; FR_Cin = FR_in * y_Cin; FR_Din = FR_in * y_Din; FR_Ein = FR_in * y_Ein; 
          FR_Fin = FR_in * y_Fin; % [mol / s]  
     
     Q_in = FR_in * T_rin / (PR_in);
          Q_Ain = Q_in * y_Ain;
          Q_Bin = Q_in * y_Bin;
          Q_Cin = Q_in * y_Cin;
          Q_Din = Q_in * y_Din;
          Q_Ein = Q_in * y_Ein;
          Q_Fin = Q_in * y_Fin;

     fr_Ain = FR_Ain * PM_A; fr_Bin = FR_Bin * PM_B; fr_Cin = FR_Cin * PM_C; fr_Din = FR_Din * PM_D; fr_Ein = FR_Ein * PM_E;
          fr_Fin = FR_Fin * PM_F; % [kg / s]

     % Fluido de Control: HP Hytherm 600 (https://www.hplubricants.in/products/specialties/thermic-fluids/hytherm-500-and-600-thermic-fluid-oil)
     f_G = 0.5; % [kg / s];
     U = 200;
     T_cin = 273; % [K]
     T_cout = 350; % [K]
     CondTerm_G = 0.090 / 1.1629999998093; % [W / m K]

     % Características de Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * di^2 * 0.25;
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 10; % [m]

     % Características del Sistema del Lecho:
     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     ap = ((dp/2)^2 / (3 * (dp/2)^3))*1/Ro_p;
     Vacio = 0.390 * 1.740 / (di / dp + 1.140)^2;
     NR_s = 4.5; % [kg / m2 s]
     Vel_s = NR_s;
     T_sin = 573; % [K]
     h_masa = 100;
     a_p = 100;
     Ro_Lecho = 1500;

% Variables: t: x, y(1): FR_A, y(2): FR_E, y(3) : Tr, y(4) : Tc, y(5): Ts, y(6) : dP 
     y = [0, 0, 0 , 0 , 0, 0];
     FR_A = y(1);
     FR_E = y(2);
     T_r = y(3);
     T_c = y(4);
     T_s = y(5);
     dP = y(6);
     PR = @(dP) (PR_in - dP);

% Coeficientes de Capacidad Calorífica. Obtenidos de "Principios Elementales de los Procesos Químicos - Felder" - Tabla B.2 (Todos menos G)
     C1_A = (+42.93 * 10^-3); C2_A = (+8.3010 * 10^-5); C3_A = (-1.8700 * 10^-8); C4_A = (-8.030 * 10^-12);
     C1_B = (+29.10 * 10^-3); C2_B = (+1.1580 * 10^-5); C3_B = (-0.6076 * 10^-8); C4_B = (+1.311 * 10^-12);
     C1_C = (+34.28 * 10^-3); C2_C = (+4.2680 * 10^-5); C3_C = (+0.0000 * 10^-8); C4_C = (-8.694 * 10^-12);
     C1_D = (+33.46 * 10^-3); C2_D = (+0.6880 * 10^-5); C3_D = (+0.7604 * 10^-8); C4_D = (-3.593 * 10^-12);
     C1_E = (+28.95 * 10^-3); C2_E = (+0.4110 * 10^-5); C3_E = (+0.3548 * 10^-8); C4_E = (-2.220 * 10^-12);
     C1_F = (+29.00 * 10^-3); C2_F = (+0.2199 * 10^-5); C3_F = (+0.5723 * 10^-8); C4_F = (-2.871 * 10^-12);
     C1_G = (0.184355); C2_G = (+1.025 * 10^-3);

% Calores Específicos y deltas CP [J /mol K]:
     Cp_A = @(T) (C1_A + C2_A * (T-273.15) + C3_A * (T-273.15)^2 + C4_A * (T-273.15)^3) * 1000;
     Cp_B = @(T) (C1_B + C2_B * (T-273.15) + C3_B * (T-273.15)^2 + C4_B * (T-273.15)^3) * 1000;
     Cp_C = @(T) (C1_C + C2_C * (T-273.15) + C3_C * (T-273.15)^2 + C4_C * (T-273.15)^3) * 1000;
     Cp_D = @(T) (C1_D + C2_D * (T-273.15) + C3_D * (T-273.15)^2 + C4_D * (T-273.15)^3) * 1000;
     Cp_E = @(T) (C1_E + C2_E * (T-273.15) + C3_E * (T-273.15)^2 + C4_E * (T-273.15)^3) * 1000;
     Cp_F = @(T) (C1_F + C2_F * (T-273.15) + C3_F * (T-273.15)^2 + C4_F * (T-273.15)^3) * 1000;
     Cp_G = @(T) C1_G + C2_G * (T); 

     DCP_1 = @(T) (C1_D + C1_C - 0.5 * C1_B - C1_A) + (C2_D + C2_C - 0.5 * C2_B - C2_A) * (T-273.15) + ...
                  (C3_D + C3_C - 0.5 * C3_B - C3_A) * (T-273.15)^2 + (C4_D + C4_C - 0.5 * C4_B - C4_A) * (T-273.15)^3;
     DCP_2 = @(T) (C1_E + C1_D - 0.5 * C1_B - C1_C) + (C2_E + C2_D - 0.5 * C2_B - C2_C) * (T-273.15) + ...
                  (C3_E + C3_D - 0.5 * C3_B - C3_C) * (T-273.15)^2 + (C4_E + C4_D - 0.5 * C4_B - C4_C) * (T-273.15)^3;

% Entalpía de Reacción [J / mol]:
     EntR1 = @(T) ((-115.90 - 241.83) - (-201.20 + 0.5 * 0)) + ...
                  ((C1_D + C1_C - 0.5 * C1_B - C1_A) * (T-273.15) + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (T-273.15)^2 + ...
                  (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (T-273.15)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (T-273.15)^4) - ...
                  (((C1_D + C1_C - 0.5 * C1_B - C1_A) * (25) + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (25)^2 + ...
                  (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (25)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (25)^4));

     EntR2 = @(T) ((-110.52 - 241.83) - (-115.90 + 0.5 * 0)) + ...
                  ((C1_E + C1_D - 0.5 * C1_B - C1_C) * (T-273.15) + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (T-273.15)^2 + ...
                  (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (T-273.15)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (T-273.15)^4) - ...
                  ((C1_E + C1_D - 0.5 * C1_B - C1_C) * (25) + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (25)^2 + ...
                  (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (25)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (25)^4);

% Cinéticas de Reacción: 
     k1 = @(T) A_k1 * exp(-Eact_k1 / (ConstR * T));
     K1A = @(T) A_K1A * exp(-Eact_K1A / (ConstR * T));
     k2 = @(T) A_k2 * exp(-Eact_k2 / (ConstR * T));
     K2D = @(T) A_K2D * exp(-Eact_K2D / (ConstR * T));

% Relaciones Estequiométricas y Flujos Molares:
     % Variables Independientes: FR_A, FR_E:
     E1 = @(FR_A) FR_Ain - FR_A; % Avance de Reacción 1 en función de Metanol (A)
     E2 = @(FR_E) FR_Ein - FR_E; % Avance de Reacción 2 en función de Monóxido de Carbono (E)
     FR_B = @(FR_A, FR_E) FR_Bin - 0.5 * E1(FR_A) - 0.5 * E2(FR_E); 
     FR_C = @(FR_A, FR_E) FR_Cin + E1(FR_A) - E2(FR_E); % Concentración de Formaldehído (C)
     FR_D = @(FR_A, FR_E) FR_Din + E1(FR_A) + E2(FR_E); % Concentración de Agua (D)
     FR_F = FR_Fin; % Concentración de Nitrógeno (F)
     FR = @(FR_A, FR_E) FR_A + FR_B(FR_A, FR_E) + FR_C(FR_A, FR_E) + FR_D(FR_A, FR_E) + FR_E + FR_F;

% Flujos (f), Fluxes (N) y Concentraciones (x) Másicos:
     fR_A = @(FR_A) FR_A * PM_A; fR_B = @(FR_A, FR_E) FR_B(FR_A, FR_E) * PM_B; fR_C = @(FR_A, FR_E) FR_C(FR_A, FR_E) * PM_C; 
          fR_D = @(FR_A, FR_E) FR_D(FR_A, FR_E) * PM_D; fR_E = @(FR_E) FR_E * PM_E; fR_F = FR_F * PM_F;
     fR = @(FR_A, FR_E) fR_A(FR_A) + fR_B(FR_A, FR_E) + fR_C(FR_A, FR_E) + fR_D(FR_A, FR_E) + fR_E(FR_E) + fR_F;

     NR_A = @(FR_A) fR_A(FR_A) / A_transi; NR_B = @(FR_A, FR_E) fR_B(FR_A, FR_E) / A_transi; NR_C = @(FR_A, FR_E) fR_C(FR_A, FR_E) / A_transi; 
          NR_D = @(FR_A, FR_E) fR_D(FR_A, FR_E) / A_transi; NR_E = @(FR_E) fR_E(FR_E) / A_transi; NR_F = fR_F / A_transi; 
     NR = @(FR_A, FR_E) NR_A(FR_A) + NR_B(FR_A, FR_E) + NR_C(FR_A, FR_E) + NR_D(FR_A, FR_E) + NR_E(FR_E) + NR_F; 

     x_A = @(FR_A, FR_E) fR_A(FR_A) / fR(FR_A, FR_E); x_B = @(FR_A, FR_E) fR_B(FR_A, FR_E) / fR(FR_A, FR_E); 
          x_C = @(FR_A, FR_E) fR_C(FR_A, FR_E) / fR(FR_A, FR_E); x_D = @(FR_A, FR_E) fR_D(FR_A, FR_E) / fR(FR_A, FR_E); 
          x_E = @(FR_A, FR_E) fR_E(FR_E) / fR(FR_A, FR_E); x_F = @(FR_A, FR_E) fR_F / fR(FR_A, FR_E);

     Q = @(FR_A, FR_E, T_r, dP) Q_in * (FR(FR_A, FR_E) / FR_in) * (T_r / T_rin) * (PR_in / PR(dP));

% Concentraciones (y) molares de gas, caudales parciales, presiones parciales y densidades:
     y_A = @(FR_A, FR_E) FR_A / FR(FR_A, FR_E); y_B = @(FR_A, FR_E) FR_B(FR_A, FR_E) / FR(FR_A, FR_E); 
          y_C = @(FR_A, FR_E) FR_C(FR_A, FR_E) / FR(FR_A, FR_E); y_D = @(FR_A, FR_E) FR_D(FR_A, FR_E) / FR(FR_A, FR_E); 
          y_E = @(FR_A, FR_E) FR_E / FR(FR_A, FR_E); y_F = @(FR_A, FR_E) FR_F / FR(FR_A, FR_E);
     
     Q_A = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_A(FR_A, FR_E); Q_B = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_B(FR_A, FR_E); 
          Q_C = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_C(FR_A, FR_E); Q_D = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_D(FR_A, FR_E); 
          Q_E = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_E(FR_A, FR_E); Q_F = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) * y_F(FR_A, FR_E);

     PR_A = @(FR_A, FR_E, dP) PR(dP) * y_A(FR_A, FR_E); PR_B = @(FR_A, FR_E, dP) PR(dP) * y_B(FR_A, FR_E); 
          PR_C = @(FR_A, FR_E, dP) PR(dP) * y_C(FR_A, FR_E); PR_D = @(FR_A, FR_E, dP) PR(dP) * y_D(FR_A, FR_E); 
          PR_E = @(FR_A, FR_E, dP) PR(dP) * y_E(FR_A, FR_E); PR_F = @(FR_A, FR_E, dP) PR(dP) * y_F(FR_A, FR_E);
     
     Ro_A = @(FR_A, FR_E, T_r, dP) fR_A(FR_A) / Q_A(FR_A, FR_E, T_r, dP); Ro_B = @(FR_A, FR_E, T_r, dP) fR_B(FR_A, FR_E) / Q_B(FR_A, FR_E, T_r, dP); 
          Ro_C = @(FR_A, FR_E, T_r, dP) fR_C(FR_A, FR_E) / Q_C(FR_A, FR_E, T_r, dP); Ro_D = @(FR_A, FR_E, T_r, dP) fR_D(FR_A, FR_E) / Q_D(FR_A, FR_E, T_r, dP); 
          Ro_E = @(FR_A, FR_E, T_r, dP) fR_E(FR_E) / Q_E(FR_A, FR_E, T_r, dP); Ro_F = @(FR_A, FR_E, T_r, dP) fR_F / Q_F(FR_A, FR_E, T_r, dP);
     Ro_R = 1 / (x_A(FR_A, FR_E) / Ro_A(FR_A, FR_E, T_r, dP) + x_B(FR_A, FR_E) / Ro_B(FR_A, FR_E, T_r, dP) + ...
                 x_C(FR_A, FR_E) / Ro_C(FR_A, FR_E, T_r, dP) + x_D(FR_A, FR_E) / Ro_D(FR_A, FR_E, T_r, dP) + ...
                 x_E(FR_A, FR_E) / Ro_E(FR_A, FR_E, T_r, dP) + x_F(FR_A, FR_E) / Ro_F(FR_A, FR_E, T_r, dP));

     Vel_R = @(FR_A, FR_E, T_r, dP) Q(FR_A, FR_E, T_r, dP) / A_transi;

% Conversión de Flujos a Concentraciones [mol / m3]:
     CR_A = @(FR_A, FR_E, T_r, dP) FR_A / Q(FR_A, FR_E, T_r, dP); CR_B = @(FR_A, FR_E, T_r, dP) FR_B(FR_A, FR_E) / Q(FR_A, FR_E, T_r, dP); 
          CR_C = @(FR_A, FR_E, T_r, dP) FR_C(FR_A, FR_E) / Q(FR_A, FR_E, T_r, dP); CR_D = @(FR_A, FR_E, T_r, dP) FR_D(FR_A, FR_E) / Q(FR_A, FR_E, T_r, dP); 
          CR_E = @(FR_A, FR_E, T_r, dP) FR_E / Q(FR_A, FR_E, T_r, dP); CR_F = @(FR_A, FR_E, T_r, dP) FR_F(FR_A, FR_E) / Q(FR_A, FR_E, T_r, dP);

% Viscosidades [Pa * s]:
     Mu_A = @(T_r) 0.00003771 * T_r + (-0.002267); Mu_B = @(T_r) 0.00004980 * T_r + (+0.005360); Mu_C = @(T_r) 0.00003771 * T_r + (-0.002267);
     Mu_D = @(T_r) 0.00003710 * T_r + (-0.001140); Mu_E = @(T_r) 0.00004160 * T_r + (+0.005100); Mu_F = @(T_r) 0.00004160 * T_r + (+0.005100);
     Mu_R = @(FR_A, FR_E, T_r) 1 / (x_A(FR_A, FR_E) / Mu_A(T_r) + x_B(FR_A, FR_E) / Mu_B(T_r) + x_C(FR_A, FR_E) / Mu_C(T_r) + ...
               x_D(FR_A, FR_E) / Mu_D(T_r) + x_E(FR_A, FR_E) / Mu_E(T_r) + x_F(FR_A, FR_E) / Mu_F(T_r));

     Mu_G = 0.0279;

% Caída de presión:
     f = @(FR_A, FR_E, T_r) (1 - Vacio) / (Vacio^3) * (1.75 + (150 * (1 - Vacio)) / (dp * NR(FR_A, FR_E) / Mu_R(FR_A, FR_E, T_r)));

% Paquete de Ecuaciones: El Diferencial se encuentra despejado para ser solucionado:
     Ec = [0; 0; 0; 0; 0; 0];
     % Ecuación 1: Balance de Materia de Metanol (A) - Se encuentra dFA / dx igualado a la expresión:
     Ec(1) = (k1(y(3)) * (PR_in * y_Ain * (Q_Ain / Q_A(y(1), y(2), T_r, dP)) * (y(1) / FR_Ain) * (y(3)/T_rin))) / (1 + K1A(y(3)) * (PR_in * y_Ain * ...
             (Q_Ain / Q_A(y(1), y(2), T_r, dP)) * (y(1) / FR_Ain) * (y(3)/T_rin)));
     
     % Ecuación 2: Balance de Materia de Monóxido de Carbono (E) - Se encuentra dFE / dx igualado a la expresión:
     Ec(2) = (k2(y(3)) * (PR_in * y_Cin * (Q_Cin / Q_C(y(1), y(2), T_r, dP)) * FR_C(y(1), y(2)) / FR_Cin * (y(3)/T_rin))) / (1 + K2D(y(3)) * (PR_in * y_Din * ...
             (Q_Din / Q_D(y(1), y(2), T_r, dP)) * (FR_D(y(1), y(2)) / FR_Din) * (y(3)/T_rin)));
     
     % Ecuación 3: Balance de Energía de Fluido de Reacción - Se encuentra dTr / dx igualado a la expresión:
     Ec(3) = (h_masa * a_p * Ro_Lecho * (y(3) - y(5)) + U * pi() * di / A_transi) / (y(1) * Cp_A(y(3)) + FR_B(y(1), y(2)) * Cp_B(y(3)) + ...
             FR_C(y(1), y(2)) * Cp_C(y(3)) + FR_D(y(1), y(2)) * Cp_D(y(3)) + y(2) * Cp_E(y(3)) + FR_F * Cp_F(y(3)));

     % Ecuación 4: Balance de Energía Interfacial - Se encuentra 0 igualado a la expresión para manejar dTs / dx:
     n1 = 1; n2 = 1; x1 = 1; x2 = 1; x3 = 1; x4 = 1; x5 = 1;
     Ec(4) = h_masa * a_p * Ro_p * (y(5) - y(3)) - ...
             n1 * EntR1(y(5)) * (x1 * y(5)^2 + x2 * y(5) + x3) * ((k1(y(5)) * (PR_in * y_Ain * (Q_Ain / Q_A(y(1), y(2), T_r, dP)) * (y(1) / FR_Ain) * (y(3)/T_rin))) / ...
             (1 + K1A(y(5)) * (PR_in * y_Ain * (Q_Ain / Q_A(y(1), y(2), T_r, dP)) * (y(1) / FR_Ain) * (y(3)/T_rin)))) - ...
             n2 * EntR2(y(5)) * (x3 * y(5)^2 + x4 * y(5) + x5) * ((k2(y(5)) * (PR_in * y_Cin * (Q_Cin / Q_C(y(1), y(2), T_r, dP)) * FR_C(y(1), y(2)) / FR_Cin * (y(3)/T_rin))) / ...
             (1 + K2D(y(5)) * (PR_in * y_Din * (Q_Din / Q_D(y(1), y(2), T_r, dP)) * (FR_D(y(1), y(2)) / FR_Din) * (y(3)/T_rin))));

     % Ecuación 5: Balance de Energía de Fluido de Control: - Se encuentra dTc / dx igualado a la expresión:
     Ec(5) = (-U * pi() * do * (y(4) - T_cin)) / (f_G * Cp_G(y(3)));

     % Ecuación 6: Balance de Presión  - Se encuentra dPR / dz igualado a la expresión
     Ec(6) = -f(y(1), y(2), y(3)) * Ro_R * Vel_s / dp;

     % Conjunto de Ecuaciones:
     Ec = [Ec(1); Ec(2); Ec(3); Ec(4); Ec(5); Ec(6)];

end

% Pendientes:
% - Coeficientes y parámetros adicionales
% - Balances adicionales

