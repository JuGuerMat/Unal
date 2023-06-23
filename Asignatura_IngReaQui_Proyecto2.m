
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
     
% Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539; Eact_k1 = 62661.3851;
     A_K1A = 62661.3851; Eact_K1A = 133441.5273;
     A_k2 = 0.017971; Eact_k2 = 32820.4093;
     A_K2D = 1.1538; Eact_K2D = -16817.5167;

% Parámetros de Ingreso (En x = 0):
     % Condiciones generales iniciales para fluido de reacción:
     T_rin = 573; % [K]
     PR_in = 2.5 * 101.3 * 1000; % [Pa]
     ConstR = 8.31446261815324; % [J / mol K]
     dP_in = 0.1;
     
     % Concentraciones molares iniciales:
     y_Ain = 0.2; y_Cin = 0.001; y_Din = 0.002; y_Ein = 0.008; 
     y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21; y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;

     % Flujos másicos (fr [kg / s]) y molares (FR [mol / s]) iniciales:
     fR_in = 0.6; 

     FR_in = fR_in / (y_Ain * PM_A + y_Bin * PM_B + y_Cin * PM_C + y_Din * PM_D + y_Ein * PM_E + y_Fin * PM_F); 
     FR_Ain = FR_in * y_Ain; FR_Bin = FR_in * y_Bin; FR_Cin = FR_in * y_Cin; FR_Din = FR_in * y_Din; FR_Ein = FR_in * y_Ein; 
     FR_Fin = FR_in * y_Fin;

     fr_Ain = FR_Ain * PM_A; fr_Bin = FR_Bin * PM_B; fr_Cin = FR_Cin * PM_C; 
     fr_Din = FR_Din * PM_D; fr_Ein = FR_Ein * PM_E; fr_Fin = FR_Fin * PM_F; 
     
     % Caudales (Q) iniciales:
     Q_in = FR_in * T_rin / (PR_in);
     Q_Ain = Q_in * y_Ain; Q_Bin = Q_in * y_Bin; Q_Cin = Q_in * y_Cin;
     Q_Din = Q_in * y_Din; Q_Ein = Q_in * y_Ein; Q_Fin = Q_in * y_Fin;

% Fluido de Control: 
     % Referencia: HP Hytherm 600 (https://www.hplubricants.in/products/specialties/thermic-fluids/hytherm-500-and-600-thermic-fluid-oil)
     f_G = 0.5; % [kg / s];
     U = 200;
     T_cin = 300; % [K]
     T_cout = 400; % [K]
     CondTerm_G = 0.090 / 1.1629999998093; % [W / m K]
     Mu_G = 0.0279;

% Características de Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * di^2 * 0.25;
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 5; % [m]

% Características del Sistema del Lecho:
     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     ap = ((dp/2)^2 / (3 * (dp/2)^3))*1/Ro_p;
     Vacio = 0.390 * 1.740 / (di / dp + 1.140)^2;
     NR_s = 4.5; % [kg / m2 s]
     Vel_s = NR_s * A_transi;
     T_sin = 600; % [K]
     h_masa = 100;
     Ro_Lecho = 1500;

% ESTABLECIMIENTO RESOLUCIÓN DE SISTEMA DE ECUACIONES CON SOLVER ODE15S: https://www.mathworks.com/help/matlab/ref/ode15s.html
     y = [0, 0, 0, 0, 0, 0];
     Lmax = 5;
     tspan = linspace(0, Lmax, 2000); % Intervalo de evaluación
     % CI = [FR_Ain; FR_Ein; T_rin; T_sin; T_cin; PR_in]; % Condiciones Iniciales
     CI = [FR_Ain; FR_Ein; T_rin; T_cin; PR_in]; % Condiciones Iniciales
     y = CI;
     %CI = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]';
     M = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]; % Matriz recomendada
     options = odeset('Mass', M, 'RelTol', 1e-2, 'AbsTol', [1e-2 1e-2 1e-2 1e-2 1e-2]); % Configuración propuesta
     [X Y] = ode15s(@Derivadas, tspan, CI, options);
          % Para ecuaciones Diferenciales: dZ / dx = ...
          % Para expresiones Algebráicas: 0 = ...

     for i = 1 : length(Y')
          VFRA(i) = Y(i, 1);
          VFRE(i) = Y(i, 2);
          VTR(i) = Y(i, 3);
          VTC(i) = Y(i, 4);
          VPR(i) = Y(i, 5);
     end

     subplot(3, 2, 1);
     plot(X, VFRA);
     hold on;
     subplot(3, 2, 2);
     plot(X, VFRE);
     hold on;
     subplot(3, 2, 3);
     plot(X, VTR);
     hold on;
     subplot(3, 2, 4);
     plot(X, VTC);
     hold on;
     subplot(3, 2, 5);
     plot(X, VTS);
     hold on;
     subplot(3, 2, 6);
     plot(X, VPR);
     hold on; 
     

% PLANTEAMIENTO DE SOLUCIÓN DE ECUACIONES:
function Ec = Derivadas(x, y)

% VALORES CONSTANTES / INICIALES:

% Pesos moleculares [kg / mol]:
     PM_A = 32.04 / 1000; PM_B = 32.00 / 1000; PM_C = 30.03 / 1000; PM_D = 18.00 / 1000; PM_E = 28.01 / 1000; PM_F = 28.01 / 1000;
     
% Parámetros obtenidas con las expresiones de Cinética B (Denominador más corto), [Eact = J/mol]:
     A_k1 = 121.5539; Eact_k1 = 62661.3851;
     A_K1A = 62661.3851; Eact_K1A = 133441.5273;
     A_k2 = 0.017971; Eact_k2 = 32820.4093;
     A_K2D = 1.1538; Eact_K2D = -16817.5167;

% Parámetros de Ingreso (En x = 0):
     % Condiciones generales iniciales para fluido de reacción:
     T_rin = 573; % [K]
     PR_in = 2.5 * 101.325 * 1000; % [Pa]
     ConstR = 8.31446261815324; % [J / mol K]
     
     % Concentraciones molares iniciales:
     y_Ain = 0.2; y_Cin = 0.001; y_Din = 0.02; y_Ein = 0.008; 
     y_Bin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.21; y_Fin = (1 - y_Ain - y_Cin - y_Din - y_Ein) * 0.79;

     % Flujos másicos (fr [kg / s]) y molares (FR [mol / s]) iniciales:
     fR_in = 0.6; 

     FR_in = fR_in / (y_Ain * PM_A + y_Bin * PM_B + y_Cin * PM_C + y_Din * PM_D + y_Ein * PM_E + y_Fin * PM_F); 
     FR_Ain = FR_in * y_Ain; FR_Bin = FR_in * y_Bin; FR_Cin = FR_in * y_Cin; FR_Din = FR_in * y_Din; FR_Ein = FR_in * y_Ein; 
     FR_Fin = FR_in * y_Fin;

     fr_Ain = FR_Ain * PM_A; fr_Bin = FR_Bin * PM_B; fr_Cin = FR_Cin * PM_C; 
     fr_Din = FR_Din * PM_D; fr_Ein = FR_Ein * PM_E; fr_Fin = FR_Fin * PM_F; 
     
     % Caudales (Q) iniciales:
     Q_in = (FR_in * T_rin * ConstR) / (PR_in);
     Q_Ain = Q_in * y_Ain; Q_Bin = Q_in * y_Bin; Q_Cin = Q_in * y_Cin;
     Q_Din = Q_in * y_Din; Q_Ein = Q_in * y_Ein; Q_Fin = Q_in * y_Fin;

% Fluido de Control: 
     % Referencia: HP Hytherm 600 (https://www.hplubricants.in/products/specialties/thermic-fluids/hytherm-500-and-600-thermic-fluid-oil)
     f_G = 0.5; % [kg / s]
     U = 200; % [W / m2 s]
     T_cin = 300; % [K]
     T_cout = 350; % [K]
     CondTerm_G = 0.090 / 1.1629999998093; % [W / m K]
     Mu_G = 0.0279; %[Pa s]

% Características Dadas para Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * ((di^2) / 4); % [m]
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     Lmax = 5; % [m]

% Características del Sistema del Lecho:
     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     ap = (6 / dp) * 1/Ro_p; % [m2 particula / m3 lecho]
     Vacio = 0.390 + (1.740 / ((di/dp) + 1.140)^2);
     NR_s = 4.5; % [kg / m2 s]
     f_s = NR_s * A_transi; % [kg / s]
     T_sin = 300; % [K]
     h_calor = 0.01;
     Ro_Lecho = (1 - Vacio) * Ro_p; % [kg / m3]

% PARÁMETROS VARIABLES
% Variables: t: x, y(1): FR_A, y(2): FR_E, y(3) : Tr, y(4) : Tc, y(5): Ts, y(6) : P 
     FR_A = y(1);
     FR_C = y(2);
     T_r = y(3);
     T_c = y(4);
     dP = y(5);
     PR = dP;

     T_s = 580;

% Coeficientes de Capacidad Calorífica. Obtenidos de "Principios Elementales de los Procesos Químicos - Felder" - Tabla B.2 (Todos menos G)
     C1_A = (+42.93 * 10^-3); C2_A = (+8.3010 * 10^-5); C3_A = (-1.8700 * 10^-8); C4_A = (-8.030 * 10^-12); % [KJ / mol K]
     C1_B = (+29.10 * 10^-3); C2_B = (+1.1580 * 10^-5); C3_B = (-0.6076 * 10^-8); C4_B = (+1.311 * 10^-12); % [KJ / mol K]
     C1_C = (+34.28 * 10^-3); C2_C = (+4.2680 * 10^-5); C3_C = (+0.0000 * 10^-8); C4_C = (-8.694 * 10^-12); % [KJ / mol K]
     C1_D = (+33.46 * 10^-3); C2_D = (+0.6880 * 10^-5); C3_D = (+0.7604 * 10^-8); C4_D = (-3.593 * 10^-12); % [KJ / mol K]
     C1_E = (+28.95 * 10^-3); C2_E = (+0.4110 * 10^-5); C3_E = (+0.3548 * 10^-8); C4_E = (-2.220 * 10^-12); % [KJ / mol K]
     C1_F = (+29.00 * 10^-3); C2_F = (+0.2199 * 10^-5); C3_F = (+0.5723 * 10^-8); C4_F = (-2.871 * 10^-12); % [KJ / mol K]
     C1_G = (0.184355); C2_G = (+1.025 * 10^-3); % [Kg / mol K]

% Calores Específicos y deltas CP [J /mol K]:
     Cp_A = (C1_A + C2_A * (T_r-273.15) + C3_A * (T_r-273.15)^2 + C4_A * (T_r-273.15)^3) * 1000;
     Cp_B = (C1_B + C2_B * (T_r-273.15) + C3_B * (T_r-273.15)^2 + C4_B * (T_r-273.15)^3) * 1000;
     Cp_C = (C1_C + C2_C * (T_r-273.15) + C3_C * (T_r-273.15)^2 + C4_C * (T_r-273.15)^3) * 1000;
     Cp_D = (C1_D + C2_D * (T_r-273.15) + C3_D * (T_r-273.15)^2 + C4_D * (T_r-273.15)^3) * 1000;
     Cp_E = (C1_E + C2_E * (T_r-273.15) + C3_E * (T_r-273.15)^2 + C4_E * (T_r-273.15)^3) * 1000;
     Cp_F = (C1_F + C2_F * (T_r-273.15) + C3_F * (T_r-273.15)^2 + C4_F * (T_r-273.15)^3) * 1000;
     Cp_G =  C1_G + C2_G * (y(5)); 

     DCP_1 = (C1_D + C1_C - 0.5 * C1_B - C1_A) + (C2_D + C2_C - 0.5 * C2_B - C2_A) * (T_r-273.15) + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A) * (T_r-273.15)^2 + (C4_D + C4_C - 0.5 * C4_B - C4_A) * (T_r-273.15)^3;
     DCP_2 = (C1_E + C1_D - 0.5 * C1_B - C1_C) + (C2_E + C2_D - 0.5 * C2_B - C2_C) * (T_r-273.15) + ...
             (C3_E + C3_D - 0.5 * C3_B - C3_C) * (T_r-273.15)^2 + (C4_E + C4_D - 0.5 * C4_B - C4_C) * (T_r-273.15)^3;

% Entalpía de Reacción [J / mol]:
     EntR1 = ((-115.90 - 241.83) - (-201.20 + 0.5 * 0)) + ...
             ((C1_D + C1_C - 0.5 * C1_B - C1_A) * (T_r-273.15) + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (T_r-273.15)^2 + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (T_r-273.15)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (T_r-273.15)^4) - ...
             (((C1_D + C1_C - 0.5 * C1_B - C1_A) * (25) + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (25)^2 + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (25)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (25)^4));

     EntR2 = ((-110.52 - 241.83) - (-115.90 + 0.5 * 0)) + ...
             ((C1_E + C1_D - 0.5 * C1_B - C1_C) * (T_r-273.15) + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (T_r-273.15)^2 + ...
             (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (T_r-273.15)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (T_r-273.15)^4) - ...
             ((C1_E + C1_D - 0.5 * C1_B - C1_C) * (25) + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (25)^2 + ...
             (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (25)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (25)^4);

% Cinéticas de Reacción: 
     k1 = A_k1 * exp(-Eact_k1 / (ConstR * T_s));
     K1A = A_K1A * exp(-Eact_K1A / (ConstR * T_s));
     k2 = A_k2 * exp(-Eact_k2 / (ConstR * T_s));
     K2D = A_K2D * exp(-Eact_K2D / (ConstR * T_s));

% Relaciones Estequiométricas y Flujos Molares:
     % Variables Independientes: FR_A, FR_E:
     FR_A = y(1);
     FR_E = y(2);

     E1 = FR_Ain - FR_A; % Avance de Reacción 1 en función de Metanol (A)
     E2 = FR_Ein - FR_E; % Avance de Reacción 2 en función de Monóxido de Carbono (
     FR_B = FR_Bin - 0.5 * E1 - 0.5 * E2; % Concentración de Oxígeno (B)
     FR_C = FR_Cin + E1 - E2; % Concentración de Formaldehído (C)
     FR_D = FR_Din + E1 + E2; % Concentración de Agua (D)  
     FR_F = FR_Fin; % Concentración de Nitrógeno (F)
     FR = FR_A + FR_B + FR_C + FR_D + FR_E + FR_F;

% Flujos másicos de reacción (fR) [kg / s]:
     fR_A = FR_A * PM_A; 
     fR_B = FR_B * PM_B; 
     fR_C = FR_C * PM_C; 
     fR_D = FR_D * PM_D; 
     fR_E = FR_E * PM_E; 
     fR_F = FR_F * PM_F;
     fR = fR_A + fR_B + fR_C + fR_D + fR_E + fR_F;

% Fluxes Másicos de reacción (NR) [kg / m2 s]:
     NR_A = fR_A / A_transi; 
     NR_B = fR_B / A_transi; 
     NR_C = fR_C / A_transi; 
     NR_D = fR_D / A_transi; 
     NR_E = fR_E / A_transi; 
     NR_F = fR_F / A_transi; 
     NR = NR_A + NR_B + NR_C + NR_D + NR_E + NR_F; 

% Concentraciones másicas de reacción (x):
     x_A = fR_A / fR; 
     x_B = fR_B / fR; 
     x_C = fR_C / fR; 
     x_D = fR_D / fR; 
     x_E = fR_E / fR; 
     x_F = fR_F / fR;

% Concentraciones molares de reacción (y):
     y_A = FR_A / FR; 
     y_B = FR_B / FR; 
     y_C = FR_C / FR; 
     y_D = FR_D / FR; 
     y_E = FR_E / FR; 
     y_F = FR_F / FR;

% Caudales de reacción (Q) [m3 / s]:
     Q = Q_in * (FR / FR_in) * (T_r / T_rin) * (PR_in / PR);
     
     Q_A = Q * y_A; 
     Q_B = Q * y_B; 
     Q_C = Q * y_C; 
     Q_D = Q * y_D; 
     Q_E = Q * y_E; 
     Q_F = Q * y_F;

     % Velocidad de Flujo de Reacción:
     Vel_R = Q / A_transi;

% Presiones parciales de reacción (PR):
     PR_A = PR * y_A; 
     PR_B = PR * y_B; 
     PR_C = PR * y_C;
     PR_D = PR * y_D; 
     PR_E = PR * y_E; 
     PR_F = PR * y_F;

% Densidades de reacción (Ro)
     Ro_A = fR_A / Q_A; 
     Ro_B = fR_B / Q_B; 
     Ro_C = fR_C / Q_C;
     Ro_D = fR_D / Q_D; 
     Ro_E = fR_E / Q_E; 
     Ro_F = fR_F / Q_F;

     Ro_R = 1 / (x_A / Ro_A + x_B / Ro_B + x_C / Ro_C + x_D / Ro_D + x_E / Ro_E + x_F / Ro_F);

     Vel_s = NR_s * Ro_R;

% Conversión de Flujos Molares a Concentraciones [mol / m3]:
     CR_A = FR_A / Q; 
     CR_B = FR_B / Q; 
     CR_C = FR_C / Q; 
     CR_D = FR_D / Q; 
     CR_E = FR_E / Q; 
     CR_F = FR_F / Q;

% Viscosidades de Flujo de Reacción [Pa * s]:
     Mu_A = 0.00003771 * T_r + (-0.002267); 
     Mu_B = 0.00004980 * T_r + (+0.005360); 
     Mu_C = 0.00003771 * T_r + (-0.002267);
     Mu_D = 0.00003710 * T_r + (-0.001140); 
     Mu_E = 0.00004160 * T_r + (+0.005100); 
     Mu_F = 0.00004160 * T_r + (+0.005100);
     
     Mu_R = 1 / (x_A / Mu_A + x_B / Mu_B + x_C/ Mu_C + ...
               x_D / Mu_D + x_E / Mu_E + x_F / Mu_F);

% Factor de Caída de Presión:
     f = (1 - Vacio) / (Vacio^3) * (1.75 + (150 * (1 - Vacio)) / (dp * NR / Mu_R));

% PAQUETE DE ECUACIONES DIFERENCIALES Y ALGEBRÁICAS:

% Resolución de expresión algebráica:
     n1 = 1; n2 = 1; 

     syms Balance(T_s)

     Balance(T_s) = h_calor * ap * Ro_p * (T_s - T_r) - ...
           n1 * EntR1 * (k1 * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) / ...
           (1 + K1A * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) - ...
           n2 * EntR2 * ((k2) * (PR_in * y_Cin * (Q_Cin / Q_C) * (FR_C / FR_Cin) * (y(3)/T_rin))) / ...
           (1 + K2D * (PR_in * y_Din * (Q_Din / Q_D) * (FR_D / FR_Din) * (y(3)/T_rin)));

     T_s = vpasolve(Balance);

% Paquete de Ecuaciones: El Diferencial se encuentra despejado para ser solucionado:
% Variables: t: x, y(1): FR_A, y(2): FR_E, y(3) : Tr, y(4) : Tc, y(5): Ts, y(6) : P 
     Ec = zeros(5,1);
     % Ecuación 1: Balance de Materia de Metanol (A) - Se encuentra dFA / dx igualado a la expresión:
     Ec1 = (k1 * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) / ...
             (1 + K1A * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) * Ro_R * A_transi;

     % Ecuación 2: Balance de Materia de Monóxido de Carbono (E) - Se encuentra dFE / dx igualado a la expresión:
     Ec2 = (k2 * (PR_in * y_Cin * (Q_Cin / Q_C) * (FR_C / FR_Cin) * (y(3)/T_rin))) / ...
             (1 + K2D * (PR_in * y_Din * (Q_Din / Q_D) * (FR_D / FR_Din) * (y(3)/T_rin))) * Ro_R * A_transi;

     % Ecuación 3: Balance de Energía de Fluido de Reacción - Se encuentra dTr / dx igualado a la expresión:
     Ec3 = (-h_calor * ap * Ro_Lecho * (y(3) - y(5)) * A_transi + U * pi() * do * (y(4) - y(3))) / ...
             (FR_A * Cp_A + FR_B * Cp_B + FR_C * Cp_C + FR_D * Cp_D + FR_E * Cp_E + FR_F * Cp_F);

     % Ecuación 4: Balance de Energía de Fluido de Control: - Se encuentra dTc / dx igualado a la expresión:
     Ec4 = (U * pi() * do * (y(4) - y(3))) / (f_G * Cp_G);

     % Ecuación 5: Balance de Presión - Se encuentra dPR / dz igualado a la expresión
     Ec5 = -f * Ro_R * Vel_s / dp;

     % Conjunto de Ecuaciones:
     Ec = [Ec1; Ec2; Ec3; Ec4; Ec5];

end

