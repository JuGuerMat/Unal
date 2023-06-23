% Institución: Universidad Nacional de Colombia
%
% Nombres: - María Paula Bravo Castillo, mbravoc@unal.edu.co 
%          - Juan Pablo Chitiva Arteaga, jchitiva@unal.edu.co
%          - Paula Alejandra Duarte Salamanca, pduarte@unal.edu.co
%          - Juan José Guerrero Acosta, juguerrero@unal.edu.co
%          - Daniel Andrés López Brijaldo, dlopezbr@unal.edu.co
%
% Asignatura: Ingeniería de las Reacciones Químicas
%
% Docente: Hugo Martín Galiendo Valbuena
%
% Título: Proyecto Final "Diseño de Reactor" - Código en MatLab de Solución de Ecuaciones
%
% Fecha: 23 de junio de 2023

% 1) NOMENCLATURA DE SUSTANCIAS QUÍMICAS:
     % A: Metanol (CH3OH)
     % B: Oxígeno (O2)
     % C: Formaldehído (CH2O)
     % D: Agua (H2O)
     % E: Monóxido de Carbono (CO)
     % F: Nitrógeno (N2)
     % G: Fluido de Control

% 2) DECLARACIÓN DE VALORES CONSTANTES E INICIALES:

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
     Lmax = 20; % [m]

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

% 3) ESTABLECIMIENTO RESOLUCIÓN DE SISTEMA DE ECUACIONES CON SOLVER ODE15S: 
     % https://www.mathworks.com/help/matlab/ref/ode15s.html
     y = [0, 0, 0, 0, 0, 0];
     Lmax = 2.5;
     tspan = linspace(0, Lmax, 20000); % Intervalo de evaluación
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
          VFR_A(i) = Y(i, 1); % Flujo de Metanol
          VFR_E(i) = Y(i, 2); % Flujo de CO
          VT_r(i) = Y(i, 3); % Temperatura de R
          VT_c(i) = Y(i, 4); % Temperatura de C
          VPR(i) = Y(i, 5) / 1000; % Presión del FLuido de Reaccion
     end

% 4) MANIPULACIÓN DE VECTORES RESULTANTES Y ELABORACIÓN DE GRÁFICAS:
     % Vectores de Variables Independientes: FR_A, FR_E:
     FR_A = y(1);
     FR_E = y(2);

     % Vectores de Flujos, Concentraciones, Temperaturas y Presión:
     for i = 1 : length (VFR_A)
          VE1(i) = FR_Ain - VFR_A(i); % Avance de Reacción 1 en función de Metanol (A)
          VE2(i) = FR_Ein - VFR_E(i); % Avance de Reacción 2 en función de Monóxido de Carbono (
          VFR_B(i) = FR_Bin - 0.5 * VE1(i) - 0.5 * VE2(i); % Concentración de Oxígeno (B)
          VFR_C(i) = FR_Cin + VE1(i) - VE2(i); % Concentración de Formaldehído (C)
          VFR_D(i) = FR_Din + VE1(i) + VE2(i); % Concentración de Agua (D)  
          VFR_F(i) = FR_Fin; % Concentración de Nitrógeno (F)
          VFR_R(i) = VFR_A(i) + VFR_B(i) + VFR_C(i) + VFR_D(i) + VFR_E(i) + VFR_F(i);

          VQ(i) = Q_in * (VFR_R(i) / FR_in) * (VT_r(i) / T_rin) * (PR_in / VPR(i));

          VConv_A(i) = (FR_Ain - VFR_A(i)) / FR_Ain;
          VRend_C(i) = (VE1(i) - VE2(i)) / (VE1(i));

          VConc_A(i) = VFR_A(i) / VQ(i);
          VConc_B(i) = VFR_B(i) / VQ(i);
          VConc_C(i) = VFR_C(i) / VQ(i);
          VConc_D(i) = VFR_D(i) / VQ(i);
          VConc_E(i) = VFR_E(i) / VQ(i);
          VConc_F(i) = VFR_F(i) / VQ(i);
          VConc_R(i) = VFR_R(i) / VQ(i);
     end 

     % Elección de tipo de gráficas a elaborar: (1 Cuantificación de Materia, 2 Temperaturas y Presion): 
     Graficas = 2;
     switch(Graficas)
          case 1
               subplot(2, 2, 1);
               plot(X, VConv_A, '-', 'linewidth', 2);
               hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Conversión")
               title("Conversión de Metanol Respecto a Longitud de Reactor")
               grid on;
          
               subplot(2, 2, 2);
               plot(X, VRend_C, '-', 'linewidth', 2);
               hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Rendimiento")
               title("Rendimiento de Formaldehído Respecto a Longitud de Reactor")
               grid on;
          
               subplot(2, 2, 3);
               plot(X, VFR_A, '-', 'linewidth', 2); hold on;
               plot(X, VFR_B, '-', 'linewidth', 2); hold on;
               plot(X, VFR_C, '-', 'linewidth', 2); hold on;
               plot(X, VFR_D, '-', 'linewidth', 2); hold on;
               plot(X, VFR_E, '-', 'linewidth', 2); hold on;
               plot(X, VFR_F, '-', 'linewidth', 2); hold on;
               plot(X, VFR_R, '-', 'linewidth', 2); hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Flujo molar [mol /s]")
               title("Flujos Molares de Reacción Respecto a Longitud de Reactor")
               legend("Metanol", "Oxígeno", "Formaldehído", "Agua", "Monóxido de Carbono", "Nitrógeno", "Flujo Total")
               grid on;
          
               subplot(2, 2, 4);
               plot(X, VConc_A, '-', 'linewidth', 2); hold on;
               plot(X, VConc_B, '-', 'linewidth', 2); hold on;
               plot(X, VConc_C, '-', 'linewidth', 2); hold on;
               plot(X, VConc_D, '-', 'linewidth', 2); hold on;
               plot(X, VConc_E, '-', 'linewidth', 2); hold on;
               plot(X, VConc_F, '-', 'linewidth', 2); hold on;
               plot(X, VConc_R, '-', 'linewidth', 2); hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Concentracion molar [mol / m3]")
               title("Concentraciones de Flujo de Reacción Respecto a Longitud de Reactor")
               legend("Metanol", "Oxígeno", "Formaldehído", "Agua", "Monóxido de Carbono", "Nitrógeno", "Flujo Total")
               grid on;
          case 2
               subplot(1, 2, 1)
               plot(X, VT_r, 'linewidth', 2); hold on;
               plot(X, VT_c, 'linewidth', 2); hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Temperatura [K]")
               title("Temperaturas de Flujos de Reacción y de Control Respecto a Longitud de Reactor")
               legend("Fluido de Reacción", "Fluido de Control")
               grid on;

               subplot(1, 2, 2)
               plot(X, VPR, 'linewidth', 2); hold on;
               xlabel("Distancia de Reactor [m]")
               ylabel("Presión [kpa]")
               title("Presión en las Tuberías de Flujo de Reacción Respecto a Longitud de Reactor")
               grid on;
     end

% 5. PLANTEAMIENTO DE SOLUCIÓN DE ECUACIONES:
function Ec = Derivadas(x, y)

% 5.A VALORES CONSTANTES / INICIALES (NUEVAMENTE PARA EL INTERIOR DE LA FUNCIÓN):

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
     T_cin = 300; % [K]
     T_cout = 350; % [K]
     CondTerm_G = 0.090 / 1.1629999998093; % [W / m K]
     Mu_G = 0.0279; %[Pa s]

% Características Dadas para Diseño de Reactor:
     di = 2.0574 / 100; % [m]
     A_transi = pi() * ((di^2) / 4); % [m]
     do = 2.54 / 100; % [m]
     d = do - di; % [m]
     MLDd = (do - di) / log(do / di);
     Lmax = 5; % [m]
     Censu = 0.01 / 5.67826334106816;
     CondTerm_Tubo = 16.20; % [W / m2 K]

% Características del Sistema del Lecho:
     dp = 0.4 / 100; % [m]
     Ro_p = 1.892 * 1000; % [kg / m3]
     ap = (6 / dp) * 1/Ro_p; % [m2 particula / m3 lecho]
     Vacio = 0.390 + (1.740 / ((di/dp) + 1.140)^2);
     NR_G = 4.5; % [kg / m2 s]
     f_s = NR_G * A_transi; % [kg / s]
     T_sin = 300; % [K]
     h_calor = 0.01;
     Ro_Lecho = (1 - Vacio) * Ro_p; % [kg / m3]

% 5.B ESTABLECIMIENTO DE VARIABLES:
% Variables: t: x, y(1): FR_A, y(2): FR_E, y(3) : Tr, y(4) : Tc, y(5): Ts, y(6) : P 
     FR_A = y(1);
     FR_C = y(2);
     T_r = y(3);
     T_c = y(4);
     dP = y(5);
     PR = dP;

     T_s = 580;

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
     NR_R = NR_A + NR_B + NR_C + NR_D + NR_E + NR_F; 

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
     PM_R = y_A * PM_A + y_B * PM_B + y_C * PM_C  + y_D * PM_D + y_E * PM_E + y_F * PM_F;

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

     Vel_s = NR_G * Ro_R;

% Conversión de Flujos Molares a Concentraciones [mol / m3]:
     CR_A = FR_A / Q; 
     CR_B = FR_B / Q; 
     CR_C = FR_C / Q; 
     CR_D = FR_D / Q; 
     CR_E = FR_E / Q; 
     CR_F = FR_F / Q;

% Viscosidades de Flujo de Reacción [Pa * s]:
     Mu_A = 3 * 10^-8 * T_r + 5 * 10^-7;
     Mu_B = 2 * 10^-9 * T_r^2 - 2 * 10-7 * T_r + 1 * 10^-5; 
     Mu_C = 3 * 10^-8 * T_r + 5 * 10^-7;
     Mu_D = 4 * 10^-9 * T_r - 2 * 10^-6; 
     Mu_E = 1 * 10^-9 * T_r^2 - 2 * 10^-7 * T_r + 9 * 10^-6; 
     Mu_F = 2 * 10^-9 * T_r^2 - 2 * 10^-7 * T_r + 1 * 10^-5;

     Mu_R = 1 / (x_A / Mu_A + x_B / Mu_B + x_C/ Mu_C + x_D / Mu_D + x_E / Mu_E + x_F / Mu_F);

% Conductividades:
     CondTerm_A = 5 * 10^-7 * T_r^2 - 0.0002 * T_r + 0.0284;
     CondTerm_B = 1 * 10^-5 * T_r^2 - 0.0018 * T_r + 0.0758;
     CondTerm_C = 5 * 10^-7 * T_r^2 - 0.0002 * T_r + 0.0284;
     CondTerm_D = 1 * 10^-6 * T_r^2 - 0.0011 * T_r + 0.2011;
     CondTerm_E = 9 * 10^-6 * T_r^2 - 0.0014 * T_r + 0.0588;
     CondTerm_F = 7 * 10^-6 * T_r^2 - 0.0010 * T_r + 0.0462;

     CondTerm_R = x_A / CondTerm_A + x_B / CondTerm_B + x_C / CondTerm_C + x_D / CondTerm_D + x_E / CondTerm_E + x_F / CondTerm_F;

% Factor de Caída de Presión:
     f = (1 - Vacio) / (Vacio^3) * (1.75 + (150 * (1 - Vacio)) / (dp * NR_R / Mu_R));
     Re_R = (dp * NR_R / Mu_R);

% Coeficientes de Capacidad Calorífica. Obtenidos de "Principios Elementales de los Procesos Químicos - Felder" - Tabla B.2 (Todos menos G)
     %{
     C1_A = (+42.93 * 10^-3); C2_A = (+8.3010 * 10^-5); C3_A = (-1.8700 * 10^-8); C4_A = (-8.030 * 10^-12); % [KJ / mol K]
     C1_B = (+29.10 * 10^-3); C2_B = (+1.1580 * 10^-5); C3_B = (-0.6076 * 10^-8); C4_B = (+1.311 * 10^-12); % [KJ / mol K]
     C1_C = (+34.28 * 10^-3); C2_C = (+4.2680 * 10^-5); C3_C = (+0.0000 * 10^-8); C4_C = (-8.694 * 10^-12); % [KJ / mol K]
     C1_D = (+33.46 * 10^-3); C2_D = (+0.6880 * 10^-5); C3_D = (+0.7604 * 10^-8); C4_D = (-3.593 * 10^-12); % [KJ / mol K]
     C1_E = (+28.95 * 10^-3); C2_E = (+0.4110 * 10^-5); C3_E = (+0.3548 * 10^-8); C4_E = (-2.220 * 10^-12); % [KJ / mol K]
     C1_F = (+29.00 * 10^-3); C2_F = (+0.2199 * 10^-5); C3_F = (+0.5723 * 10^-8); C4_F = (-2.871 * 10^-12); % [KJ / mol K]
     C1_G = (0.184355); C2_G = (+1.025 * 10^-3); % [Kg / mol K]
     %}
     
     C1_A = (+26.84900); C2_A = (+0.075000); C3_A = (-0.000010); C4_A = (0); C5_A = (0); % [J / mol K]
     C1_B = (+31.32234); C2_B = (-20.23531); C3_B = (+57.86644); C4_B = (-36.50624); C5_B = (-0.007374); % [J / mol K]
     C1_C = (+5.193767); C2_C = (+93.23249); C3_C = (-44.85457); C4_C = (+7.882279); C5_C = (+0.551175); % [J / mol K]
     C1_D = (+30.09200); C2_D = (+6.832514); C3_D = (+6.793435); C4_D = (-2.534480); C5_D = (+0.082139); % [J / mol K]
     C1_E = (+25.56759); C2_E = (+6.096130); C3_E = (+4.054656); C4_E = (-2.671021); C5_E = (+0.131021);% [J / mol K]
     C1_F = (+19.50583); C2_F = (+19.88705); C3_F = (-8.598535); C4_F = (+1.369784); C5_F = (+0.527601); % [J / mol K]
     C1_G = (+0.184355); C2_G = (+1.025 * 10^-3); % [J / kg K]

     
% Calores Específicos y deltas CP [J /mol K]:
     
     %{
     Cp_A = (C1_A + C2_A * (T_r-273.15) + C3_A * (T_r-273.15)^2 + C4_A * (T_r-273.15)^3) * 1000;
     Cp_B = (C1_B + C2_B * (T_r-273.15) + C3_B * (T_r-273.15)^2 + C4_B * (T_r-273.15)^3) * 1000;
     Cp_C = (C1_C + C2_C * (T_r-273.15) + C3_C * (T_r-273.15)^2 + C4_C * (T_r-273.15)^3) * 1000;
     Cp_D = (C1_D + C2_D * (T_r-273.15) + C3_D * (T_r-273.15)^2 + C4_D * (T_r-273.15)^3) * 1000;
     Cp_E = (C1_E + C2_E * (T_r-273.15) + C3_E * (T_r-273.15)^2 + C4_E * (T_r-273.15)^3) * 1000;
     Cp_F = (C1_F + C2_F * (T_r-273.15) + C3_F * (T_r-273.15)^2 + C4_F * (T_r-273.15)^3) * 1000;
     Cp_G =  C1_G + C2_G * (y(5));
     Cp_R = 1 / (y_A / Cp_A + y_B / Cp_B + y_C / Cp_C + y_D / Cp_D + y_E / Cp_E + y_F / Cp_F);

     DCP_1 = (C1_D + C1_C - 0.5 * C1_B - C1_A) + (C2_D + C2_C - 0.5 * C2_B - C2_A) * (T_r-273.15) + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A) * (T_r-273.15)^2 + (C4_D + C4_C - 0.5 * C4_B - C4_A) * (T_r-273.15)^3;
     DCP_2 = (C1_E + C1_D - 0.5 * C1_B - C1_C) + (C2_E + C2_D - 0.5 * C2_B - C2_C) * (T_r-273.15) + ...
             (C3_E + C3_D - 0.5 * C3_B - C3_C) * (T_r-273.15)^2 + (C4_E + C4_D - 0.5 * C4_B - C4_C) * (T_r-273.15)^3;
     %}

     Cp_A = (C1_A + C2_A * (T_r / 1000) + C3_A * (T_r / 1000)^2 + C4_A * (T_r / 1000)^3 + C5_A * (T_r / 1000)^-2);
     Cp_B = (C1_B + C2_B * (T_r / 1000) + C3_B * (T_r / 1000)^2 + C4_B * (T_r / 1000)^3 + C5_B * (T_r / 1000)^-2);
     Cp_C = (C1_C + C2_C * (T_r / 1000) + C3_C * (T_r / 1000)^2 + C4_C * (T_r / 1000)^3 + C5_C * (T_r / 1000)^-2);
     Cp_D = (C1_D + C2_D * (T_r / 1000) + C3_D * (T_r / 1000)^2 + C4_D * (T_r / 1000)^3 + C5_D * (T_r / 1000)^-2);
     Cp_E = (C1_E + C2_E * (T_r / 1000) + C3_E * (T_r / 1000)^2 + C4_E * (T_r / 1000)^3 + C5_E * (T_r / 1000)^-2);
     Cp_F = (C1_F + C2_F * (T_r / 1000) + C3_F * (T_r / 1000)^2 + C4_F * (T_r / 1000)^3 + C5_F * (T_r / 1000)^-2);
     Cp_G =  C1_G + C2_G * (T_r);
     Cp_R = 1 / (y_A / Cp_A + y_B / Cp_B + y_C / Cp_C + y_D / Cp_D + y_E / Cp_E + y_F / Cp_F);

     DCP_1 = (C1_D + C1_C - 0.5 * C1_B - C1_A) * (T_r / 1000)^0 + (C2_D + C2_C - 0.5 * C2_B - C2_A) * (T_r / 1000)^1 + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A) * (T_r / 1000)^2 + (C4_D + C4_C - 0.5 * C4_B - C4_A) * (T_r / 1000)^3 + ... 
             (C5_D + C5_C - 0.5 * C5_B - C5_A) * (T_r / 1000)^-2;
     DCP_2 = (C1_E + C1_D - 0.5 * C1_B - C1_C) + (T_r / 1000)^0 + (C2_E + C2_D - 0.5 * C2_B - C2_C) * (T_r / 1000)^1 + ... 
             (C3_E + C3_D - 0.5 * C3_B - C3_C) * (T_r / 1000)^2 + (C4_E + C4_D - 0.5 * C4_B - C4_C) * (T_r / 1000)^3 + ...
             (C5_E + C5_D - 0.5 * C5_B - C5_C) * (T_r / 1000)^-2;

% Entalpía de Reacción [J / mol]:
     EntR1 = (((-115.90 - 241.83) - (-201.20 + 0.5 * 0))  + ...
            ((C1_D + C1_C - 0.5 * C1_B - C1_A)/1 * (T_r / 1000)^1 + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (T_r / 1000)^2 + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (T_r / 1000)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (T_r / 1000)^4 + ... 
            ((C5_D + C5_C - 0.5 * C5_B - C5_A)/-1 * (T_r / 1000)^-1) - ...
             (C1_D + C1_C - 0.5 * C1_B - C1_A)/1 * (T_r / 1000)^1 + (C2_D + C2_C - 0.5 * C2_B - C2_A)/2 * (T_r / 1000)^2 + ...
             (C3_D + C3_C - 0.5 * C3_B - C3_A)/3 * (T_r / 1000)^3 + (C4_D + C4_C - 0.5 * C4_B - C4_A)/4 * (T_r / 1000)^4 + ... 
             (C5_D + C5_C - 0.5 * C5_B - C5_A)/-1 * (T_r / 1000)^-1)) * 1000;

     EntR2 = (((-110.52 - 241.83) - (-115.90 + 0.5 * 0)) * 1000 + ...
            ((C1_E + C1_D - 0.5 * C1_B - C1_C)/1 + (T_r / 1000)^1 + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (T_r / 1000)^2 + ... 
             (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (T_r / 1000)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (T_r / 1000)^4 + ...
             (C5_E + C5_D - 0.5 * C5_B - C5_C)/-1 * (T_r / 1000)^-1) - ...
            ((C1_E + C1_D - 0.5 * C1_B - C1_C)/1 + (T_r / 1000)^1 + (C2_E + C2_D - 0.5 * C2_B - C2_C)/2 * (T_r / 1000)^2 + ... 
             (C3_E + C3_D - 0.5 * C3_B - C3_C)/3 * (T_r / 1000)^3 + (C4_E + C4_D - 0.5 * C4_B - C4_C)/4 * (T_r / 1000)^4 + ...
             (C5_E + C5_D - 0.5 * C5_B - C5_C)/-1 * (T_r / 1000)^-1)) * 1000;

% Aspectos de calor y transferencia de masa 
     Pr_R = Mu_R * Cp_R / CondTerm_R;
     h_CalorPelI = (0.4 * Re_R^0.5 + 0.2 * Re_R ^ (2/3)) * Pr_R^0.4 * (1 - Vacio) / Vacio * CondTerm_R / dp;
     h_CalorPelE = 0.2 * ((do * NR_G) / Mu_G)^0.6 * (Cp_G * Mu_G)^0.33 * CondTerm_R / do;
     h_CalorR = 0.458 / Vacio * ((dp * NR_R / Mu_R) ^ -0.407) * (Cp_R * FR / PM_R) * ((Cp_R * Mu_R) / (CondTerm_R * PM_R)) ^ (-2/3);
     U = 1 / (1 / h_CalorPelI * (d / CondTerm_Tubo) + (di / MLDd) * (1 / h_CalorPelE) + Censu); % [W / m2 s]

% 5.C RESOLUCIÓN DE ECUACIONES ALGEBRÁICAS

% Resolución de expresión algebráica:
     n1 = 1; n2 = 1; 
     syms Balance(T_s)
     Balance(T_s) = h_CalorR * ap * Ro_p * (T_s - T_r) - ...
           n1 * EntR1 * (k1 * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) / ...
           (1 + K1A * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) - ...
           n2 * EntR2 * ((k2) * (PR_in * y_Cin * (Q_Cin / Q_C) * (FR_C / FR_Cin) * (y(3)/T_rin))) / ...
           (1 + K2D * (PR_in * y_Din * (Q_Din / Q_D) * (FR_D / FR_Din) * (y(3)/T_rin)));
     T_s = vpasolve(Balance);
    
% Paquete de Ecuaciones: El Diferencial se encuentra despejado para ser solucionado:
% Variables: t: x, y(1): FR_A, y(2): FR_E, y(3) : Tr, y(4) : Tc, y(5): Ts, y(6) : P 
     Ec = zeros(5,1);
          Tubos = 294;
     % Ecuación 1: Balance de Materia de Metanol (A) - Se encuentra dFA / dx igualado a la expresión:
     Ec1 = -(k1 * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) / ...
             (1 + K1A * (PR_in * y_Ain * (Q_Ain / Q_A) * (y(1) / FR_Ain) * (y(3)/T_rin))) * Ro_R * A_transi * Tubos;

     % Ecuación 2: Balance de Materia de Monóxido de Carbono (E) - Se encuentra dFE / dx igualado a la expresión:
     Ec2 = (k2 * (PR_in * y_Cin * (Q_Cin / Q_C) * (FR_C / FR_Cin) * (y(3)/T_rin))) / ...
             (1 + K2D * (PR_in * y_Din * (Q_Din / Q_D) * (FR_D / FR_Din) * (y(3)/T_rin))) * Ro_R * A_transi * Tubos;

     % Ecuación 3: Balance de Energía de Fluido de Reacción - Se encuentra dTr / dx igualado a la expresión:
     Ec3 = -(h_calor * ap * Ro_Lecho * (y(3) - y(5)) * A_transi + U * pi() * do * (y(4) - y(3))) / ...
             (FR_A * Cp_A + FR_B * Cp_B + FR_C * Cp_C + FR_D * Cp_D + FR_E * Cp_E + FR_F * Cp_F);

     % Ecuación 4: Balance de Energía de Fluido de Control: - Se encuentra dTc / dx igualado a la expresión:
     Ec4 = -(U * pi() * do * (y(4) - y(3))) / (f_G * Cp_G);

     % Ecuación 5: Balance de Presión - Se encuentra dPR / dz igualado a la expresión
     Ec5 = -f * Ro_R * Vel_s / dp;

     % Conjunto de Ecuaciones:
     Ec = [Ec1; Ec2; Ec3; Ec4; Ec5];

end

% Fin :)
