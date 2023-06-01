% Datos y Parámetros:
     PM543 = [0.0793143 0.0766408 0.0614908 0.0583717 0.0993656 0.0966921 0.0891172 0.0917907 0.0877804 0.1082770];
     PF543 = [0.0445586 0.0454498 0.0704026 0.0663923 0.0289631 0.0280719 0.0392116 0.0409939 0.0499056 0.0151499];
     PW543 = [0.0463409 0.0472321 0.0744128 0.0708481 0.0302998 0.0289631 0.0405483 0.0423307 0.0516880 0.0151499];
     r1543 = [8.9062E-06 8.6458E-06 6.7708E-06 6.6146E-06 1.1198E-05 1.0677E-05 1.0156E-05 1.0208E-05 9.8958E-06 1.1875E-05];
     r2543 = [1.5625E-07 1.5625E-07 1.9792E-07 2.0833E-07 1.9271E-07 1.7187E-07 1.6667E-07 1.4062E-07 1.7708E-07 1.0937E-07];
     
     PM602 = [0.0988002 0.0163020 0.0345801 0.0350741 0.0237121 0.0251941 0.0568101 0.0829922];
     PF602 = [0.0464361 0.1185600 0.1081860 0.1244880 0.0953422 0.0968242 0.0943542 0.0597741];
     PW602 = [0.0474241 0.1289340 0.1185600 0.1373320 0.1052220 0.1067040 0.1012700 0.0627382];
     r1602 = [3.1406E-05 7.1875E-06 1.2708E-05 1.3125E-05 9.8958E-06 9.8958E-06 2.1406E-05 2.7396E-05];
     r2602 = [3.9583E-07 5.2604E-07 5.2604E-07 6.4583E-07 4.8958E-07 4.6875E-07 7.6042E-07 6.0417E-07];

     % x1: k1 (Constante de reacción 1)
     % x2: K1M (Constante de Equilibrio para M)
     % x3: K1W (Constante de Equilibrio para W)
     % x4: k2 (Constante de reacción 2)
     % x5: K2F (Constante de Equilibrio para F)
     % x6: K2W (Constante de Equilibrio para W)

% Configuración previa para ciclos de iteración:
     ExpTol = -12; 
     Tol = 10 ^ - abs(ExpTol);
     nmax = 10000;
     h = 0.00000001;
     [S1, S2, S3] = deal(1, 1, 1);

Temper = input("Temperatura a revisar (543 o 602): ");

switch(Temper)
     case (543)
% Etapa 1: Alternativa A a 543 K (10 Ecuaciones):
     % Configuración de Vectores:
     PM = PM543; PF = PF543; PW = PW543; r1 = r1543; r2 = r2543;

     % Funciones para Cinética 1 (A - 543 K - 10 ecuaciones):
     r1n1 = @(x1, x2, x3) (x1^2*PM(1)) / (1 + x2^2*PM(1) + x3^2*PW(1)) - r1(1);
     r1n2 = @(x1, x2, x3) (x1^2*PM(2)) / (1 + x2^2*PM(2) + x3^2*PW(2)) - r1(2);
     r1n3 = @(x1, x2, x3) (x1^2*PM(3)) / (1 + x2^2*PM(3) + x3^2*PW(3)) - r1(3);
     r1n4 = @(x1, x2, x3) (x1^2*PM(4)) / (1 + x2^2*PM(4) + x3^2*PW(4)) - r1(4);
     r1n5 = @(x1, x2, x3) (x1^2*PM(5)) / (1 + x2^2*PM(5) + x3^2*PW(5)) - r1(5);
     r1n6 = @(x1, x2, x3) (x1^2*PM(6)) / (1 + x2^2*PM(6) + x3^2*PW(6)) - r1(6);
     r1n7 = @(x1, x2, x3) (x1^2*PM(7)) / (1 + x2^2*PM(7) + x3^2*PW(7)) - r1(7);
     r1n8 = @(x1, x2, x3) (x1^2*PM(8)) / (1 + x2^2*PM(8) + x3^2*PW(8)) - r1(8);
     r1n9 = @(x1, x2, x3) (x1^2*PM(9)) / (1 + x2^2*PM(9) + x3^2*PW(9)) - r1(9);
     r1n10 = @(x1, x2, x3) (x1^2*PM(10)) / (1 + x2^2*PM(10) + x3^2*PW(10)) - r1(10);
     
     % Derivadas para Cinética 1 (A - 543 K - 10 ecuaciones):
     dr1n1 = @(x1, hx1, x2, hx2, x3, hx3) (r1n1(x1 + hx1, x2 + hx2, x3 + hx3) - r1n1(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n2 = @(x1, hx1, x2, hx2, x3, hx3) (r1n2(x1 + hx1, x2 + hx2, x3 + hx3) - r1n2(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n3 = @(x1, hx1, x2, hx2, x3, hx3) (r1n3(x1 + hx1, x2 + hx2, x3 + hx3) - r1n3(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n4 = @(x1, hx1, x2, hx2, x3, hx3) (r1n4(x1 + hx1, x2 + hx2, x3 + hx3) - r1n4(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n5 = @(x1, hx1, x2, hx2, x3, hx3) (r1n5(x1 + hx1, x2 + hx2, x3 + hx3) - r1n5(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n6 = @(x1, hx1, x2, hx2, x3, hx3) (r1n6(x1 + hx1, x2 + hx2, x3 + hx3) - r1n6(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n7 = @(x1, hx1, x2, hx2, x3, hx3) (r1n7(x1 + hx1, x2 + hx2, x3 + hx3) - r1n7(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n8 = @(x1, hx1, x2, hx2, x3, hx3) (r1n8(x1 + hx1, x2 + hx2, x3 + hx3) - r1n8(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n9 = @(x1, hx1, x2, hx2, x3, hx3) (r1n9(x1 + hx1, x2 + hx2, x3 + hx3) - r1n9(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n10 = @(x1, hx1, x2, hx2, x3, hx3) (r1n10(x1 + hx1, x2 + hx2, x3 + hx3) - r1n10(x1, x2, x3)) / (hx1 + hx2 + hx3);
    
     % Funciones para Cinética 2 (A - 543 K - 10 ecuaciones):
     r2n1 = @(x4, x5, x6) (x4^2*PF(1)) / (1 + x5^2*PF(1) + x6^2*PW(1)) - r2(1);
     r2n2 = @(x4, x5, x6) (x4^2*PF(2)) / (1 + x5^2*PF(2) + x6^2*PW(2)) - r2(2);
     r2n3 = @(x4, x5, x6) (x4^2*PF(3)) / (1 + x5^2*PF(3) + x6^2*PW(3)) - r2(3);
     r2n4 = @(x4, x5, x6) (x4^2*PF(4)) / (1 + x5^2*PF(4) + x6^2*PW(4)) - r2(4);
     r2n5 = @(x4, x5, x6) (x4^2*PF(5)) / (1 + x5^2*PF(5) + x6^2*PW(5)) - r2(5);
     r2n6 = @(x4, x5, x6) (x4^2*PF(6)) / (1 + x5^2*PF(6) + x6^2*PW(6)) - r2(6);
     r2n7 = @(x4, x5, x6) (x4^2*PF(7)) / (1 + x5^2*PF(7) + x6^2*PW(7)) - r2(7);
     r2n8 = @(x4, x5, x6) (x4^2*PF(8)) / (1 + x5^2*PF(8) + x6^2*PW(8)) - r2(8);
     r2n9 = @(x4, x5, x6) (x4^2*PF(9)) / (1 + x5^2*PF(9) + x6^2*PW(9)) - r2(9);
     r2n10 = @(x4, x5, x6) (x4^2*PF(10)) / (1 + x5^2*PF(10) + x6^2*PW(10)) - r2(10);
     
     % Derivadas para Cinética 2 (A - 543 K - 10 ecuaciones):
     dr2n1 = @(x4, hx4, x5, hx5, x6, hx6) (r2n1(x4 + hx4, x5 + hx5, x6 + hx6) - r2n1(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n2 = @(x4, hx4, x5, hx5, x6, hx6) (r2n2(x4 + hx4, x5 + hx5, x6 + hx6) - r2n2(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n3 = @(x4, hx4, x5, hx5, x6, hx6) (r2n3(x4 + hx4, x5 + hx5, x6 + hx6) - r2n3(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n4 = @(x4, hx4, x5, hx5, x6, hx6) (r2n4(x4 + hx4, x5 + hx5, x6 + hx6) - r2n4(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n5 = @(x4, hx4, x5, hx5, x6, hx6) (r2n5(x4 + hx4, x5 + hx5, x6 + hx6) - r2n5(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n6 = @(x4, hx4, x5, hx5, x6, hx6) (r2n6(x4 + hx4, x5 + hx5, x6 + hx6) - r2n6(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n7 = @(x4, hx4, x5, hx5, x6, hx6) (r2n7(x4 + hx4, x5 + hx5, x6 + hx6) - r2n7(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n8 = @(x4, hx4, x5, hx5, x6, hx6) (r2n8(x4 + hx4, x5 + hx5, x6 + hx6) - r2n8(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n9 = @(x4, hx4, x5, hx5, x6, hx6) (r2n9(x4 + hx4, x5 + hx5, x6 + hx6) - r2n9(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n10 = @(x4, hx4, x5, hx5, x6, hx6) (r2n10(x4 + hx4, x5 + hx5, x6 + hx6) - r2n10(x4, x5, x6)) / (hx4 + hx5 + hx6);

     % Establecimiento del ciclo de iteración para Cinética 1:
     VTol1 = [Tol, Tol, Tol];
     [R1, R2, R3] = deal(S1, S2, S3); 
     RaicesA1 = [R1 R2 R3];
     RaicesB1 = [R1 R2 R3]; 
     Error1 = [1, 1, 1];
     n = 0;

     while max(Error1) > max(VTol1)
          % Matriz de derivadas parciales evaluadas:
          Jacob1 = [dr1n1(R1, h, R2, 0, R3, 0), dr1n1(R1, 0, R2, h, R3, 0), dr1n1(R1, 0, R2, 0, R3, h);
                    dr1n2(R1, h, R2, 0, R3, 0), dr1n2(R1, 0, R2, h, R3, 0), dr1n2(R1, 0, R2, 0, R3, h);
                    dr1n3(R1, h, R2, 0, R3, 0), dr1n3(R1, 0, R2, h, R3, 0), dr1n3(R1, 0, R2, 0, R3, h);
                    dr1n4(R1, h, R2, 0, R3, 0), dr1n4(R1, 0, R2, h, R3, 0), dr1n4(R1, 0, R2, 0, R3, h);
                    dr1n5(R1, h, R2, 0, R3, 0), dr1n5(R1, 0, R2, h, R3, 0), dr1n5(R1, 0, R2, 0, R3, h);
                    dr1n6(R1, h, R2, 0, R3, 0), dr1n6(R1, 0, R2, h, R3, 0), dr1n6(R1, 0, R2, 0, R3, h);
                    dr1n7(R1, h, R2, 0, R3, 0), dr1n7(R1, 0, R2, h, R3, 0), dr1n7(R1, 0, R2, 0, R3, h);
                    dr1n8(R1, h, R2, 0, R3, 0), dr1n8(R1, 0, R2, h, R3, 0), dr1n8(R1, 0, R2, 0, R3, h);
                    dr1n9(R1, h, R2, 0, R3, 0), dr1n9(R1, 0, R2, h, R3, 0), dr1n9(R1, 0, R2, 0, R3, h);
                    dr1n10(R1, h, R2, 0, R3, 0), dr1n10(R1, 0, R2, h, R3, 0), dr1n10(R1, 0, R2, 0, R3, h)];
        % Vector de funciones evaluadas:
             Vec1 = [r1n1(R1, R2, R3) r1n2(R1, R2, R3) r1n3(R1, R2, R3) r1n4(R1, R2, R3) r1n5(R1, R2, R3), ...
                     r1n6(R1, R2, R3) r1n7(R1, R2, R3) r1n8(R1, R2, R3) r1n9(R1, R2, R3) r1n10(R1, R2, R3)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta1 = Jacob1 \ -Vec1;
        % X(i + 1) = X(i) + deltaX:
             RaicesB1 = RaicesA1 + RaicesDelta1';
        % Evaluación del Error y continuación:
             Error1 = abs(RaicesB1 - RaicesA1);
             RaicesA1 = RaicesB1;
             [R1, R2, R3] = deal(RaicesA1(1), RaicesA1(2), RaicesA1(3));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error1 = [0, 0, 0];
             end
     end 
     x1A543 = R1^2;
     x2A543 = R2^2;
     x3A543 = R3^2;

     % Establecimiento del ciclo de iteración para Cinética 2:
     Jacob1 = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
     VTol1 = [Tol, Tol, Tol];
     [R1, R2, R3] = deal(S1, S2, S3); 
     RaicesA1 = [R1 R2 R3];
     RaicesB1 = [R1 R2 R3]; 
     Error1 = [1, 1, 1];
     n = 0;

     while max(Error1) > max(VTol1)
        % Matriz de derivadas parciales evaluadas:
             Jacob1 = [dr2n1(R1, h, R2, 0, R3, 0), dr2n1(R1, 0, R2, h, R3, 0), dr2n1(R1, 0, R2, 0, R3, h);
                       dr2n2(R1, h, R2, 0, R3, 0), dr2n2(R1, 0, R2, h, R3, 0), dr2n2(R1, 0, R2, 0, R3, h);
                       dr2n3(R1, h, R2, 0, R3, 0), dr2n3(R1, 0, R2, h, R3, 0), dr2n3(R1, 0, R2, 0, R3, h);
                       dr2n4(R1, h, R2, 0, R3, 0), dr2n4(R1, 0, R2, h, R3, 0), dr2n4(R1, 0, R2, 0, R3, h);
                       dr2n5(R1, h, R2, 0, R3, 0), dr2n5(R1, 0, R2, h, R3, 0), dr2n5(R1, 0, R2, 0, R3, h);
                       dr2n6(R1, h, R2, 0, R3, 0), dr2n6(R1, 0, R2, h, R3, 0), dr2n6(R1, 0, R2, 0, R3, h);
                       dr2n7(R1, h, R2, 0, R3, 0), dr2n7(R1, 0, R2, h, R3, 0), dr2n7(R1, 0, R2, 0, R3, h);
                       dr2n8(R1, h, R2, 0, R3, 0), dr2n8(R1, 0, R2, h, R3, 0), dr2n8(R1, 0, R2, 0, R3, h);
                       dr2n9(R1, h, R2, 0, R3, 0), dr2n9(R1, 0, R2, h, R3, 0), dr2n9(R1, 0, R2, 0, R3, h);
                       dr2n10(R1, h, R2, 0, R3, 0), dr2n10(R1, 0, R2, h, R3, 0), dr2n10(R1, 0, R2, 0, R3, h)];
        % Vector de funciones evaluadas:
             Vec1 = [r2n1(R1, R2, R3) r2n2(R1, R2, R3) r2n3(R1, R2, R3) r2n4(R1, R2, R3) r2n5(R1, R2, R3), ...
                     r2n6(R1, R2, R3) r2n7(R1, R2, R3) r2n8(R1, R2, R3) r2n9(R1, R2, R3) r2n10(R1, R2, R3)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta1 = Jacob1 \ -Vec1;
        % X(i + 1) = X(i) + deltaX:
             RaicesB1 = RaicesA1 + RaicesDelta1';
        % Evaluación del Error y continuación:
             Error1 = abs(RaicesB1 - RaicesA1);
             RaicesA1 = RaicesB1;
             [R1, R2, R3] = deal(RaicesA1(1), RaicesA1(2), RaicesA1(3));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error1 = [0, 0, 0];
             end
     end
     x4A543 = R1^2;
     x5A543 = R2^2;
     x6A543 = R3^2;

% Etapa 2: Alternativa B a 543 K (10 Ecuaciones):         

     % Funciones para Cinética 1 (B - 543 K - 10 ecuaciones):
     r1n1 = @(x1, x2) (x1^2*PM(1)) / (1 + x2^2*PM(1)) - r1(1);
     r1n2 = @(x1, x2) (x1^2*PM(2)) / (1 + x2^2*PM(2)) - r1(2);
     r1n3 = @(x1, x2) (x1^2*PM(3)) / (1 + x2^2*PM(3)) - r1(3);
     r1n4 = @(x1, x2) (x1^2*PM(4)) / (1 + x2^2*PM(4)) - r1(4);
     r1n5 = @(x1, x2) (x1^2*PM(5)) / (1 + x2^2*PM(5)) - r1(5);
     r1n6 = @(x1, x2) (x1^2*PM(6)) / (1 + x2^2*PM(6)) - r1(6);
     r1n7 = @(x1, x2) (x1^2*PM(7)) / (1 + x2^2*PM(7)) - r1(7);
     r1n8 = @(x1, x2) (x1^2*PM(8)) / (1 + x2^2*PM(8)) - r1(8);
     r1n9 = @(x1, x2) (x1^2*PM(9)) / (1 + x2^2*PM(9)) - r1(9);
     r1n10 = @(x1, x2) (x1^2*PM(10)) / (1 + x2^2*PM(10)) - r1(10);
     
     % Derivadas para Cinética 1 (B - 543 K - 10 ecuaciones):
     dr1n1 = @(x1, hx1, x2, hx2) (r1n1(x1 + hx1, x2 + hx2) - r1n1(x1, x2)) / (hx1 + hx2);
     dr1n2 = @(x1, hx1, x2, hx2) (r1n2(x1 + hx1, x2 + hx2) - r1n2(x1, x2)) / (hx1 + hx2);
     dr1n3 = @(x1, hx1, x2, hx2) (r1n3(x1 + hx1, x2 + hx2) - r1n3(x1, x2)) / (hx1 + hx2);
     dr1n4 = @(x1, hx1, x2, hx2) (r1n4(x1 + hx1, x2 + hx2) - r1n4(x1, x2)) / (hx1 + hx2);
     dr1n5 = @(x1, hx1, x2, hx2) (r1n5(x1 + hx1, x2 + hx2) - r1n5(x1, x2)) / (hx1 + hx2);
     dr1n6 = @(x1, hx1, x2, hx2) (r1n6(x1 + hx1, x2 + hx2) - r1n6(x1, x2)) / (hx1 + hx2);
     dr1n7 = @(x1, hx1, x2, hx2) (r1n7(x1 + hx1, x2 + hx2) - r1n7(x1, x2)) / (hx1 + hx2);
     dr1n8 = @(x1, hx1, x2, hx2) (r1n8(x1 + hx1, x2 + hx2) - r1n8(x1, x2)) / (hx1 + hx2);
     dr1n9 = @(x1, hx1, x2, hx2) (r1n9(x1 + hx1, x2 + hx2) - r1n9(x1, x2)) / (hx1 + hx2);
     dr1n10 = @(x1, hx1, x2, hx2) (r1n10(x1 + hx1, x2 + hx2) - r1n10(x1, x2)) / (hx1 + hx2);
     
     % Funciones para Cinética 2 (B - 543 K - 10 ecuaciones):
     r2n1 = @(x4, x6) (x4^2*PF(1)) / (1 + x6^2*PW(1)) - r2(1);
     r2n2 = @(x4, x6) (x4^2*PF(2)) / (1 + x6^2*PW(2)) - r2(2);
     r2n3 = @(x4, x6) (x4^2*PF(3)) / (1 + x6^2*PW(3)) - r2(3);
     r2n4 = @(x4, x6) (x4^2*PF(4)) / (1 + x6^2*PW(4)) - r2(4);
     r2n5 = @(x4, x6) (x4^2*PF(5)) / (1 + x6^2*PW(5)) - r2(5);
     r2n6 = @(x4, x6) (x4^2*PF(6)) / (1 + x6^2*PW(6)) - r2(6);
     r2n7 = @(x4, x6) (x4^2*PF(7)) / (1 + x6^2*PW(7)) - r2(7);
     r2n8 = @(x4, x6) (x4^2*PF(8)) / (1 + x6^2*PW(8)) - r2(8);
     r2n9 = @(x4, x6) (x4^2*PF(9)) / (1 + x6^2*PW(9)) - r2(9);
     r2n10 = @(x4, x6) (x4^2*PF(10)) / (1 + x6^2*PW(10)) - r2(10);
     
     % Funciones para Cinética 2 (B - 543 K - 10 ecuaciones):
     dr2n1 = @(x4, hx4, x6, hx6) (r2n1(x4 + hx4, x6 + hx6) - r2n1(x4, x6)) / (hx4 + hx6);
     dr2n2 = @(x4, hx4, x6, hx6) (r2n2(x4 + hx4, x6 + hx6) - r2n2(x4, x6)) / (hx4 + hx6);
     dr2n3 = @(x4, hx4, x6, hx6) (r2n3(x4 + hx4, x6 + hx6) - r2n3(x4, x6)) / (hx4 + hx6);
     dr2n4 = @(x4, hx4, x6, hx6) (r2n4(x4 + hx4, x6 + hx6) - r2n4(x4, x6)) / (hx4 + hx6);
     dr2n5 = @(x4, hx4, x6, hx6) (r2n5(x4 + hx4, x6 + hx6) - r2n5(x4, x6)) / (hx4 + hx6);
     dr2n6 = @(x4, hx4, x6, hx6) (r2n6(x4 + hx4, x6 + hx6) - r2n6(x4, x6)) / (hx4 + hx6);
     dr2n7 = @(x4, hx4, x6, hx6) (r2n7(x4 + hx4, x6 + hx6) - r2n7(x4, x6)) / (hx4 + hx6);
     dr2n8 = @(x4, hx4, x6, hx6) (r2n8(x4 + hx4, x6 + hx6) - r2n8(x4, x6)) / (hx4 + hx6);
     dr2n9 = @(x4, hx4, x6, hx6) (r2n9(x4 + hx4, x6 + hx6) - r2n9(x4, x6)) / (hx4 + hx6);
     dr2n10 = @(x4, hx4, x6, hx6) (r2n10(x4 + hx4, x6 + hx6) - r2n10(x4, x6)) / (hx4 + hx6);

     % Establecimiento del ciclo de iteración para Cinética 1:
     VTol2 = [Tol, Tol];
     [R1, R2] = deal(S1, S2); 
     RaicesA2 = [R1 R2];
     RaicesB2 = [R1 R2]; 
     Error2 = [1, 1];
     n = 0;

     while max(Error2) > max(VTol2)
          % Matriz de derivadas parciales evaluadas:
          Jacob2 = [dr1n1(R1, h, R2, 0), dr1n1(R1, 0, R2, h);
                    dr1n2(R1, h, R2, 0), dr1n2(R1, 0, R2, h);
                    dr1n3(R1, h, R2, 0), dr1n3(R1, 0, R2, h);
                    dr1n4(R1, h, R2, 0), dr1n4(R1, 0, R2, h);
                    dr1n5(R1, h, R2, 0), dr1n5(R1, 0, R2, h);
                    dr1n6(R1, h, R2, 0), dr1n6(R1, 0, R2, h);
                    dr1n7(R1, h, R2, 0), dr1n7(R1, 0, R2, h);
                    dr1n8(R1, h, R2, 0), dr1n8(R1, 0, R2, h);
                    dr1n9(R1, h, R2, 0), dr1n9(R1, 0, R2, h);
                    dr1n10(R1, h, R2, 0), dr1n10(R1, 0, R2, h)];
        % Vector de funciones evaluadas:
             Vec2 = [r1n1(R1, R2) r1n2(R1, R2) r1n3(R1, R2) r1n4(R1, R2) r1n5(R1, R2), ...
                     r1n6(R1, R2) r1n7(R1, R2) r1n8(R1, R2) r1n9(R1, R2) r1n10(R1, R2)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta2 = Jacob2 \ -Vec2;
        % X(i + 1) = X(i) + deltaX:
             RaicesB2 = RaicesA2 + RaicesDelta2';
        % Evaluación del Error y continuación:
             Error2 = abs(RaicesB2 - RaicesA2);
             RaicesA2 = RaicesB2;
             [R1, R2] = deal(RaicesA2(1), RaicesA2(2));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error2 = [0, 0];
             end
     end 
     x1B543 = R1^2;
     x2B543 = R2^2;

     % Establecimiento del ciclo de iteración para Cinética 2:
     Jacob2 = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
     VTol2 = [Tol, Tol];
     [R1, R2] = deal(S1, S2); 
     RaicesA2 = [R1 R2];
     RaicesB2 = [R1 R2]; 
     Error2 = [1, 1];
     n = 0;

     while max(Error2) > max(VTol2)
        % Matriz de derivadas parciales evaluadas:
             Jacob2 = [dr2n1(R1, h, R2, 0), dr2n1(R1, 0, R2, h);
                       dr2n2(R1, h, R2, 0), dr2n2(R1, 0, R2, h);
                       dr2n3(R1, h, R2, 0), dr2n3(R1, 0, R2, h);
                       dr2n4(R1, h, R2, 0), dr2n4(R1, 0, R2, h);
                       dr2n5(R1, h, R2, 0), dr2n5(R1, 0, R2, h);
                       dr2n6(R1, h, R2, 0), dr2n6(R1, 0, R2, h);
                       dr2n7(R1, h, R2, 0), dr2n7(R1, 0, R2, h);
                       dr2n8(R1, h, R2, 0), dr2n8(R1, 0, R2, h);
                       dr2n9(R1, h, R2, 0), dr2n9(R1, 0, R2, h);
                       dr2n10(R1, h, R2, 0), dr2n10(R1, 0, R2, h)];
        % Vector de funciones evaluadas:
             Vec2 = [r2n1(R1, R2) r2n2(R1, R2) r2n3(R1, R2) r2n4(R1, R2) r2n5(R1, R2), ...
                     r2n6(R1, R2) r2n7(R1, R2) r2n8(R1, R2) r2n9(R1, R2) r2n10(R1, R2)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta2 = Jacob2 \ -Vec2;
        % X(i + 1) = X(i) + deltaX:
             RaicesB2 = RaicesA2 + RaicesDelta2';
        % Evaluación del Error y continuación:
             Error2 = abs(RaicesB2 - RaicesA2);
             RaicesA2 = RaicesB2;
             [R1, R2] = deal(RaicesA2(1), RaicesA2(2));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error2 = [0, 0];
             end
     end
     x4B543 = R1^2;
     x6B543 = R2^2;

% Respuestas:
     disp("A 543 K se tienen los siguientes parámetros por alternativas (A o B) bajo una tolerancia de " + Tol + " :")
     disp(" - k1 (A): " + x1A543);
     disp(" - K1M (A): " + x2A543);
     disp(" - K1W (A): " + x3A543);
     disp(" - k2 (A): " + x4A543);
     disp(" - K2F (A): " + x5A543);
     disp(" - K2W (A): " + x6A543);
     disp(" - k1 (B): " + x1B543);
     disp(" - K1M (B): " + x2B543);
     disp(" - k2 (B): " + x4B543);
     disp(" - K2W (B): " + x6B543);

case (602)
% Etapa 3: Alternativa A a 602 K (8 Ecuaciones):
     % Configuración de Vectores:
     PM = PM602; PF = PF602; PW = PW602; r1 = r1602; r2 = r2602;

     % Funciones para Cinética 1 (A - 543 K - 10 ecuaciones):
     r1n1 = @(x1, x2, x3) (x1^2*PM(1)) / (1 + x2^2*PM(1) + x3^2*PW(1)) - r1(1);
     r1n2 = @(x1, x2, x3) (x1^2*PM(2)) / (1 + x2^2*PM(2) + x3^2*PW(2)) - r1(2);
     r1n3 = @(x1, x2, x3) (x1^2*PM(3)) / (1 + x2^2*PM(3) + x3^2*PW(3)) - r1(3);
     r1n4 = @(x1, x2, x3) (x1^2*PM(4)) / (1 + x2^2*PM(4) + x3^2*PW(4)) - r1(4);
     r1n5 = @(x1, x2, x3) (x1^2*PM(5)) / (1 + x2^2*PM(5) + x3^2*PW(5)) - r1(5);
     r1n6 = @(x1, x2, x3) (x1^2*PM(6)) / (1 + x2^2*PM(6) + x3^2*PW(6)) - r1(6);
     r1n7 = @(x1, x2, x3) (x1^2*PM(7)) / (1 + x2^2*PM(7) + x3^2*PW(7)) - r1(7);
     r1n8 = @(x1, x2, x3) (x1^2*PM(8)) / (1 + x2^2*PM(8) + x3^2*PW(8)) - r1(8);
     
     % Derivadas para Cinética 1 (A - 543 K - 10 ecuaciones):
     dr1n1 = @(x1, hx1, x2, hx2, x3, hx3) (r1n1(x1 + hx1, x2 + hx2, x3 + hx3) - r1n1(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n2 = @(x1, hx1, x2, hx2, x3, hx3) (r1n2(x1 + hx1, x2 + hx2, x3 + hx3) - r1n2(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n3 = @(x1, hx1, x2, hx2, x3, hx3) (r1n3(x1 + hx1, x2 + hx2, x3 + hx3) - r1n3(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n4 = @(x1, hx1, x2, hx2, x3, hx3) (r1n4(x1 + hx1, x2 + hx2, x3 + hx3) - r1n4(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n5 = @(x1, hx1, x2, hx2, x3, hx3) (r1n5(x1 + hx1, x2 + hx2, x3 + hx3) - r1n5(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n6 = @(x1, hx1, x2, hx2, x3, hx3) (r1n6(x1 + hx1, x2 + hx2, x3 + hx3) - r1n6(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n7 = @(x1, hx1, x2, hx2, x3, hx3) (r1n7(x1 + hx1, x2 + hx2, x3 + hx3) - r1n7(x1, x2, x3)) / (hx1 + hx2 + hx3);
     dr1n8 = @(x1, hx1, x2, hx2, x3, hx3) (r1n8(x1 + hx1, x2 + hx2, x3 + hx3) - r1n8(x1, x2, x3)) / (hx1 + hx2 + hx3);
    
     % Funciones para Cinética 2 (A - 543 K - 10 ecuaciones):
     r2n1 = @(x4, x5, x6) (x4^2*PF(1)) / (1 + x5^2*PF(1) + x6^2*PW(1)) - r2(1);
     r2n2 = @(x4, x5, x6) (x4^2*PF(2)) / (1 + x5^2*PF(2) + x6^2*PW(2)) - r2(2);
     r2n3 = @(x4, x5, x6) (x4^2*PF(3)) / (1 + x5^2*PF(3) + x6^2*PW(3)) - r2(3);
     r2n4 = @(x4, x5, x6) (x4^2*PF(4)) / (1 + x5^2*PF(4) + x6^2*PW(4)) - r2(4);
     r2n5 = @(x4, x5, x6) (x4^2*PF(5)) / (1 + x5^2*PF(5) + x6^2*PW(5)) - r2(5);
     r2n6 = @(x4, x5, x6) (x4^2*PF(6)) / (1 + x5^2*PF(6) + x6^2*PW(6)) - r2(6);
     r2n7 = @(x4, x5, x6) (x4^2*PF(7)) / (1 + x5^2*PF(7) + x6^2*PW(7)) - r2(7);
     r2n8 = @(x4, x5, x6) (x4^2*PF(8)) / (1 + x5^2*PF(8) + x6^2*PW(8)) - r2(8);
     
     % Derivadas para Cinética 2 (A - 543 K - 10 ecuaciones):
     dr2n1 = @(x4, hx4, x5, hx5, x6, hx6) (r2n1(x4 + hx4, x5 + hx5, x6 + hx6) - r2n1(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n2 = @(x4, hx4, x5, hx5, x6, hx6) (r2n2(x4 + hx4, x5 + hx5, x6 + hx6) - r2n2(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n3 = @(x4, hx4, x5, hx5, x6, hx6) (r2n3(x4 + hx4, x5 + hx5, x6 + hx6) - r2n3(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n4 = @(x4, hx4, x5, hx5, x6, hx6) (r2n4(x4 + hx4, x5 + hx5, x6 + hx6) - r2n4(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n5 = @(x4, hx4, x5, hx5, x6, hx6) (r2n5(x4 + hx4, x5 + hx5, x6 + hx6) - r2n5(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n6 = @(x4, hx4, x5, hx5, x6, hx6) (r2n6(x4 + hx4, x5 + hx5, x6 + hx6) - r2n6(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n7 = @(x4, hx4, x5, hx5, x6, hx6) (r2n7(x4 + hx4, x5 + hx5, x6 + hx6) - r2n7(x4, x5, x6)) / (hx4 + hx5 + hx6);
     dr2n8 = @(x4, hx4, x5, hx5, x6, hx6) (r2n8(x4 + hx4, x5 + hx5, x6 + hx6) - r2n8(x4, x5, x6)) / (hx4 + hx5 + hx6);

     % Establecimiento del ciclo de iteración para Cinética 1:
     VTol3 = [Tol, Tol, Tol];
     [R1, R2, R3] = deal(S1, S2, S3); 
     RaicesA3 = [R1 R2 R3];
     RaicesB3 = [R1 R2 R3]; 
     Error3 = [1, 1, 1];
     n = 0;

     while max(Error3) > max(VTol3)
          % Matriz de derivadas parciales evaluadas:
          Jacob3 = [dr1n1(R1, h, R2, 0, R3, 0), dr1n1(R1, 0, R2, h, R3, 0), dr1n1(R1, 0, R2, 0, R3, h);
                    dr1n2(R1, h, R2, 0, R3, 0), dr1n2(R1, 0, R2, h, R3, 0), dr1n2(R1, 0, R2, 0, R3, h);
                    dr1n3(R1, h, R2, 0, R3, 0), dr1n3(R1, 0, R2, h, R3, 0), dr1n3(R1, 0, R2, 0, R3, h);
                    dr1n4(R1, h, R2, 0, R3, 0), dr1n4(R1, 0, R2, h, R3, 0), dr1n4(R1, 0, R2, 0, R3, h);
                    dr1n5(R1, h, R2, 0, R3, 0), dr1n5(R1, 0, R2, h, R3, 0), dr1n5(R1, 0, R2, 0, R3, h);
                    dr1n6(R1, h, R2, 0, R3, 0), dr1n6(R1, 0, R2, h, R3, 0), dr1n6(R1, 0, R2, 0, R3, h);
                    dr1n7(R1, h, R2, 0, R3, 0), dr1n7(R1, 0, R2, h, R3, 0), dr1n7(R1, 0, R2, 0, R3, h);
                    dr1n8(R1, h, R2, 0, R3, 0), dr1n8(R1, 0, R2, h, R3, 0), dr1n8(R1, 0, R2, 0, R3, h)];
        % Vector de funciones evaluadas:
             Vec3 = [r1n1(R1, R2, R3) r1n2(R1, R2, R3) r1n3(R1, R2, R3) r1n4(R1, R2, R3) r1n5(R1, R2, R3), ...
                     r1n6(R1, R2, R3) r1n7(R1, R2, R3) r1n8(R1, R2, R3)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta3 = Jacob3 \ -Vec3;
        % X(i + 1) = X(i) + deltaX:
             RaicesB3 = RaicesA3 + RaicesDelta3';
        % Evaluación del Error y continuación:
             Error3 = abs(RaicesB3 - RaicesA3);
             RaicesA3 = RaicesB3;
             [R1, R2, R3] = deal(RaicesA3(1), RaicesA3(2), RaicesA3(3));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error3 = [0, 0, 0];
             end
     end 
     x1A602 = R1^2;
     x2A602 = R2^2;
     x3A602 = R3^2;

     % Establecimiento del ciclo de iteración para Cinética 2:
     Jacob3 = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
     VTol3 = [Tol, Tol, Tol];
     [R1, R2, R3] = deal(S1, S2, S3); 
     RaicesA3 = [R1 R2 R3];
     RaicesB3 = [R1 R2 R3]; 
     Error3 = [1, 1, 1];
     n = 0;

     while max(Error3) > max(VTol3)
        % Matriz de derivadas parciales evaluadas:
             Jacob3 = [dr2n1(R1, h, R2, 0, R3, 0), dr2n1(R1, 0, R2, h, R3, 0), dr2n1(R1, 0, R2, 0, R3, h);
                       dr2n2(R1, h, R2, 0, R3, 0), dr2n2(R1, 0, R2, h, R3, 0), dr2n2(R1, 0, R2, 0, R3, h);
                       dr2n3(R1, h, R2, 0, R3, 0), dr2n3(R1, 0, R2, h, R3, 0), dr2n3(R1, 0, R2, 0, R3, h);
                       dr2n4(R1, h, R2, 0, R3, 0), dr2n4(R1, 0, R2, h, R3, 0), dr2n4(R1, 0, R2, 0, R3, h);
                       dr2n5(R1, h, R2, 0, R3, 0), dr2n5(R1, 0, R2, h, R3, 0), dr2n5(R1, 0, R2, 0, R3, h);
                       dr2n6(R1, h, R2, 0, R3, 0), dr2n6(R1, 0, R2, h, R3, 0), dr2n6(R1, 0, R2, 0, R3, h);
                       dr2n7(R1, h, R2, 0, R3, 0), dr2n7(R1, 0, R2, h, R3, 0), dr2n7(R1, 0, R2, 0, R3, h);
                       dr2n8(R1, h, R2, 0, R3, 0), dr2n8(R1, 0, R2, h, R3, 0), dr2n8(R1, 0, R2, 0, R3, h)];
        % Vector de funciones evaluadas:
             Vec3 = [r2n1(R1, R2, R3) r2n2(R1, R2, R3) r2n3(R1, R2, R3) r2n4(R1, R2, R3) r2n5(R1, R2, R3), ...
                     r2n6(R1, R2, R3) r2n7(R1, R2, R3) r2n8(R1, R2, R3)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta3 = Jacob3 \ -Vec3;
        % X(i + 1) = X(i) + deltaX:
             RaicesB3 = RaicesA3 + RaicesDelta3';
        % Evaluación del Error y continuación:
             Error3 = abs(RaicesB3 - RaicesA3);
             RaicesA3 = RaicesB3;
             [R1, R2, R3] = deal(RaicesA3(1), RaicesA3(2), RaicesA3(3));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error3 = [0, 0, 0];
             end
     end
     x4A602 = R1^2;
     x5A602 = R2^2;
     x6A602 = R3^2;

% Etapa 4: Alternativa B a 602 K (8 Ecuaciones):         

     % Funciones para Cinética 1 (B - 543 K - 10 ecuaciones):
     r1n1 = @(x1, x2) (x1^2*PM(1)) / (1 + x2^2*PM(1)) - r1(1);
     r1n2 = @(x1, x2) (x1^2*PM(2)) / (1 + x2^2*PM(2)) - r1(2);
     r1n3 = @(x1, x2) (x1^2*PM(3)) / (1 + x2^2*PM(3)) - r1(3);
     r1n4 = @(x1, x2) (x1^2*PM(4)) / (1 + x2^2*PM(4)) - r1(4);
     r1n5 = @(x1, x2) (x1^2*PM(5)) / (1 + x2^2*PM(5)) - r1(5);
     r1n6 = @(x1, x2) (x1^2*PM(6)) / (1 + x2^2*PM(6)) - r1(6);
     r1n7 = @(x1, x2) (x1^2*PM(7)) / (1 + x2^2*PM(7)) - r1(7);
     r1n8 = @(x1, x2) (x1^2*PM(8)) / (1 + x2^2*PM(8)) - r1(8);
     
     % Derivadas para Cinética 1 (B - 543 K - 10 ecuaciones):
     dr1n1 = @(x1, hx1, x2, hx2) (r1n1(x1 + hx1, x2 + hx2) - r1n1(x1, x2)) / (hx1 + hx2);
     dr1n2 = @(x1, hx1, x2, hx2) (r1n2(x1 + hx1, x2 + hx2) - r1n2(x1, x2)) / (hx1 + hx2);
     dr1n3 = @(x1, hx1, x2, hx2) (r1n3(x1 + hx1, x2 + hx2) - r1n3(x1, x2)) / (hx1 + hx2);
     dr1n4 = @(x1, hx1, x2, hx2) (r1n4(x1 + hx1, x2 + hx2) - r1n4(x1, x2)) / (hx1 + hx2);
     dr1n5 = @(x1, hx1, x2, hx2) (r1n5(x1 + hx1, x2 + hx2) - r1n5(x1, x2)) / (hx1 + hx2);
     dr1n6 = @(x1, hx1, x2, hx2) (r1n6(x1 + hx1, x2 + hx2) - r1n6(x1, x2)) / (hx1 + hx2);
     dr1n7 = @(x1, hx1, x2, hx2) (r1n7(x1 + hx1, x2 + hx2) - r1n7(x1, x2)) / (hx1 + hx2);
     dr1n8 = @(x1, hx1, x2, hx2) (r1n8(x1 + hx1, x2 + hx2) - r1n8(x1, x2)) / (hx1 + hx2);
     
     % Funciones para Cinética 2 (B - 543 K - 10 ecuaciones):
     r2n1 = @(x4, x6) (x4^2*PF(1)) / (1 + x6^2*PW(1)) - r2(1);
     r2n2 = @(x4, x6) (x4^2*PF(2)) / (1 + x6^2*PW(2)) - r2(2);
     r2n3 = @(x4, x6) (x4^2*PF(3)) / (1 + x6^2*PW(3)) - r2(3);
     r2n4 = @(x4, x6) (x4^2*PF(4)) / (1 + x6^2*PW(4)) - r2(4);
     r2n5 = @(x4, x6) (x4^2*PF(5)) / (1 + x6^2*PW(5)) - r2(5);
     r2n6 = @(x4, x6) (x4^2*PF(6)) / (1 + x6^2*PW(6)) - r2(6);
     r2n7 = @(x4, x6) (x4^2*PF(7)) / (1 + x6^2*PW(7)) - r2(7);
     r2n8 = @(x4, x6) (x4^2*PF(8)) / (1 + x6^2*PW(8)) - r2(8);
     
     % Funciones para Cinética 2 (B - 543 K - 10 ecuaciones):
     dr2n1 = @(x4, hx4, x6, hx6) (r2n1(x4 + hx4, x6 + hx6) - r2n1(x4, x6)) / (hx4 + hx6);
     dr2n2 = @(x4, hx4, x6, hx6) (r2n2(x4 + hx4, x6 + hx6) - r2n2(x4, x6)) / (hx4 + hx6);
     dr2n3 = @(x4, hx4, x6, hx6) (r2n3(x4 + hx4, x6 + hx6) - r2n3(x4, x6)) / (hx4 + hx6);
     dr2n4 = @(x4, hx4, x6, hx6) (r2n4(x4 + hx4, x6 + hx6) - r2n4(x4, x6)) / (hx4 + hx6);
     dr2n5 = @(x4, hx4, x6, hx6) (r2n5(x4 + hx4, x6 + hx6) - r2n5(x4, x6)) / (hx4 + hx6);
     dr2n6 = @(x4, hx4, x6, hx6) (r2n6(x4 + hx4, x6 + hx6) - r2n6(x4, x6)) / (hx4 + hx6);
     dr2n7 = @(x4, hx4, x6, hx6) (r2n7(x4 + hx4, x6 + hx6) - r2n7(x4, x6)) / (hx4 + hx6);
     dr2n8 = @(x4, hx4, x6, hx6) (r2n8(x4 + hx4, x6 + hx6) - r2n8(x4, x6)) / (hx4 + hx6);

     % Establecimiento del ciclo de iteración para Cinética 1:
     VTol4 = [Tol, Tol];
     [R1, R2] = deal(S1, S2); 
     RaicesA4 = [R1 R2];
     RaicesB4 = [R1 R2]; 
     Error4 = [1, 1];
     n = 0;

     while max(Error4) > max(VTol4)
          % Matriz de derivadas parciales evaluadas:
          Jacob4 = [dr1n1(R1, h, R2, 0), dr1n1(R1, 0, R2, h);
                    dr1n2(R1, h, R2, 0), dr1n2(R1, 0, R2, h);
                    dr1n3(R1, h, R2, 0), dr1n3(R1, 0, R2, h);
                    dr1n4(R1, h, R2, 0), dr1n4(R1, 0, R2, h);
                    dr1n5(R1, h, R2, 0), dr1n5(R1, 0, R2, h);
                    dr1n6(R1, h, R2, 0), dr1n6(R1, 0, R2, h);
                    dr1n7(R1, h, R2, 0), dr1n7(R1, 0, R2, h);
                    dr1n8(R1, h, R2, 0), dr1n8(R1, 0, R2, h)];
        % Vector de funciones evaluadas:
             Vec4 = [r1n1(R1, R2) r1n2(R1, R2) r1n3(R1, R2) r1n4(R1, R2) r1n5(R1, R2), ...
                     r1n6(R1, R2) r1n7(R1, R2) r1n8(R1, R2)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta4 = Jacob4 \ -Vec4;
        % X(i + 1) = X(i) + deltaX:
             RaicesB4 = RaicesA4 + RaicesDelta4';
        % Evaluación del Error y continuación:
             Error4 = abs(RaicesB4 - RaicesA4);
             RaicesA4 = RaicesB4;
             [R1, R2] = deal(RaicesA4(1), RaicesA4(2));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error4 = [0, 0];
             end
     end 
     x1B602 = R1^2;
     x2B602 = R2^2;

     % Establecimiento del ciclo de iteración para Cinética 2:
     Jacob4 = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0;];
     VTol4 = [Tol, Tol];
     [R1, R2] = deal(S1, S2); 
     RaicesA4 = [R1 R2];
     RaicesB4 = [R1 R2]; 
     Error4 = [1, 1];
     n = 0;

     while max(Error4) > max(VTol4)
        % Matriz de derivadas parciales evaluadas:
             Jacob4 = [dr2n1(R1, h, R2, 0), dr2n1(R1, 0, R2, h);
                       dr2n2(R1, h, R2, 0), dr2n2(R1, 0, R2, h);
                       dr2n3(R1, h, R2, 0), dr2n3(R1, 0, R2, h);
                       dr2n4(R1, h, R2, 0), dr2n4(R1, 0, R2, h);
                       dr2n5(R1, h, R2, 0), dr2n5(R1, 0, R2, h);
                       dr2n6(R1, h, R2, 0), dr2n6(R1, 0, R2, h);
                       dr2n7(R1, h, R2, 0), dr2n7(R1, 0, R2, h);
                       dr2n8(R1, h, R2, 0), dr2n8(R1, 0, R2, h)];
        % Vector de funciones evaluadas:
             Vec4 = [r2n1(R1, R2) r2n2(R1, R2) r2n3(R1, R2) r2n4(R1, R2) r2n5(R1, R2), ...
                     r2n6(R1, R2) r2n7(R1, R2) r2n8(R1, R2)]';
        % Resolución del sistema de ecuaciones:
             RaicesDelta4 = Jacob4 \ -Vec4;
        % X(i + 1) = X(i) + deltaX:
             RaicesB4 = RaicesA4 + RaicesDelta4';
        % Evaluación del Error y continuación:
             Error4 = abs(RaicesB4 - RaicesA4);
             RaicesA4 = RaicesB4;
             [R1, R2] = deal(RaicesA4(1), RaicesA4(2));
             n = n + 1;
        % En caso de emergencia:
             if n == nmax
                  Error4 = [0, 0];
             end
     end
     x4B602 = R1^2;
     x6B602 = R2^2;

     % Respuestas:
          disp("A 602 K se tienen los siguientes parámetros por alternativas (A o B) bajo una tolerancia de " + Tol + " :")
          disp(" - k1 (A): " + x1A602);
          disp(" - K1M (A): " + x2A602);
          disp(" - K1W (A): " + x3A602);
          disp(" - k2 (A): " + x4A602);
          disp(" - K2F (A): " + x5A602);
          disp(" - K2W (A): " + x6A602);
          disp(" - k1 (B): " + x1B602);
          disp(" - K1M (B): " + x2B602);
          disp(" - k2 (B): " + x4B602);
          disp(" - K2W (B): " + x6B602);
end



