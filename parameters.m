% Parámetros del sistema.
L1=0.278;
L2=0.300;
l1=0.150;
l2=0.148;
m1=0.300;
m2=0.075;
J1=2.48e-2;
J2=3.86e-3;
b1=1.00e-4;
b2=2.80e-4;
L=0.005;
R=7.80;
Km=0.090;
g=9.81;

% Parámetros de las consideraciones
J_hat_1=J1+m1*l1^2;
J_hat_2=J2+m2*l2^2;
J_hat_0=J_hat_1+m2*L1^2;

% Parámetros de espacio de estados
Aijdenominador=J_hat_0*J_hat_2-m2^2*L1^2*l2^2;
A31=0;
A32=g*m2^2*l2^2*L1/Aijdenominador;
A33=-b1*J_hat_2/Aijdenominador;
A34=-b2*m2*l2*L1/Aijdenominador;
A41=0;
A42=g*m2*l2*J_hat_0/Aijdenominador;
A43=-b1*m2*l2*L1/Aijdenominador;
A44=-b2*J_hat_0/Aijdenominador;
B31=J_hat_2/Aijdenominador;
B32=m2*L1*l2/Aijdenominador;
B41=m2*L1*l2/Aijdenominador;
B42=J_hat_0/Aijdenominador;

% Modelo linealizado, 180 grados
A1=[0 0 1 0 0; 0 0 0 1 0; A31 A32 A33 A34 B31*Km; A41 A42 A43 A44 B41*Km; 0 0 -Km/L 0 -R/L];
B1=[0 0; 0 0; 0 B32; 0 B42; 1/L 0]; % Incluye tau2
%B1=[0; 0; 0 ; 0 ; 1/L];

% Modelo linealizado, 0 grados
A2=[0 0 1 0 0; 0 0 0 1 0; A31 A32 A33 -A34 B31*Km; A41 -A42 -A43 A44 -B41*Km; 0 0 -Km/L 0 -R/L];
B2=[0 0; 0 0; 0 -B32; 0 B42; 1/L 0]; % Incluye Tau2
%B2=[0; 0; 0; 0; 1/L];

% C
C=[1 0 0 0 0; 0 1 0 0 0];

% n = orden de la planta
% p = número de entradas
% q = número de salidas
n=length(A1);
p=width(B1);
q=height(C);

% D
D=zeros(q,p);

% Condiciones iniciales, caso 1
%CI1=[0 pi 0 0 0];
CI1=[0 pi*180/180 0 0 0];

% Condiciones iniciales, caso 2
%CI2=[0 0 0 0 0];
CI2=[0 pi*5/180 0 0 0];

% Condiciones iniciales para dinámica error (L), caso 1
CIL1=[0 pi*170/180 0 0 0];

% Condiciones iniciales para dinámica error (L), caso 2
CIL2=[0 pi*10/180 0 0 0];