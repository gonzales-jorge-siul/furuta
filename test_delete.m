
%% Problema1 _ Lab5
%[A,B2,C,D]=linmod('motor_lab5')
%B=B2(:,1);Bw=B2(:,2);
%% condiciones de dise�o
mp=0.9/100;
ts=3;
z=abs(log(mp))/sqrt(pi^2+log(mp)^2);
wn=1/(z*ts); %% f�rmula del 2%
alpha=4; %% Modelo de referencia de orden superior a 2
p=roots(conv([1 2*z*wn wn^2], [ 1 alpha*wn]));
k=acker(A,B,p);
N=1/dcgain(A-B*k,B,C,0) ;  % 1/(C*inv(B*K-A)*B)
step(A-B*k,N*B,C,0),
%% Observador de estados
%mp=0.9/100;
ts=1; % mas r�pido que el controlador
%z=abs(log(mp))/sqrt(pi^2+log(mp)^2);
wn=4/(z*ts);
alpha=4; %% Modelo de referencia de orden superior a 2
po=roots(conv([1 2*z*wn wn^2], [ 1 alpha*wn]));
L=acker(A',C',po)';% Ganancia de observaci�n
% LSR
% Q=C'*C;
% k= lqr(A,B,C'*C,1);
% 
% L=lqr(A',C',C'*C,1)'

%% Efecto integral
Aa=[A zeros(3,1); -C 0];
Ba=[B;0];
pa=[p ;-100];

Ka=place(Aa,Ba,pa);
%Ka=lqr(Aa,Ba,[Q zeros(3,1); 0 0 0  1],1);
K=Ka(1:3);
Ki=-Ka(end);

