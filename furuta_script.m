%% Parameters
clc;
clear;
parameters;
fig_count = 0;

%% Plant
% I'll just examinate case 1 (arm 2 position at 180 degrees), case 2 is trivial
A=A1;
B=B1;
Plant=ss(A,B,C,D);

%% Controllability test
% Note that i'm using sym
rank(sym(ctrb(A,B)))

%% Observability test
% Note that i'm using sym
rank(sym(obsv(A,C)))

%% Full state feedback
% 1. Desired poles, three cases: bessel poles, ITAE poles, dominant pole
w0 = 5;

% 1.1 Bessel desired poles
pb1 = -0.9264;
pb2 = -0.5906+0.9072i;
pb3 = -0.5906-0.9072i;
pb4 = -0.8516+0.4427i;
pb5 = -0.8516-0.4427i;
DPc_bessel = [pb1 pb2 pb3 pb4 pb5]*w0;

% 1.2 ITAE desired poles
pi1 = -0.8955;
pi2 = -0.3764+1.2920i;
pi3 = -0.3764-1.2920i;
pi4 = -0.5758+0.5339i;
pi5 = -0.5758-0.5339i;
DPc_itae = [pi1 pi2 pi3 pi4 pi5]*w0;

% 1.3 Desired dominant poles
mp=0.04/100; % Overshoot
ts=0.6; % Settling time
z=abs(log(mp))/sqrt(pi^2+log(mp)^2); % Damping ratio
wn=1/(z*ts);  % Natural frequency
alpha=20; % ¿?
DPc_dominant = roots(conv([1 2*z*wn wn^2], poly([-alpha*wn+2i -alpha*wn-2i -alpha*wn])));

% 2. 'K'
Kbessel = place(A1,B1,DPc_bessel);
Kitae = place(A1,B1,DPc_itae);
Kdominant = place(A1,B1,DPc_dominant);

% 3. Pre-filter
Nbbessel=-(C*(A-B*Kbessel)^-1*B)\eye(q,q);
Nbitae=-(C*(A-B*Kitae)^-1*B)\eye(q,q);
Nbdominant=-(C*(A-B*Kdominant)^-1*B)\eye(q,q);

% 4. Plant with controller
Dc=zeros(q,q);
Pc_bessel = ss(A-B*Kbessel,B*Nbbessel,C,Dc);
Pc_itae = ss(A-B*Kitae,B*Nbitae,C,Dc);
Pc_dominant = ss(A-B*Kdominant,B*Nbdominant,C,Dc);

CI=[-5*pi/180 -pi*5/180 0 0 0]; % Note that CI are different cause linearization

% 5. Simulate
t=0:0.01:5;
r=[zeros(size(t)); (pi*180/180)*ones(size(t))];

[Ybessel,~,X]=lsim(Pc_bessel, r, t, CI);
U_bessel = Nbbessel*r - Kbessel*X';
[Yitae,~,X]=lsim(Pc_itae, r, t, CI);
U_itae = Nbitae*r - Kitae*X';
[Ydominant,T,X]=lsim(Pc_dominant, r, t, CI);
U_dominant = Nbdominant*r - Kdominant*X';

% 6. Plot
% To degrees and add operation point to output cause linearization
Ybessel=(Ybessel+r')*180/pi;
Yitae=(Yitae+r')*180/pi;
Ydominant=(Ydominant+r')*180/pi;
r=(r')*180/pi;
T = [T T];

fig_count=fig_count+1;
figure(fig_count)

subplot(2,1,1)
plot(T,Ybessel,'r',T,Yitae,'g',T,Ydominant,'b',T,r,'--k')
legend('X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant','X1-reference','X2-reference')
title('Output')
xlabel('Time(s)')
ylabel('Output(degrees)')

T=T(:,1);
subplot(2,1,2)
plot(T,U_bessel,'r',T,U_itae,'g',T,U_dominant,'b')
legend('U-Bessel','U-ITAE','U-Dominant')
title('Control signal')
xlabel('Time(s)')
ylabel('Input(V)')

%fig_count=fig_count+1;
%figure(fig_count)
%bodeplot(Plant,'k',Pc_bessel,'r',Pc_itae,'g',Pc_dominant,'b')

%% Observer
rap = 4; % Observer desired poles are 'rap' faster than controller

% 1. Desired poles
% 1.1 Bessel
DPo_bessel = real(DPc_bessel(:,:))*rap + imag(DPc_bessel(:,:))*1i;
% 1.2 ITAE
DPo_itae = real(DPc_itae(:,:))*rap + imag(DPc_itae(:,:))*1i;
% 1.3 Dominant
DPo_dominant = real(DPc_dominant(:,:))*rap + imag(DPc_dominant(:,:))*1i;

% 2. 'L'
L_bessel = place(A',C',DPo_bessel)';
L_itae = place(A',C',DPo_itae)';
L_dominant = place(A',C',DPo_dominant)';

% 3. Observer
D=zeros(q,p+q);
Po_bessel = ss(A-L_bessel*C,[B L_bessel],C,D);
Po_itae = ss(A-L_itae*C,[B L_itae],C,D);
Po_dominant = ss(A-L_dominant*C,[B L_dominant],C,D);

% 4. Plot
t=0:0.01:5;
u=t;
u(1)=1e6; % impulse
[Y,~,~]=lsim(Plant, u, t);

u=[u; Y'];
[Y_bessel,~,~]=lsim(Po_bessel, u, t);
[Y_itae,~,~]=lsim(Po_itae, u, t);
[Y_dominant,T,~]=lsim(Po_dominant, u, t);

fig_count=fig_count+1;
figure(fig_count)

subplot(2,1,1)
plot(T,Y,'k',T,Y_bessel,'r',T,Y_itae,'g',T,Y_dominant,'b')
legend('X1','X2','X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant')
title('Output')
xlabel('Time(s)')
ylabel('Output(rad)')

% Since scale is big, i plot error in porcentual terms
E_bessel = (Y-Y_bessel)*100./Y;
E_itae = (Y-Y_itae)*100./Y;
E_dominant = (Y-Y_dominant)*100./Y;

subplot(2,1,2)
plot(T,E_bessel,'r',T,E_itae,'g',T,E_dominant,'b')
legend('X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant')
title('Error')
xlabel('Time(s)')
ylabel('Error(%)')

%% Controller + Observer, K+L

