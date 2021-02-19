%% Parameters
clc;
clear;
parameters;

%% Plant
% I'll just examinate case 1 (arm 2 position at 180), case 2 is trivial
A=A1;
B=B1;
Plant=ss(A,B,C,D)

%% Controllability test
% Note that i'm using sym
rank(sym(ctrb(A,B)))

%% Observability test
% Note that i'm using sym
rank(sym(obsv(A,C)))

%% Full state feedback
% 1. Desired poles, three cases: bessel poles, ITAE poles, dominant pole
w0 = 5;

% 1.1 Bessel poles
pb1 = -0.9264;
pb2 = -0.5906+0.9072i;
pb3 = -0.5906-0.9072i;
pb4 = -0.8516+0.4427i;
pb5 = -0.8516-0.4427i;
Pbessel = [pb1 pb2 pb3 pb4 pb5]*w0

% 1.2 ITAE poles
pi1 = -0.8955;
pi2 = -0.3764+1.2920i;
pi3 = -0.3764-1.2920i;
pi4 = -0.5758+0.5339i;
pi5 = -0.5758-0.5339i;
Pitae = [pi1 pi2 pi3 pi4 pi5]*w0

% 1.3 Dominant poles
mp=0.04/100; % Overshoot
ts=0.6; % Settling time
z=abs(log(mp))/sqrt(pi^2+log(mp)^2); % Damping ratio
wn=1/(z*ts);  % Natural frequency
alpha=20; % ¿?
Pdominant = roots(conv([1 2*z*wn wn^2], poly([-alpha*wn+2i -alpha*wn-2i -alpha*wn])))

% 2. Gain
Kbessel = place(A1,B1,Pbessel);
Kitae = place(A1,B1,Pitae);
Kdominant = place(A1,B1,Pdominant);

% 3. Pre filter
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
t=0:0.01:15;
r=[zeros(size(t)); (pi*180/180)*ones(size(t))];

[Ybessel,~,~]=lsim(Pc_bessel, r, t, CI);
[Yitae,~,~]=lsim(Pc_itae, r, t, CI);
[Ydominant,T,X]=lsim(Pc_dominant, r, t, CI);

% 6. Plot
% To degrees and add operation point cause linearization
Ybessel=(Ybessel+r')*180/pi;
Yitae=(Yitae+r')*180/pi;;
Ydominant=(Ydominant+r')*180/pi;
r=(r')*180/pi
T = [T T];

figure(1)
plot(T,Ybessel,'r',T,Yitae,'g',T,Ydominant,'b',T,r,'--k')
legend('X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant','X1-reference','X2-reference')
