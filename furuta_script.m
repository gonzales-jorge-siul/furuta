%% Parameters
clc;
clear;
parameters;
fig_count = 0;
PLOT_CONTROLLER=true;
PLOT_OBSERVER=false;
PLOT_OBERVER_CONTROLLER=true;
PLOT_REDUCED_OBERVER_CONTROLLER=true;


%% Plant
% I'll just examinate case 1 (arm 2 position at 180 degrees), case 2 is trivial
A=A1;
B=B1;
Plant=ss(A,B,C,D);

%% Controllability test
% Note that i'm using sym
fprintf('Expected controllability matrix rank: %i\n',n)
fprintf('System controllability matrix rank: %i\n',rank(sym(ctrb(A,B))))

%% Observability test
% Note that i'm using sym
fprintf('Expected observability matrix rank: %i\n',n)
fprintf('System observability matrix rank: %i\n',rank(sym(obsv(A,C))))

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
DPc_dominant = roots(conv([1 2*z*wn wn^2], poly([-alpha*z*wn+2i -alpha*z*wn-2i -alpha*z*wn])));

% 2. 'K'
K_bessel = place(A1,B1,DPc_bessel);
K_itae = place(A1,B1,DPc_itae);
K_dominant = place(A1,B1,DPc_dominant);

% 3. Pre-filter
Nb_bessel=-(C*(A-B*K_bessel)^-1*B)\eye(q,q);
Nb_itae=-(C*(A-B*K_itae)^-1*B)\eye(q,q);
Nb_dominant=-(C*(A-B*K_dominant)^-1*B)\eye(q,q);

% 4. Plant with controller
Dc=zeros(q,q);
Pc_bessel = ss(A-B*K_bessel,B*Nb_bessel,C,Dc);
Pc_itae = ss(A-B*K_itae,B*Nb_itae,C,Dc);
Pc_dominant = ss(A-B*K_dominant,B*Nb_dominant,C,Dc);

CI=[-5*pi/180 -pi*5/180 0 0 0]; % Note that CI are different cause linearization

if PLOT_CONTROLLER 
   % 5. Simulate
    t=0:0.01:5;
    r=[zeros(size(t)); (pi*180/180)*ones(size(t))];

    [Ybessel,~,X]=lsim(Pc_bessel, r, t, CI);
    U_bessel = Nb_bessel*r - K_bessel*X';
    [Yitae,~,X]=lsim(Pc_itae, r, t, CI);
    U_itae = Nb_itae*r - K_itae*X';
    [Ydominant,T,X]=lsim(Pc_dominant, r, t, CI);
    U_dominant = Nb_dominant*r - K_dominant*X';

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
end

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

if PLOT_OBSERVER
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
end

%% Controller + Observer, K+L
% K, L previously computed remain same
% N previously computed remain same if N=BM
% 1. Close loop plant
% 1.1 Bessel
Acl_bessel=[(A-B*K_bessel) B*K_bessel; zeros(n,n) (A-L_bessel*C)];
Bcl_bessel=[B*Nb_bessel; zeros(n,q)];
Ccl_bessel=[C zeros(q,n)];
Dcl_bessel=zeros(q,q);

Poc_bessel=ss(Acl_bessel,Bcl_bessel,Ccl_bessel,Dcl_bessel);

% 1.2 Itae
Acl_itae=[(A-B*K_itae) B*K_itae; zeros(n,n) (A-L_itae*C)];
Bcl_itae=[B*Nb_itae; zeros(n,q)];
Ccl_itae=[C zeros(q,n)];
Dcl_itae=zeros(q,q);

Poc_itae=ss(Acl_itae,Bcl_itae,Ccl_itae,Dcl_itae);

% 1.3 Dominant
Acl_dominant=[(A-B*K_dominant) B*K_dominant; zeros(n,n) (A-L_dominant*C)];
Bcl_dominant=[B*Nb_dominant; zeros(n,q)];
Ccl_dominant=[C zeros(q,n)];
Dcl_dominant=zeros(q,q);

Poc_dominant=ss(Acl_dominant,Bcl_dominant,Ccl_dominant,Dcl_dominant);

if PLOT_OBERVER_CONTROLLER
    % 2. Simulate
    CI=[-5*pi/180 -pi*5/180 0 0 0 zeros(1,n)]; % Note that CI are different cause linearization
    t=0:0.01:5;
    r=[zeros(size(t)); (pi*180/180)*ones(size(t))];

    [Yoc_bessel,~,X]=lsim(Poc_bessel, r, t, CI);
    X=X';
    Uoc_bessel=Nb_bessel*r-K_bessel*(X(1:5,:)-X(6:10,:));

    [Yoc_itae,~,X]=lsim(Poc_itae, r, t, CI);
    X=X';
    Uoc_itae=Nb_itae*r-K_itae*(X(1:5,:)-X(6:10,:));

    [Yoc_dominant,T,X]=lsim(Poc_dominant, r, t, CI);
    X=X';
    Uoc_dominant = Nb_dominant*r - K_dominant*(X(1:5,:)-X(6:10,:));

    % 3. Plot
    % To degrees and add operation point to output cause linearization
    Yoc_bessel=(Yoc_bessel+r')*180/pi;
    Yoc_itae=(Yoc_itae+r')*180/pi;
    Yoc_dominant=(Yoc_dominant+r')*180/pi;
    r=(r')*180/pi;
    T = [T T];

    fig_count=fig_count+1;
    figure(fig_count)

    subplot(2,1,1)
    plot(T,Yoc_bessel,'r',T,Yoc_itae,'g',T,Yoc_dominant,'b',T,r,'--k')
    legend('X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant','X1-reference','X2-reference')
    title('Output')
    xlabel('Time(s)')
    ylabel('Output(degrees)')

    T=T(:,1);
    subplot(2,1,2)
    plot(T,Uoc_bessel,'r',T,Uoc_itae,'g',T,Uoc_dominant,'b')
    legend('U-Bessel','U-ITAE','U-Dominant')
    title('Control signal')
    xlabel('Time(s)')
    ylabel('Input(V)')
end

%% Reduced-order observer
% Measured states: x1, x2 => p=2
% Estimated state: x3, x4, x5 => r=3
% |A11 A12||X1|-> Estimated
% |A21 A22||X2|-> Measured
% | 0   I ||X1|
%          |X2|
% 1. Transform the system to get measured states to the right.
p=2;
r=n-p;
V=[0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0];  % Arbitrary linear independet matrix
P=[V; C]^-1; % Transform matrix
Ab=P\A*P;
Bb=P\B;
Cb=C*P;

% 2. Separate matrix into components
A11=Ab(1:r,1:r);
A12=Ab(1:r,r+1:n);
A21=Ab(r+1:n,1:r);
A22=Ab(r+1:n,r+1:n);
B1=Bb(1:r);
B2=Bb(r+1:n);

Ar=A11;
Cr=A21;

% 3. Observability test for reduced system
fprintf('Expected observability matrix rank: %i\n',r)
fprintf('Reduced system observability matrix rank: %i\n',rank(obsv(Ar,Cr)))

% 4. Desired poles, three cases: bessel poles, ITAE poles, dominant pole
w0 = 5;

% 4.1 Bessel desired poles
pb1 = -0.7081;
pb2 = -0.5210+1.068i;
pb3 = -0.5210-1.068i;
DPr_bessel = [pb1 pb2 pb3]*w0;

% 4.2 ITAE desired poles
pi1 = -0.9420;
pi2 = -0.7455+0.7112i;
pi3 = -0.7455-0.7112i;
DPr_itae = [pi1 pi2 pi3]*w0;

% 4.3 Desired dominant poles
mp=0.04/100; % Overshoot
ts=0.6; % Settling time
z=abs(log(mp))/sqrt(pi^2+log(mp)^2); % Damping ratio
wn=1/(z*ts);  % Natural frequency
alpha=20; % Factor that indicates how far are non-dominant poles from dominant ones
DPr_dominant = roots(conv([1 2*z*wn wn^2], poly(-alpha*z*wn)));

% 5. L computation
Lr_bessel=place(Ar',Cr',DPr_bessel)';
Lr_itae=place(Ar',Cr',DPr_itae)';
Lr_dominant=place(Ar',Cr',DPr_dominant)';

%F_bessel=A11-Lr_bessel*A21;
%G_bessel=B1-Lr_bessel*B2;
%H_bessel=(A11-Lr_bessel*A21)*Lr_bessel+A12-Lr_bessel*A22;
%F_itae=A11-Lr_itae*A21;
%G_itae=B1-Lr_itae*B2;
%H_itae=(A11-Lr_itae*A21)*Lr_itae+A12-Lr_itae*A22;
%F_dominant=A11-Lr_dominant*A21;
%G_dominant=B1-Lr_dominant*B2;
%H_dominant=(A11-Lr_dominant*A21)*Lr_dominant+A12-Lr_dominant*A22;

% 6. Close loop plant
% 6.1 Bessel
helpMatrix = zeros(n,n);
tempA22=(Ar-Lr_bessel*Cr);
for i=1:height(tempA22), for j=1:width(tempA22), helpMatrix(i,j)=tempA22(i,j); end; end
Arcl_bessel=[(Ab-Bb*K_bessel*P) Bb*K_bessel*P; zeros(n,n) helpMatrix];
Brcl_bessel=[Bb*Nb_bessel; zeros(n,q)];
Crcl_bessel=[Cb zeros(q,n)];
Drcl_bessel=zeros(q,q);

Proc_bessel=ss(Arcl_bessel,Brcl_bessel,Crcl_bessel,Drcl_bessel);

% 6.2 Itae
helpMatrix = zeros(n,n);
tempA22=(Ar-Lr_itae*Cr);
for i=1:height(tempA22), for j=1:width(tempA22), helpMatrix(i,j)=tempA22(i,j); end; end
Arcl_itae=[(Ab-Bb*K_itae*P) Bb*K_itae*P; zeros(n,n) helpMatrix];
Brcl_itae=[Bb*Nb_itae; zeros(n,q)];
Crcl_itae=[Cb zeros(q,n)];
Drcl_itae=zeros(q,q);

Proc_itae=ss(Arcl_itae,Brcl_itae,Crcl_itae,Drcl_itae);

% 6.3 Dominant
helpMatrix = zeros(n,n);
tempA22=(Ar-Lr_dominant*Cr);
for i=1:height(tempA22), for j=1:width(tempA22), helpMatrix(i,j)=tempA22(i,j); end; end
Arcl_dominant=[(Ab-Bb*K_dominant*P) Bb*K_dominant*P; zeros(n,n) helpMatrix];
Brcl_dominant=[Bb*Nb_dominant; zeros(n,q)];
Crcl_dominant=[Cb zeros(q,n)];
Drcl_dominant=zeros(q,q);

Proc_dominant=ss(Arcl_dominant,Brcl_dominant,Crcl_dominant,Drcl_dominant);

if PLOT_REDUCED_OBERVER_CONTROLLER
    % 2. Simulate
    CI=[-5*pi/180 -pi*5/180 0 0 0 zeros(1,n)]; % Note that CI are different cause linearization
    t=0:0.01:5;
    r=[zeros(size(t)); (pi*180/180)*ones(size(t))];

    [Yroc_bessel,~,X]=lsim(Proc_bessel, r, t, CI);
    X=X';
    Uroc_bessel=Nb_bessel*r-K_bessel*P*(X(1:5,:)-X(6:10,:));

    [Yroc_itae,~,X]=lsim(Proc_itae, r, t, CI);
    X=X';
    Uroc_itae=Nb_itae*r-K_itae*P*(X(1:5,:)-X(6:10,:));

    [Yroc_dominant,T,X]=lsim(Proc_dominant, r, t, CI);
    X=X';
    Uroc_dominant = Nb_dominant*r - K_dominant*P*(X(1:5,:)-X(6:10,:));

    % 3. Plot
    % To degrees and add operation point to output cause linearization
    Yroc_bessel=(Yroc_bessel+r')*180/pi;
    Yroc_itae=(Yroc_itae+r')*180/pi;
    Yroc_dominant=(Yroc_dominant+r')*180/pi;
    r=(r')*180/pi;
    T = [T T];

    fig_count=fig_count+1;
    figure(fig_count)

    subplot(2,1,1)
    plot(T,Yroc_bessel,'r',T,Yroc_itae,'g',T,Yroc_dominant,'b',T,r,'--k')
    legend('X1-Bessel','X2-Bessel','X1-ITAE','X2-ITAE','X1-Dominant','X2-Dominant','X1-reference','X2-reference')
    title('Output')
    xlabel('Time(s)')
    ylabel('Output(degrees)')

    T=T(:,1);
    subplot(2,1,2)
    plot(T,Uroc_bessel,'r',T,Uroc_itae,'g',T,Uroc_dominant,'b')
    legend('U-Bessel','U-ITAE','U-Dominant')
    title('Control signal')
    xlabel('Time(s)')
    ylabel('Input(V)')
end
