%%
%% %%%%%%%%%  Lower Limb Rehabilitation Robot - Decoupling Based  %%%%%%%%%%%%%
%%

%% Clean-Up
clear; close all; clc;

%% Change Default Font to Cambria
fontname = 'Cambria'                 ;
fontsize = 12                        ;
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);

%% Laplace Variable
s = tf('s');

%% Parameters
m1 = 4.0  ; m2 = 3    ;
M1 = 9.847; M2 = 3.086;
M3 = 1.93 ; 
g  = 9.8  ;
l1 = 0.496; l2 = 0.396;
L1 = 0.265; L2 = 0.229;
I1 = 0.164; I2 = 0.028;
b1 = 1e-1 ; b2 = 1e-1 ;

% Please fill out the ... lines!!!
a11 =((m1+M1)*(L1^2))+((m2+M2)*(l1^2))+((m2+M2)*(L2^2))+(2*(m2+M2)*(l1*L2))+(M3*((l1+l2)^2))+I1; 
a12 = -(((m2+M2)*((L2^2)+(l1*l2)))+(M3*((l2^2)+(l1+l2))));
a21 =a12;
a22 =((m2+M2)*(L2^2))+(M3*(l2^2))+I2;

J  = [a11 a12; a21 a22] ;



B  = [b1 0; 0 b2];

c11 = g*(((m1+M1)*L1+((m2+M2)*(l1+L2))+(M3*(l1+l2))));
c12 = (-g)*((L2*(m2+M2))+(M3*l2));
c21 = (-g)*((L2*(m2+M2)+(M3*l2)));
c22 = g*((L2*(m2+M2)+(M3*l2)));

C  = [c11 c12; c21 c22]    ;

O   = zeros(2);
I   = eye(2)  ;

Gam = eye(2);
Sgm = eye(2);

%% Plant
Ap = [O I; -inv(J)*C -inv(J)*B];
Bp = [O; inv(J)*Gam]           ;
Cp = [Sgm O]                   ;
Dp = O                         ;

% Here you can define the corresponding models: State-space and transfer
% function representations

P_tf  = ss(Ap,Bp,Cp,Dp)   ; zpk(P_tf);

P.StateName   = {'theta1','theta2','theta1_dot','theta2_dot'};
P.InputName   = {'u_1','u_2'};
P.OutputName  = {'y_1','y_2'};
 

%Check the bode plot of the original plant
figure
bode(P_tf)
%% Decoupled Plant
%  Define the decoupled plant (review the MIMO_Part2 lecture) here and obtain P1_tf and P2_tf
gamma = 1/(a11*a22+(a21^2));
phi_1 = gamma*a22*b1;
phi_2 = gamma*a11*b2;
delta_1=gamma*-a21*b2;
delta_2=gamma*-a21*b1;

P1_tf   = (1/(s*(s+phi_1)));
P2_tf   = (1/(s*(s+phi_2)));
figure
bode(P1_tf,P2_tf)
P1_ss= ss(P1_tf);
P2_ss=ss(P2_tf);
%% Controllers for Decoupled Plant
% Design the controller for P1
beta=1/max(norm(delta_1),norm(delta_2));
%Reference motion p1 part 1
% y1_0     = 0;
% y2_0     = 0;
% y1_dot_0 = 0;
% y2_dot_0 = 0;
% 
% t1i  = 0          ;
% t1m  = t1i + 2.5  ;
% t1f  = t1m + 2.0  ;
% 
% t2i  = 0          ;
% t2m  = t2i + 2.5  ;
% t2f  = t2m + 2.0  ;
% tf   = 6          ;
% dt   = 0.01       ;
% 
% r1i  = 0         ;
% r2i  = 0         ;
% r1m  = 4.5*pi/18 ;    
% r2m  = 7.0*pi/18 ; 
% r1f  = 4.0*pi/18 ;
% r2f  = 4.0*pi/18 ;

%section 13.3.3
omega_c1 = 20;
%lecture notes slide 20
alpha1 = 0.1;
zeta1 = 0.2;

ki1=omega_c1*(sqrt((omega_c1^2)+(phi_1^2)));
PM1=rad2deg(asin((1-alpha1)/(1+alpha1))-atan(zeta1)+atan(phi_1/omega_c1));

C1_tf=(ki1*((s+(omega_c1*sqrt(alpha1)))/((sqrt(alpha1)*s)+omega_c1))*(s+(zeta1*omega_c1))/(s*sqrt(1+(zeta1^2))));
% Check the Small Gain Condition
figure
bode(P1_tf/(1+(P1_tf*C1_tf)), (beta/s))
legend
%Check controller tuning
OL1_tf=C1_tf*P1_tf;
figure
margin(OL1_tf)

% Design the controller for P2
omega_c2 = 20;

alpha2 = 0.1;
zeta2 = 0.2;

ki2 = omega_c2*(sqrt((omega_c2^2)+(phi_2^2)));
PM2 = rad2deg(asin((1-alpha2)/(1+alpha2))-atan(zeta2)+atan(phi_2/omega_c2));

C2_tf = (ki2*((s+(omega_c2*sqrt(alpha2)))/(omega_c2+(s*sqrt(alpha2))))*((s+(zeta2*omega_c2))/(s*sqrt(1+(zeta2^2)))));

%L2_tf = ;

% Check the Small Gain Condition
figure
bode(P2_tf/(1+(P2_tf*C2_tf)), (beta/s))
legend

%Check controller tuning
OL2_tf=C2_tf*P2_tf;
figure
margin(OL2_tf)

% Define the decoupled controller here! 
Cd_tf =[(C1_tf) 0 ; 0  (C2_tf)];

%The following plotting commands might be helpful for freq. domain cont.
%design techiques. You can modify them accordingly.  

%figure(2); margin(L1_tf,'m--'); hold on; margin(L2_tf,'m-.'); grid on; 
%figure(9); subplot(221); nyqlog(tf(L1_tf));subplot(223); nyqlog(tf(L2_tf));
%subplot(122); nyquist(L1_tf,'m--',L2_tf,'m-.'); grid off;
%legend('Loop 1','Loop 2');
%axis([-15 0 -5 5]);


% SISO Loops
L1 = P1_tf*C1_tf ;
S1 = (1+L1)^(-1) ;
T1 = L1*S1 ;
N1 = S1*C1_tf;
R1 = S1*P1_tf;

L2 = P2_tf*C2_tf;
S2 = (1+L2)^(-1) ;
T2 = L2*S2 ;
N2 = S2*C2_tf;
R2 = S2*P2_tf;

%% Bode Diagrams of SISO Loops
opt = bodeoptions;
w_min = 1; w_max = 1000;
opt.YLimMode = 'manual';
opt.XlimMode = 'manual';
opt.XLim     = {[w_min w_max]};

opt.YLim = {[-50 10]};
figure(3);subplot(221); 
bodemag(S1,'m--',S2,'m-.',opt)   ; grid on;
title('Sensitivity')             ;
legend('Loop 1','Loop 2','Location','SouthEast');

opt.YLim = {[-50 10]};
figure(3); subplot(224); 
bodemag(T1,'m--',T2 ,'m-.',opt)   ; grid on;
title('Complementary Sensitivity');  

opt.YLim = {[-60 0]};
figure(3); subplot(222);
bodemag(beta/s,'r-',R1,'m--',R2,'m-.',opt); grid on;
title('Load Sensitivity')                ; 
legend('\beta/s','Location','Best');

opt.YLim = {[20 80]};
figure(3); subplot(223); 
bodemag(N1,'m--',N2,'m-.',opt)    ; grid on;
title('Noise Sensitivity')        ;



%% MIMO Controller

C_tf = inv(Gam)*(J*Cd_tf - (C))*inv(Sgm) ; 
%F_tf = inv(Gam)*ss(C)*inv(Sgm); 

[A1,B1,C1,D1] = ssdata(C_tf)              ;
%[A2,B2,C2,D2] = ssdata(F_tf)              ;


%% Check of Internal Stability of Closed-Loop
[A3,B3,C3,D3] = ssdata(P_tf)           ;
[A4,B4,C4,D4] = ssdata(C_tf)           ;

% .... continue to check internal stability!





%% MIMO Loop

L_tf = P_tf*C_tf;
S_tf = 1/(1+L_tf);
T_tf = L_tf*S_tf;
N_tf = C_tf*S_tf;
R_tf = S_tf*P_tf;
I = S_tf + T_tf;

L_in = C_tf*P_tf ;
S_in = inv(I+C_tf*P_tf); 
T_in = L_in*S_in;

T_er = S_tf;
T_yr = T_tf;
T_ur = N_tf;

%% Repeated Stability Check (You can check it again for safety!)


%% Bode Diagrams: Closed-Loop
opt = bodeoptions;
w_min = 0.01; w_max = 1000;
opt.YLimMode = 'manual';
opt.XlimMode = 'manual';
opt.XLim     = {[w_min w_max]};

opt.YLim = {[-50 10]};
%figure(4);subplot(221); 
figure(4); bodemag(T_er,'b--',S_tf,'m',opt) ; grid on;
title('Sensitivity')                                 ;
lg=legend('$$\mathcal{T}\_{er}$$')                   ;
set(lg,'Interpreter','latex','Location','SouthEast') ;

opt.YLim = {[-50 10]};
%figure(4); subplot(224); 
figure(5); bodemag(T_yr,'b--',T_tf ,'m',opt) ; grid on;
title('Complementary Sensitivity')                    ;  
lg=legend('$$\mathcal{T}\_{yr}$$')                    ;
set(lg,'Interpreter','latex','Location','SouthWest')  ;

opt.YLim = {[-60 0]};
%figure(4); subplot(222);
figure(6); bodemag(R_tf,'m',opt)             ; grid on;
title('Load Sensitivity')                    ; 

opt.YLim = {[20 80]};
%figure(4); subplot(223); 
figure(7); bodemag(T_ur,'b--',N_tf,'m',opt)  ; grid on;
title('Noise Sensitivity')                            ;
lg=legend('$$\mathcal{T}\_{ur}$$')                    ;
set(lg,'Interpreter','latex','Location','SouthEast')  ;


%evalfr(zpk(...),...);

%% Simulations
% Initial Conditions
kp = [400 0 ;0 600];
kd =[10 0;0 10];
y1_0     = 0;
y2_0     = 0;
y1_dot_0 = 0;
y2_dot_0 = 0;

t1i  = 0          ;
t1m  = t1i + 2.5  ;
t1f  = t1m + 2.0  ;

t2i  = 0          ;
t2m  = t2i + 2.5  ;
t2f  = t2m + 2.0  ;
tf   = 6          ;
dt   = 0.01       ;

r1i  = 0         ;
r2i  = 0         ;
r1m  = 4.5*pi/18 ;    
r2m  = 7.0*pi/18 ; 
r1f  = 4.0*pi/18 ;
r2f  = 4.0*pi/18 ;

%l1*cos(r1m) + l2*cos(r2m-r1m),
%l1*cos(r1f) + l2,

r2d = 180/pi;

open('REHROBNL'); SIMNL = sim('REHROBNL');

t       = SIMNL.get('t'     ); 
r1      = SIMNL.get('r1'    ); 
r2      = SIMNL.get('r2'    ); 
y1      = SIMNL.get('y1'    ); 
y2      = SIMNL.get('y2'    ); 
e1      = SIMNL.get('e1'    ); 
e2      = SIMNL.get('e2'    ); 
y1_dot  = SIMNL.get('y1_dot'); 
y2_dot  = SIMNL.get('y2_dot'); 
u1      = SIMNL.get('u1'    );
u2      = SIMNL.get('u2'    );

open('REHROBLN'); SIMLN = sim('REHROBLN');

tl      = SIMLN.get('t'     ); 
r1l     = SIMLN.get('r1'    ); 
r2l     = SIMLN.get('r2'    );
y1l     = SIMLN.get('y1'    );
y2l     = SIMLN.get('y2'    ); 
e1l     = SIMLN.get('e1'    ); 
e2l     = SIMLN.get('e2'    ); 
y1l_dot = SIMLN.get('y1_dot'); 
y2l_dot = SIMLN.get('y2_dot'); 
u1l     = SIMLN.get('u1'    );
u2l     = SIMLN.get('u2'    );


y_max = 100;
figure(8); subplot(321); plot(t,r1*r2d,'g-.',t,y1*r2d,'b-',tl,y1l*r2d,'r--','LineWidth',1.5); grid on;
ylabel('$y_1$ [deg]','Interpreter','latex'); 
legend('Reference','Nonlinear Simulation','Linear Simulation','Interpreter','latex','Location','NorthEast','Fontsize',8);
axis([0 tf 0 y_max]);
set(gca,'TickLabelInterpreter','latex');

figure(8); subplot(322); plot(t,r2*r2d,'g-.',t,y2*r2d,'b-',tl,y2l*r2d,'r--','LineWidth',1.5); grid on;
ylabel('$y_2$ [deg]','Interpreter','latex'); 
legend('Reference','Nonlinear Simulation','Linear Simulation','Interpreter','latex','Location','SouthEast');
axis([0 tf 0 y_max]);
set(gca,'TickLabelInterpreter','latex');

e_max = 2;
figure(8); subplot(323); plot(t,e1*r2d,'b-',tl,e1l*r2d,'r--','LineWidth',1.5); grid on;
ylabel('$e_1$ [deg]','Interpreter','latex'); 
axis([0 tf -e_max e_max]);
set(gca,'TickLabelInterpreter','latex');

figure(8); subplot(324); plot(t,e2*r2d,'b-',tl,e2l*r2d,'r--','LineWidth',1.5); grid on;
ylabel('$e_2$ [deg]','Interpreter','latex'); 
axis([0 tf -e_max e_max]);
set(gca,'TickLabelInterpreter','latex');

u_max = 60;
figure(8); subplot(325); plot(t,u1,'b-',tl,u1l,'r--','LineWidth',1.5); grid on;
ylabel('$u_1\,[\textrm{N}\cdot\textrm{m}]$','Interpreter','latex'); 
axis([0 tf 0 u_max]);
axis([0 tf -u_max/2 u_max]);
set(gca,'TickLabelInterpreter','latex');
xlabel('Time [s]','Interpreter','latex');

figure(8); subplot(326); plot(t,u2,'b-',tl,u2l,'r--','LineWidth',1.5); grid on;
ylabel('$u_2\,[\textrm{N}\cdot\textrm{m}]$','Interpreter','latex'); 
axis([0 tf -u_max/12 u_max/6]);
axis([0 tf -u_max/2 u_max]);
set(gca,'TickLabelInterpreter','latex');
xlabel('Time [s]','Interpreter','latex');



%% Animation
figure(10); plot([-0.2 0.2 0.2 0.6],[-l1-l2 -l1-l2 -l1*cos(r1f)-l2 -l1*cos(r1f)-l2],'k-','LineWidth',2); hold on;

axis(gca,'equal');
axis([-0.2 0.6 -1 0]); grid on;

N  = length(t);
sp = 1;       % Speed-up scaling

for i=1:N-1
    x_1 =  l1*sin(y1(i));
    y_1 = -l1*cos(y1(i));
    x_2 =  l1*sin(y1(i)) - l2*sin(y2(i)-y1(i));
    y_2 = -l1*cos(y1(i)) - l2*cos(y2(i)-y1(i));
    hip =  [0 0];
    kne =  [x_1 y_1];
    ank =  [x_2 y_2];
    thi =  line([hip(1) kne(1)], [hip(2) kne(2)],'Color','b','LineWidth',10);
    cal =  line([kne(1) ank(1)], [kne(2) ank(2)],'Color','b','LineWidth',10);
    tmr =  text(0.1,-0.05,"Timer: "+num2str(t(i),2));
    pause((t(i+1)-t(i))/sp); frames(i) = getframe(gcf);
    if i < N-1
       delete(thi);
       delete(cal);
       delete(tmr);
    end
end

obj           = VideoWriter('RehRob.avi');
obj.Quality   = 100                      ;
obj.FrameRate = 40                       ;
open(obj)                                ;

for i=1:N-1
    writeVideo(obj,frames(i));
end
obj.close();