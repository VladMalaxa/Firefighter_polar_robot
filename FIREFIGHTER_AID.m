clear all, close all, clc

%q
q1_0 = 0;
q2_0 = 0;
q1 = sym('q1',[length(q0_1) 1]);        % Create the vector q with length(q0) symbolic parameters   
q2 = sym('q2',[length(q0_2) 1]);        % Create the vector q with length(q0) symbolic parameters   
dq_1 = sym('dq1',[length(q0_1) 1]);      % Create the vector dq with length(q0) symbolic parameters  
dq_2 = sym('dq2',[length(q0_2) 1]);      % Create the vector dq with length(q0) symbolic parameters  
ddq_1 = sym('ddq1',[length(q0_1) 1]);
ddq_2 = sym('ddq2',[length(q0_1) 1]);
syms t                              % Symbolic expression for time

%variables
l1 = 1;
l2 = 1;
l_c1 = 1;
m1 = 1;
m2 = 1;
alpha = 1;
g = 1;
q = [q1;q2];
dq = [dq_1,dq_2];
ddq = [ddq_1 ; ddq_2];
%M(q)
M = [l1 + l2 + (l_c1^2)*m1 + m2*(q2); 0;
    0; m2 ];

%c(q,dq)
c = [2*dq_2*dq_2*m2*q2 ;
    (-1)*(dq_1^2)*m2*q2];
%G(q)
G = [g*cos(alpha-q1)*(l_c1*m1 + m2*q2);
    (-1)*m2*sin(alpha - q1)];
%p_e(q)
p_e = [q2*cos(q1);
        q2*sin(q1)]; 
%Jacobian
%M*ddq + c + G = eye(2)*u;