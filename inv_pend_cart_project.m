%%%%%% ------------------------------------------------------------- %%%%%%
%%%%%% MEEG 829
%%%%%% Inverted Pendulum on a Cart Project
%%%%%% Karakasis Chrysostomos
%%%%%% ------------------------------------------------------------- %%%%%%
clear all;
close all;

%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 3
M = 25;                                             %Kg
m = 20;                                             %Kg
L = 9.81;                                           %m
g = 9.81;                                           %m/s^2
J = m*(L^2)/3;

b = ((M/m)+1)*(1+J/(m*(L^2)));                      %dimensionless variable 
c = (m^2)*(L^2)/(J*(m+M)+m*M*(L^2));                %dimensionless variable
d = c+1;                                            %dimensionless variable

%The linearization matrices A and B were computed by hand in the report
A=[0 1 0 0; 1 0 0 0; 0 0 0 1; -1 0 0 0 ];
B=[0; -1; 0; b];

%The desired eigenvalues for the state variable feedback control law
p1=-3;
p2=-2;
p3=-0.7-1i*0.2;
p4=-0.7+1i*0.2;
p=[p1 p2 p3 p4];

%Controllability Matrix
contr=[B A*B A*A*B A*A*A*B];
%Verification of the Controllability Test
rank(contr);

%Calculation of the gain matrix for the state feedback controller
K=place(A,B,p);
%%%%%% ------------------------------------------------------------- %%%%%%

%%%%%% ------------------------------------------------------------- %%%%%%

%Problem 5
syms x(t) [4 1] real

%Control Law
mu = -real(K)*x;

%Linear Model Closed-Loop System
x_dot_linear = (A-B*K)*x;

%Nonlinear Model Closed-Loop System
x1_dot_nl = x2;
x2_dot_nl = (-c*(x2^2)*sin(x1)*cos(x1)/(1+c*(sin(x1))^2)) + (sin(x1)/(1+c*(sin(x1))^2)) - (cos(x1)*mu/(1+c*(sin(x1))^2));
x3_dot_nl = x4;
x4_dot_nl = (d*(x2^2)*sin(x1)/(1+c*(sin(x1))^2)) - (sin(x1)*cos(x1)/(1+c*(sin(x1))^2)) + (b*mu/(1+c*(sin(x1))^2));
x_dot_nonlinear=[x1_dot_nl; x2_dot_nl; x3_dot_nl; x4_dot_nl];

%Differential Equations for Linear Model
odes = diff(x) == x_dot_linear;
%Differential Equations for Nonlinear Model
odes_nl= diff(x) == x_dot_nonlinear;

figure(1)
sgtitle('$\textbf{Figure (1) - Problem 5 - Simulation of Linear and Nonlinear Systems for Various Initial Conditions}$','Interpreter','latex')
subplot(4,1,1);
%The cart at rest but not at the origin and the bar leaning to one side,
%but also at rest
init_conds=[pi/11; 0; 0.25; 0];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Nonlinear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 1: Initial Conditions ($\frac{\pi}{11}$, 0, 0.25, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','northeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,2);
%The cart moving to one side with some velocity and the bar initially at
%rest and straight up
init_conds=[0; 0; 1.25; 1];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Nonlinear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 2: Initial Conditions (0, 0, 1.25, 1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,3);
%The cart moving to one side with some velocity and the bar leaning to one
%side not at rest
init_conds=[pi/11; pi/40; 0.2; 0.1];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Noninear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 3: Initial Conditions ($\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, 0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','northeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,4);
%The cart moving to one side with some negative velocity and the bar leaning to the                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                negative
%side not at rest
init_conds=[-pi/11; pi/40; 0.2; -0.1];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Nonlinear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 4: Initial Conditions ($-\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, -0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')
%%%%%% ------------------------------------------------------------- %%%%%%

%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 5 - Connection between Stability and the initial angle x1(0)

%The cart moving to one side with some negative velocity and the bar leaning to negative
%side not at rest
figure(2)
sgtitle('\textbf{Figure (2) - Problem 5 - Explore Stability of Nonlinear Model}','Interpreter','latex')
subplot(2,1,1);
init_conds=[0.658; 0; 0; 0];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Nonlinear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 1: Largest Value of Initial Angle for Stable Nonlinear Model - Initial Conditions: (0.658, 0, 0, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','northeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')
ylim([-15 15]);

subplot(2,1,2);
init_conds=[0.659; 0; 0; 0];
cond1 = x1(0) == init_conds(1);
cond2 = x2(0) == init_conds(2);
cond3 = x3(0) == init_conds(3);
cond4 = x4(0) == init_conds(4);
conds = [cond1; cond2; cond3; cond4];

%Solution of the ODE for the Linear Model
[x1Sol, x2Sol, x3Sol, x4Sol] = dsolve(odes,conds);

%Solution of the ODE for the Nonlinear Model
M=odeFunction(x_dot_nonlinear,x);
[t,sol] = ode45(M,[0 10],init_conds);

%Plotting of the Responses for both Close-Loop Systems
fplot(x1Sol,[0 10],'Color',[0 0.4470 0.7410])
hold on
fplot(x2Sol,[0 10],'Color',[0.8500 0.3250 0.0980])
hold on
fplot(x3Sol,[0 10],'Color',[0.9290 0.6940 0.1250])
hold on
fplot(x4Sol,[0 10],'Color',[0.4940 0.1840 0.5560])
hold on;
plot(t,sol(:,1),'m--','Color',[0 0.4470 0.7410])
hold on;
plot(t,sol(:,2),'y--','Color',[0.8500 0.3250 0.0980])
hold on;
plot(t,sol(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t,sol(:,4),'--','Color',[0.4940 0.1840 0.5560])
grid on
title('Subplot 2: Example of Unstable Behavior of the Nonlinear Model - Initial Conditions: (0.659, 0, 0, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$\phi_{L}$','$\dot{\phi}_{L}$','$\bar{s}_{L}$','$\dot{\bar{s}}_{L}$','$\phi_{NL}$','$\dot{\phi}_{NL}$','$\bar{s}_{NL}$','$\dot{\bar{s}}_{NL}$','Location','southwest','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')
ylim([-25 25]);
%%%%%% ------------------------------------------------------------- %%%%%%

%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 6

%We derive the C matrix from the analysis in the report
C=[1 0 0 0; 0 0 1 0];
%Observability Matrix
O=[C; C*A; C*A*A];
rank(O);

%The desired eigenvalues for the observer matrix A-LC
% p1=-3;
% p2=-2;
% p3=-0.7-1i*0.2;
% p4=-0.7+1i*0.2;
pobs=3*[p1 p2 p3 p4];

%Calculation of the gain matrix for the state feedback controller
Ltran=place(A',C',pobs);
L=Ltran';
%%%%%% ------------------------------------------------------------- %%%%%%

%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 7
syms x_full(t) [8 1] real
x_full_14=[x_full1;x_full2;x_full3;x_full4];
x_full_58=[x_full5;x_full6;x_full7;x_full8];

%Linear Model
x_dot_full_14 = A*x_full_14 -B*K*x_full_58;
x_dot_full_58 = (A-B*K-L*C)*x_full_58+L*C*x_full_14;
x_dot_full=[x_dot_full_14;x_dot_full_58];

% Nonlinear Model
%Appliead Control Law
mu = -real(K)*x_full_58;
x_dot_full_nl_1 = x_full2;
x_dot_full_nl_2 = (-c*(x_full2^2)*sin(x_full1)*cos(x_full1)/(1+c*(sin(x_full1))^2)) + (sin(x_full1)/(1+c*(sin(x_full1))^2)) - (cos(x_full1)*mu/(1+c*(sin(x_full1))^2));
x_dot_full_nl_3 = x_full4;
x_dot_full_nl_4 = (d*(x_full2^2)*sin(x_full1)/(1+c*(sin(x_full1))^2)) - (sin(x_full1)*cos(x_full1)/(1+c*(sin(x_full1))^2)) + (b*mu/(1+c*(sin(x_full1))^2));
x_dot_full_nl_58 = (A-B*K-L*C)*x_full_58+L*C*x_full_14;
x_dot_full_nonlinear=[x_dot_full_nl_1; x_dot_full_nl_2; x_dot_full_nl_3; x_dot_full_nl_4; x_dot_full_nl_58];

figure(3)
sgtitle('$\textbf{Figure (3) - Problem 7 - Simulation of Dynamic Controller for Linear and Nonlinear System Models}$','Interpreter','latex')

subplot(4,1,1);
%The cart at rest but not at the origin and the bar leaning to one side,
%but also at rest
init_conds=[pi/11; 0; 0.25; 0];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 15],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 15],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
grid on
title('Subplot 1: Initial Conditions ($\frac{\pi}{11}$, 0, 0.25, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,2);
%The cart moving to one side with some velocity and the bar initially at
%rest and straight up
init_conds=[0; 0; 1.25; 1];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 15],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 15],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
grid on
title('Subplot 2: Initial Conditions (0, 0, 1.25, 1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,3);
%The cart moving to one side with some velocity and the bar leaning to one
%side not at rest
init_conds=[pi/11; pi/40; 0.2; 0.1];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 15],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 15],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
grid on
title('Subplot 3: Initial Conditions ($\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, 0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,4);
%The cart moving to one side with some negative velocity and the bar leaning to the                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                negative
%side not at rest
init_conds=[-pi/11; pi/40; 0.2; -0.1];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 10],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 10],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
grid on
title('Subplot 4: Initial Conditions ($-\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, -0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

%%%%%% ------------------------------------------------------------- %%%%%%

%Problem 7 - Connection between Stability of Dynamic Controller and the initial angle x1(0)
figure(4)
sgtitle('\textbf{Figure (4) - Problem 7 - Explore Stability of Nonlinear Model for the Dynamic Controller}','Interpreter','latex')

subplot(2,1,1);
%I.C.s with the Largest Value of x1(0) with Stable Behavior
init_conds=[0.340; 0; 0; 0];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 20],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 20],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
ylim([-15 15]);
grid on
title('Subplot 1: Initial Conditions (0.340, 0, 0, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(2,1,2);
%I.C.s with the first Value of x1(0) with Unstable Behavior
init_conds=[0.341; 0; 0; 0];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 20],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 20],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
ylim([-15 15]);
grid on
title('Subplot 2: Initial Conditions (0.341, 0, 0, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')
%%%%%% ------------------------------------------------------------- %%%%%%
    
%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 8
syms x_full(t) [8 1] real
x_full_14=[x_full1;x_full2;x_full3;x_full4];
x_full_58=[x_full5;x_full6;x_full7;x_full8];

%Linear Model - Same as in Problem 7
x_dot_full_14 = A*x_full_14 -B*K*x_full_58;
%Linear Observer 
x_dot_full_58 = (A-B*K-L*C)*x_full_58+L*C*x_full_14;
x_dot_full=[x_dot_full_14;x_dot_full_58];

% Nonlinear Model - Same as in Problem 7
mu = -real(K)*x_full_58;
x_dot_full_nl_1 = x_full2;
x_dot_full_nl_2 = (-c*(x_full2^2)*sin(x_full1)*cos(x_full1)/(1+c*(sin(x_full1))^2)) + (sin(x_full1)/(1+c*(sin(x_full1))^2)) - (cos(x_full1)*mu/(1+c*(sin(x_full1))^2));
x_dot_full_nl_3 = x_full4;
x_dot_full_nl_4 = (d*(x_full2^2)*sin(x_full1)/(1+c*(sin(x_full1))^2)) - (sin(x_full1)*cos(x_full1)/(1+c*(sin(x_full1))^2)) + (b*mu/(1+c*(sin(x_full1))^2));
%Nonlinear Observer
x_dot_full_nl_5 = x_full6;
x_dot_full_nl_6 = (-c*(x_full6^2)*sin(x_full5)*cos(x_full5)/(1+c*(sin(x_full5))^2)) + (sin(x_full5)/(1+c*(sin(x_full5))^2)) - (cos(x_full5)*mu/(1+c*(sin(x_full5))^2));
x_dot_full_nl_7 = x_full8;
x_dot_full_nl_8 = (d*(x_full6^2)*sin(x_full5)/(1+c*(sin(x_full5))^2)) - (sin(x_full5)*cos(x_full5)/(1+c*(sin(x_full5))^2)) + (b*mu/(1+c*(sin(x_full5))^2));
x_dot_full_nl_58 = [x_dot_full_nl_5; x_dot_full_nl_6; x_dot_full_nl_7; x_dot_full_nl_8] +L*C*(x_full_14-x_full_58);
x_dot_full_nonlinear=[x_dot_full_nl_1; x_dot_full_nl_2; x_dot_full_nl_3; x_dot_full_nl_4; x_dot_full_nl_58];

figure(5)
sgtitle('$\textbf{Figure (5) - Problem 8 - Simulation of Dynamic Controller With Nonlinear Observer}$','Interpreter','latex')

subplot(4,1,1);
%The cart at rest but not at the origin and the bar leaning to one side,
%but also at rest
init_conds=[pi/11; 0; 0.25; 0];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 10],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
% [t_nonlinear,sol_nonlinear] = ode45(M,[0 5],init_conds);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 10],init_conds);

%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
% ylim([-15 15]);
grid on
title('Subplot 1: Initial Conditions ($\frac{\pi}{11}$, 0, 0.25, 0)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,2);
%The cart moving to one side with some velocity and the bar initially at
%rest and straight up
init_conds=[0; 0; 1.25; 1];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 10],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 10],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
% ylim([-15 15]);
grid on
title('Subplot 2: Initial Conditions (0, 0, 1.25, 1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,3);
%The cart moving to one side with some velocity and the bar leaning to one
%side not at rest
init_conds=[pi/11; pi/40; 0.2; 0.1];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 10],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 10],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
% ylim([-15 15]);
grid on
title('Subplot 3: Initial Conditions ($\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, 0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

subplot(4,1,4);
%The cart moving to one side with some negative velocity and the bar leaning to the                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                negative
%side not at rest
init_conds=[-pi/11; pi/55; 0.5; -0.25];
%Always take the initial conditions of the observer to be [0,0,0,0]
init_conds=[init_conds;0;0;0;0];
%Solution of Differential Equations for Linear Model
M=odeFunction(x_dot_full,x_full);
[t_linear,sol_linear] = ode45(M,[0 10],init_conds);
%Solution of Differential Equations for Nonlinear Model 
M=odeFunction(x_dot_full_nonlinear,x_full);
[t_nonlinear,sol_nonlinear] = ode45(M,[0 10],init_conds);
%Plotting ouputs y1(t) and y2(t)
%y1(t)=phi=x1=x_full1(t)
plot(t_linear,sol_linear(:,1),'','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_linear,sol_linear(:,3),'','Color',[0.9290 0.6940 0.1250])
hold on;
plot(t_nonlinear,sol_nonlinear(:,1),'--','Color',[0 0.4470 0.7410])
hold on;
%y2(t)=s bar=x3=x_full3(t)
plot(t_nonlinear,sol_nonlinear(:,3),'--','Color',[0.9290 0.6940 0.1250])
hold on;
% ylim([-15 15]);
grid on
title('Subplot 4: Initial Conditions ($-\frac{\pi}{11}$, $\frac{\pi}{40}$, 0.2, -0.1)','Interpreter','latex')
xlabel('Time')
lgd=legend('$y_{1}(t)=\phi_{L}$','$y_{2}(t)=\bar{s}_{L}$','$y_{1}(t)=\phi_{NL}$','$y_{2}(t)=\bar{s}_{NL}$','Location','southeast','Interpreter','latex','NumColumns',2);
title(lgd,'$\textbf{L}$ : Linear Model - $\textbf{NL}$ : Nonlinear Model','Interpreter','latex')

%%%%%% ------------------------------------------------------------- %%%%%%
%Problem 8 - Linearization about the Origin
A_augm=[0 1 0 0 0 0 0 0; 
        1 0 0 0 K(1) K(2) K(3) K(4);
        0 0 0 1 0 0 0 0;
        -1 0 0 0 -b*K(1) -b*K(2) -b*K(3) -b*K(4);
        L(1,1) 0 L(1,2) 0 -L(1,1) 1 -L(1,2) 0;
        L(2,1) 0 L(2,2) 0 1+K(1)-L(2,1) K(2) K(3)-L(2,2) K(4);
        L(3,1) 0 L(3,2) 0 -L(3,1) 0 -L(3,2) 1;
        L(4,1) 0 L(4,2) 0 -1-b*K(1)-L(4,1) -b*K(2) -b*K(3)-L(4,2) -b*K(4)]
eig(A_augm)    
%%%%%% ------------------------------------------------------------- %%%%%%


