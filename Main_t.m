clear all
close all
clc
%% Ограничения
%  constr1 < Psi(y) < constr2
%  dy_constr1 < dPsi(y) / dy < dy_constr2
dy_constr1 = -1.9;
dy_constr2 = 1;

y0 = pi / 4;
dy0 = 0.1;
yend = 3*pi/4;
dyend = 0.5;
% u0 = 0;
% uend = 0;
% ddy0 = u0 + sin(y0);
% ddyend = uend + sin(yend);

Y_0 =   [y0    dy0  ];
Y_end = [yend  dyend];

t_0 = 0;
t_end = 2;

c_t = calc_coefs(Y_0,Y_end,t_0,t_end);
c0 = c_t(1);
c1 = c_t(2);
c2 = c_t(3);
c3 = c_t(4);

%% Массив, в котором будут содержаться коэффициенты для управлений
global control_arr
control_arr = [t_0 t_end c0 c1 c2 c3 0 0 t_0 t_end];
% 1 : t_0 - начальная точка кривой;
% 2 : t_end - конечная точка кривой;
% 3-6 : коэффициенты;
% 7 : d - параметр;
% 8 :  тип кривой
% 9 : t_0_c - это значение, которое вычитается из t 
% 10 : t_end_c  

% типы будут следующие:
% 0 - кривая вида Psi_0(t) = c0 + 1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - 
% t_0_c)^2*(t - t_0_c)^2 + c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3;
% 1 - кривая вида 
% 2 - кривая вида Psi_2(t) = Psi_0(t) + d*t_norm.^2.*(3 - 2*t_norm), 
% где y_norm = (t - t_0) / (t_end - t_0);

%%
figure(1);
hold on; grid on;
title('Psi(t)');
xlabel('t');
ylabel('y');

N = 1000;
dt = t_end / N;
t = t_0:dt:t_end;

Psi = c0 + c1/(t_end - t_0)*(t - t_0) + c2/(t_end - t_0)^2*(t - t_0).^2 + c3/(t_end - t_0)^3*(t - t_0).^3;  % y(t)
dPsi = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t - t_0) + 3*c3/(t_end - t_0)^3*(t - t_0).^2;               % \dot{y}  
ddPsi = 2*c2/(t_end - t_0)^2 + 6*c3/(t_end - t_0)^3*(t - t_0);
u = ddPsi - sin(Psi);

plot(t,Psi);

figure(2);
hold on; grid on;
xlabel('t');
ylabel('dy / dt');
title('dPsi');
plot(t,dPsi,'b');

figure(3);
hold on; grid on;
xlabel('t, c');
ylabel('u(t)');
title('Управление u(t)');
plot(t,u);

%% Убираю выход за ограничение dPsi / dy < dy_constr2
[data1, data2, par, coefs] = up_dPsi_t(Y_0, Y_end, t_0, t_end, dt, dy_constr2);
t_1 = data1(1,:);
Psi_1 = data1(2,:);
dPsi_1 = data1(3,:);
ddPsi_1 = data1(4,:);
u_1 = ddPsi_1 - sin(Psi_1);

t_2 = data2(1,:);
Psi_2 = data2(2,:);
dPsi_2 = data2(3,:);
ddPsi_2 = data2(4,:);
u_2 = ddPsi_2 - sin(Psi_2);

set(0,'CurrentFigure',1);
plot(t_1, Psi_1, 'g','LineWidth',2);
plot(t_2, Psi_2, 'm','LineWidth',2);

set(0,'CurrentFigure',2);
plot(t_1, dPsi_1, 'g','LineWidth',2);
plot(t_2, dPsi_2, 'm','LineWidth',2);

set(0,'CurrentFigure',3);
plot(t_1, u_1, 'LineStyle','--','Color','b','LineWidth', 2); % color g
plot(t_2, u_2, 'LineStyle','--','Color','b','LineWidth', 2); % color m

row = [t_1(1) t_1(end) coefs(1,:) par(1) 2 par(2) par(3)];

add_control_branch_t(row);
 
row = [t_2(1) t_2(end) coefs(2,:) 0 0 t_2(1) t_2(end)];
 
add_control_branch_t(row);
%% интегрирование системы по времени
[time, traj, U] = t_modeling(Y_0,t_0,t_end);
% [time, traj, U] = t_modeling([Psi_1(1) dPsi_1(1)],t_1(1),t_1(end));

% set(0,'CurrentFigure',1);
% plot(traj(:,1), traj(:,2),'y');

figure(4);
subplot(2,1,1);
plot(time,traj(:,1));
xlabel('t, c'); ylabel('y(t)');

subplot(2,1,2);
plot(time,traj(:,2));
xlabel('t, c'); ylabel('dy / dt');

figure(5);
hold on; grid on;
xlabel('t, c');
ylabel('u(t)');
title('Стабилизирующее правление u(t)');
plot(U(:,1), U(:,2));
