function [F1, F2, params, coefs] = up_dPsi_t(Y_0, Y_end, t_0, t_end, dt, dy_constr2)
% Функция генерирует кривую, отрезок, на котором она строится, и ее
% производную таким образом, чтобы выполнялось ограничение dPsi(t) / dt < dy_constr2

% F1 соответсвует первой кривой (левой),
% F2 - правой, они имеют следующую структуру
% массив t - массив функции Psi(t) - массив dPsi / dt

% params и coefs нужны для того, чтобы затем однозначно восстановить кривую
% и построить нужное управление

c_t = calc_coefs(Y_0,Y_end,t_0,t_end);
c0 = c_t(1);
c1 = c_t(2);
c2 = c_t(3);
c3 = c_t(4);

t = t_0:dt:t_end;

Psi = c0 + c1/(t_end - t_0)*(t - t_0) + c2/(t_end - t_0)^2*(t - t_0).^2 + c3/(t_end - t_0)^3*(t - t_0).^3;
dPsi = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t - t_0) + 3*c3/(t_end - t_0)^3*(t - t_0).^2;

N_shift = 50;

N = (t_end - t_0)/dt + 1;

flag = 0;
% Ищем границы интервала, на котором нарушается ограничение
for i=1:N
    % ------- -- dPsi / dy < constr -----------------------------------
    if ((flag == 0)&&(dPsi(i) >= dy_constr2))
    % -----------------------------------------------------------------
        if ((i - N_shift) < 1)
            t_left = t_0;
            Y_left = Y_0;
        else
            t_left = t(i - N_shift);
            y_left = Psi(i - N_shift);
            dy_left = dPsi(i - N_shift);
        end
        flag = 1;
    end
    % --------- dPsi / dy < constr ------------------------------------
    if ((flag==1)&&(dPsi(i) < dy_constr2))
    % -----------------------------------------------------------------
        if ((i + N_shift) > N + 1)
            t_right = t_end;
            Y_right = Y_end;
        else
            t_right = t(i + N_shift);
            y_right = Psi(i + N_shift);
            dy_right = dPsi(i + N_shift);
        end
        break;
    end
end

t_1 = t_left:dt:t_right;
t_norm = (t_1 - t_left) / (t_right - t_left);     
h_d = -0.0001;
for d = 0:h_d:(-0.5)
    Psi_1 = c0 + c1/(t_end - t_0)*(t_1 - t_0) + c2/(t_end - t_0)^2*(t_1 - t_0).^2 + c3/(t_end - t_0)^3*(t_1 - t_0).^3 ...
            + d*(t_norm.^2).*(3 - 2*t_norm);
    dPsi_1 = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t_1 - t_0) + 3*c3/(t_end - t_0)^3*(t_1 - t_0).^2 ...
        + d*(6*(t_1 - t_left)/(t_right - t_left)^2 - 6*(t_1 - t_left).^2/(t_right - t_left)^3 );
            
    if (max(dPsi_1) < dy_constr2)
        break;
    end
end
ddPsi_1 = 2*c2/(t_end - t_0)^2 + 6*c3/(t_end - t_0)^3*(t_1 - t_0) + ...
          6*d*(1/(t_right - t_left)^2 - 2*(t_1 - t_left) / (t_right - t_left)^3);

Y_1_end = [Psi_1(end) dPsi_1(end)];

replace_part = trajectory_synthesis(Y_1_end, Y_end, t_right, t_end, dt);
coefs_2 = calc_coefs(Y_1_end, Y_end, t_right, t_end);
t_2 = t_right:dt:t_end;

Psi_2 = replace_part(1,:);
dPsi_2 = replace_part(2,:);
ddPsi_2 = 2*coefs_2(3)/(t_2(end)- t_2(1))^2 + 6*coefs_2(4)/(t_2(end) - t_2(1))^3*(t_2 - t_2(1));

F1 = [t_1; Psi_1; dPsi_1; ddPsi_1];
F2 = [t_2; Psi_2; dPsi_2; ddPsi_2];

params = [d; t_0; t_end; 0];

coefs = [c0 c1 c2 c3; coefs_2'; 0 0 0 0;0 0 0 0];