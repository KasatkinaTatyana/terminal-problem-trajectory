function F = trajectory_synthesis(Y_0, Y_end, t_0, t_end, dt)

c_t = calc_coefs(Y_0,Y_end,t_0,t_end);
c0 = c_t(1);
c1 = c_t(2);
c2 = c_t(3);
c3 = c_t(4);

t = t_0:dt:t_end;

Psi = c0 + c1/(t_end - t_0)*(t - t_0) + c2/(t_end - t_0)^2*(t - t_0).^2 + c3/(t_end - t_0)^3*(t - t_0).^3;
dPsi = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t - t_0) + 3*c3/(t_end - t_0)^3*(t - t_0).^2;
ddPsi =  2*c2/(t_end - t_0)^2 + 6*c3/(t_end - t_0)^3*(t - t_0);

F = [Psi; dPsi; ddPsi];