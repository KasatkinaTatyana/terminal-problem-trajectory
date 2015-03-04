function F = calc_coefs(Y_0,Y_end,t_0,t_end)
A = [1 0 0 0;
    1 1 1 1;
    0 1 0 0;
    0 1 2 3];
B = [Y_0(1);
    Y_end(1);
    Y_0(2)*(t_end - t_0);
    Y_end(2)*(t_end - t_0)];
F = inv(A)*B;