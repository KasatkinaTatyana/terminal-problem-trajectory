function dy = control_system(t,y)
dy = zeros(2,1);

global control_arr
global count

global Mas_u

flag = false;

% Для каждого типа кривой (8 координата в строках массива control_arr)
% управление строится по своему

for i=1:count
    t_0_i = control_arr(i,1);
    t_end_i = control_arr(i,2);
    
    if ((t < t_end_i)&&(t > t_0_i)) || (abs(t - t_0_i) < 1e-6)
        coefs = control_arr(i,3:6);
        c0 = coefs(1);
        c1 = coefs(2);
        c2 = coefs(3);
        c3 = coefs(4);
        d = control_arr(i,7);
        t_0_c = control_arr(i,9);
        t_end_c = control_arr(i,10);
                     
        if (abs(control_arr(i,8) - 2) < 1e-6)           
            
            Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 + ...
                c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3 +...
                d*((t - t_0_i)/(t_end_i - t_0_i))^2*(3 - 2*((t - t_0_i)/(t_end_i - t_0_i) ));
            
            dPsi_t = c1/(t_end_c - t_0_c) + 2*c2/(t_end_c - t_0_c)^2*(t - t_0_c) + 3*c3/(t_end_c - t_0_c)^3*(t - t_0_c)^2 + ...
                d*(6*((t - t_0_i)/(t_end_i - t_0_i)^2 - 6*(t - t_0_i)^2/(t_end_i - t_0_i)^3));
            
            ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c) + ...
                   6*d*(1/(t_end_i - t_0_i)^2 - 2*(t - t_0_i) / (t_end_i - t_0_i)^3);
               
            r1 = -0.001; r2 = -0.001;
            c1_st = -(r1 + r2);
            c0_st = r1*r2;

            u = ddPsi_t - sin(y(1)) - c1_st*(y(2) - dPsi_t) - c0_st*(y(1) - Psi_t);
        else         
%             Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 + ...
%                 c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3;
            ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c);
            u = ddPsi_t - sin(y(1));
        end 
        flag = true;
        break;
    end
    
end

% последняя точка ?
if (flag == false)    
    coefs = control_arr(count,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);

    t_0_c = control_arr(count,9);
    t_end_c = control_arr(count,10);
    
    ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c);
    
    u = ddPsi_t - sin(y(1));
end

Mas_u = [Mas_u; t u];

dy(1) = y(2);
dy(2) = sin(y(1)) + u;