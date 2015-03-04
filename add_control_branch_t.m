function add_control_branch_t(row)
% Функция добавляет в массив control_arr, описывающий фазовую кривую
% системы новую строку row
global control_arr

N = length(control_arr(:,1));

t_0_add = row(1);
t_end_add = row(2);

for i=1:N
    t_0_i = control_arr(i,1);
    t_end_i = control_arr(i,2);
    
    if  (abs(t_0_add - t_0_i) < 1e-8)&&(abs(t_end_add - t_end_i) < 1e-8)
        control_arr(i,:) = row;
        break;
    end
    if  ( (abs(t_0_add - t_0_i) < 1e-8)&&(t_end_add < t_end_i) )
        control_arr(i,:) = row;
        control_arr = [control_arr; [t_end_add t_end_i control_arr(i,3:10)]];
        break;
    end
    if  ( (t_0_add > t_0_i)&&(abs(t_end_add - t_end_i) < 1e-8) )
        control_arr(i,:) = [t_0_i t_0_add control_arr(i,3:10)];
        control_arr = [control_arr; row];
        break;
    end
    if ( (t_0_add > t_0_i)&&(t_end_add < t_end_i) )
        control_arr(i,:) = [t_0_i t_0_add control_arr(i,3:10)];
        control_arr = [control_arr; row];
        control_arr = [control_arr; [t_end_add t_end_i control_arr(i,3:10)]];
        break;
    end
end