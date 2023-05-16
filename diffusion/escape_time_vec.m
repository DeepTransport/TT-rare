function [t,x,len] = escape_time_vec(a, u)
% A technical function to compute escape time in a diffusion flux given by
% diffusivity coefficient a and pressure head u

n = round(sqrt(size(u,2)));
N = size(u, 1);

t = zeros(N,1);
x = repmat([0, 0.5], N, 1);
active_ind = (1:N)';
x_active = x;
t_active = t;
stepcnt = 0;
len = zeros(N,1);
while (~isempty(active_ind))
    F = flux_interpolate_vec(x_active, a, u);  % size N x 2
    dt = 1./((n-1)*2*sqrt(sum(F.^2,2)));
    x_active = x_active + dt.*F;
    t_active = t_active + dt;
    x(active_ind, :) = x_active;
    t(active_ind) = t_active;    
    len(active_ind) = len(active_ind) + norm(dt.*F);
    not_active = (x_active(:,1)>=1); %  | (t_active>1e2);
    active_ind(not_active) = [];
    x_active(not_active, :) = [];
    t_active(not_active) = [];
    a(not_active, :) = [];
    u(not_active, :) = [];
    stepcnt = stepcnt + 1;
    if (stepcnt>3e4)
        break;
    end
end
end
