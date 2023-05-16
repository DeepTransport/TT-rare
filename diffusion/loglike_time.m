function [L,t,xexit,len] = loglike_time(theta, beta0, beta, tmin, phil, bound, W2, Mass_summed, Q_obs, sigma_n, gamma)
% Technical function to compute log-biasing function of escape time
% exceeding tmin

[Q_cur, a, u] = diffusion_QoI(theta, phil, bound, W2, Mass_summed);
a = a.';
u = u.';  % size N x n^2

[t,xexit,len] = escape_time_vec(a, u);

f1 = log(1+exp((t-tmin)*beta0/gamma));
if (beta0==0)
    f1(:) = 0;
end
f1(isinf(f1)) = 0;
f2 = log(1+exp((t-tmin)*beta/gamma));
if (gamma==0)
    f2 = -log(double(t<tmin));
end
f2(isinf(f2)) = +inf;

L = f1-f2   -sum((Q_cur - Q_obs).^2, 2)*(beta-beta0)/(2*sigma_n);
end
