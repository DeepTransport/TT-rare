% Log-likelihood + log-sigmoid for SIR, including 2 tempering parameters
% for DIRT
function [f] = ll_logsigmoid(betagamma, data, W, sigma_n, x0, tfine, ind_obs, Imax, gammaev, beta1, beta2)
d = size(W,1);
beta = betagamma(:,1:2:2*d-1)';
beta = reshape(beta, 1, d, []);
gamma = betagamma(:,2:2:2*d)';
gamma = reshape(gamma, 1, d, []);
x0 = repmat(x0, 1, size(beta,3));

[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tfine, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
x_obs = x(ind_obs, :);
x_obs = reshape(x_obs, [], 3*d, size(beta,3));
x_obs = x_obs(:, 2+(0:d-1)*3, :);
x_obs = reshape(x_obs, [], size(beta,3));

lF = -sum((x_obs - data).^2, 1)/(2*sigma_n^2);
lF = lF(:);

x = reshape(x, [], 3*d, size(beta,3));
x = x(:, 2+(d-1)*3, :);
x = max(x,[],1); % Check for the maximum number of infected in last compart
x = x(:);

Imax1 = Imax; %  + (Imax - 65)*log10(beta1)/3;
Imax2 = Imax; %  + (Imax - 65)*log10(beta2)/3;

f1 = log(1+exp((Imax1-x)*beta1*gammaev));
if (beta1==0)
    f1(:) = 0;
end
f1(isinf(f1)) = 0;
f2 = log(1+exp((Imax2-x)*beta2*gammaev));
f2(isinf(f2)) = +inf;
f = f1-f2;

f = f + lF*(beta2-beta1);
end
